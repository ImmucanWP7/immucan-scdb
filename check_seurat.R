#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
seurat_obj = args[1] #path of seurat object
batch_var = args[2] #batch variable if known
verbose = FALSE
tidy_metadata_path <- "/home/jordi_camps/IMMUcan/tidy_metadata.xlsx"

dir <- getwd()
setwd(dir)
print(dir)
if (!dir.exists("temp")) {dir.create("temp")}
if (!dir.exists("temp/QC")) {dir.create("temp/QC")}
if (!dir.exists("out")) {dir.create("out")}
if (!dir.exists("out/plots")) {dir.create("out/plots")}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(dplyr)
  library(readxl)
  library(DescTools)
  library(tidyr)
  library(tibble)
  library(jsonlite)
  library(harmony)
  library(kBET)
  library(scater)
})

print("STEP 1: CHECKING SEURAT OBJECT")

if (is.na(seurat_obj)) {
  seurat_obj <- normalizePath(list.files(pattern = ".rds$"))
  if (length(seurat_obj) != 1) {
    stop("Specify seurat object in arguments")
  }
}
seurat <- readRDS(seurat_obj)

# Load data.json or create standard
if (file.exists("out/data.json")) {
  data <- fromJSON("out/data.json")
} else {
  data <- list()
  data$object_path = seurat_obj
  data$batch = batch_var
  data$QC_feature_min = 250 #Minimal features threshold
  data$QC_mt_max = 20 #Maximum mitochondrial content threshold
  data$pca_dims = 30 #Amount of PCA dimensions to use
  data$features_var = 2000 #Amount of variable features to select
  data$nSample = 10000
  data$cluster_resolution = seq(from = 0.4, to = 4, by = 0.1)
  data$malignant = TRUE
  data$normal_cells = NA
}

print(paste0("nCell = ", ncol(seurat), " / ", "nGene = ", nrow(seurat)))

if (sum(colnames(seurat) == rownames(seurat@meta.data)) == ncol(seurat)) {
  print("Cell IDs linked correctly")
} else {
  stop("Cell IDs linked uncorrectly")
}

gapdh <- grepl("GAPDH|Gapdh", rownames(seurat))
if (sum(grepl("\\.", seurat[["RNA"]]@counts[gapdh, ])) == 0) {
  data$norm <- FALSE
  print("Raw counts supplied")
} else if (any(seurat[["RNA"]]@counts[gapdh, ] > 100)) {
  data$norm <- FALSE
  print("Normalized counts supplied. Be careful with further interpretation")
} else {
  data$norm <- TRUE
  print("Logcounts supplied, no normalization will be done. Be careful with further interpretation")
}

seurat[["percent_mt"]] <- PercentageFeatureSet(seurat, pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-")
seurat[["RNA"]]@counts[1:5,1:5]
dplyr::glimpse(seurat@meta.data)

# Clean metadata columns
names <- tolower(colnames(seurat@meta.data))
names <- gsub("\\.", "_", names)
meta_cols <- read_excel(tidy_metadata_path)
if (any(names %in% meta_cols$col_names)) {
  change_cols <- colnames(seurat@meta.data)[colnames(seurat@meta.data) %in% meta_cols$col_names]
  for (i in change_cols) {
    hit_1 <- grepl(i, names)
    hit_2 <- meta_cols$col_names %in% i
    print(paste0("changing ", i, " to ", meta_cols$general[hit_2]))
    colnames(seurat@meta.data)[hit_1] <- meta_cols$general[hit_2]
  }
} else {
  print("No meta.data columns to tidy up")
}
data$annotation = c("seurat_clusters","annotation_CHETAH","annotation_major","annotation_immune","annotation_minor", colnames(seurat@meta.data)[grepl("authors_annotation|Authors_annotation", colnames(seurat@meta.data))])
saveRDS(seurat, seurat_obj)

# Batch
print("STEP 2a: ESTIMATING BATCH VARIABLES")

## Remove bad quality cells
seurat <- CreateSeuratObject(counts = seurat[["RNA"]]@counts, meta.data = seurat@meta.data, min.cells = 10, min.features = 200)
seurat <- subset(seurat, subset = nFeature_RNA > data$QC_feature_min & percent_mt < data$QC_mt_max)

## Subsample datasets larger than 20k cells
if (ncol(seurat) > 50000) {
  subset_size <- 0.1 #subsample to 10% of the data
  subset_id <- sample.int(n = ncol(seurat), size = floor(subset_size * ncol(seurat)), replace=FALSE)
  seurat <- seurat[, subset_id]
}

## Select potential batch columns from meta.data
seurat@meta.data <- seurat@meta.data[, sapply(sapply(seurat@meta.data, unique), length) != 1, drop = FALSE] #Remove all columns that have only one variable
#seurat@meta.data <- seurat@meta.data[, sapply(sapply(seurat@meta.data, unique), length) != nrow(seurat@meta.data), drop = FALSE] #Remove columns with only unique values
meta <- seurat@meta.data[, sapply(seurat@meta.data, class) %in% c("character", "factor")] #Select all columns that are factor or character
meta <- meta[,!grepl("Cluster|cluster|author|Author|Annotation|annotation|Cell_type|cell_type|cell|Cell|barcode|Barcode", colnames(meta))]
if (!"metadata" %in% names(data)) {data$metadata = colnames(meta)} #save metadata columns to data.json
if (length(data$batch) > 1) {
  batch <- data$batch
} else if (is.na(data$batch)) {
  temp <- meta[, apply(meta, 2, function(x) !any(is.na(x))), drop = FALSE] #Remove all columns with NAs
  batch <- colnames(temp)
} else {
  batch <- data$batch
}
print("Possible batches:")
print(paste0(batch))

## Create nearest neighbour graph
if (data$norm == FALSE) {
  seurat <- Seurat::NormalizeData(seurat, verbose = verbose)
} else {seurat[["RNA"]]@data <- seurat[["RNA"]]@counts}

seurat <- seurat %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = data$features_var, verbose = verbose) %>% 
  ScaleData(verbose = verbose) %>% 
  RunPCA(pc.genes = seurat@var.genes, npcs = data$pca_dims+20, verbose = verbose) %>%
  RunUMAP(dims = 1:data$pca_dims, a = .5, b = 1.2, verbose = verbose) %>%
  FindNeighbors(dims = 1:2, k.param = 30, reduction = "umap", verbose = verbose)

p_lbw <- ElbowPlot(seurat, ndims = data$pca_dims+20) + geom_vline(xintercept = data$pca_dims, color = "red") + ylab("STDEV PCA") + theme(axis.title.x = element_blank())

## Compute variance explained
#print("Compute explained variance")
#temp <- seurat@meta.data[, apply(seurat@meta.data, 2, function(x) !any(is.na(x))), drop = FALSE] #Remove all columns with NAs
#sce <- SingleCellExperiment(assays = list(logcounts = seurat[["RNA"]]@data), colData = temp)
#vars_exp <- getVarianceExplained(x = sce, variables = colnames(temp))
#p_var <- plotExplanatoryVariables(vars_exp)
#ggsave(plot = p_var, filename = paste0("temp/QC/Batch_variance_explained.png"))

#print("Compute kBET score")
## kBET
#kbet <- list()
#mtx <- t(as.matrix(seurat[["RNA"]]@data))
#for (b in batch) {
#  btch <- seurat@meta.data[, b]
#  kbet.estimate <- kBET(df = mtx, batch = btch, plot = FALSE)
#  kbet[[b]] <- 1 - kbet.estimate$stats$kBET.observed
#}
##Plot kbet scores
#kbet <- do.call(cbind, kbet)
#p_kbet <- kbet %>%
#  as.data.frame() %>%
#  gather("var", "Acceptance_rate") %>%
#  ggplot(aes(x = var, y = Acceptance_rate)) +
#  geom_boxplot()
#ggsave(plot = p_kbet, filename = "temp/QC/batch_kBET.png")

## Compute the percentage of batch in cell neighbors
print("Compute batch entropy")
batch_entropy <- list()
for (b in batch) {
  neighbors <- list()
  for (i in as.factor(unique(seurat@meta.data[, b]))) {
    temp <- rownames(seurat@meta.data[seurat@meta.data[ , b] == i, ])
    neighbors[[i]] <- rowSums(as.matrix(seurat@graphs$RNA_nn[, temp]))/30
  }
  neighbors <- as.data.frame(neighbors)
  ## Compute entropy per cell
  optimum <- table(seurat@meta.data[, b]) / ncol(seurat)
  batch_entropy[[b]] <- apply(neighbors, 1, Entropy)
  batch_entropy[[b]] <- batch_entropy[[b]] / Entropy(optimum)
}
batch_entropy <- do.call(cbind, batch_entropy)

## Save batch variables with entropy < 2
#batch_var <- list()
#for (i in colnames(batch_entropy)) {
#  if (median(as.numeric(batch_entropy[, i])) < 2) {
#    batch_var[[i]] <- median(as.numeric(batch_entropy[, i]))
#  }
#}
#batch_var <- batch_var[!duplicated(batch_var)]
#batches <- paste(names(batch_var), sep = ", ")
#print(paste0("Possible batch(es): ", batches))

# Run harmony

if (length(batch >= 1)) {
  print("STEP 2b: RUN HARMONY")
  batch_harmony <- list()
  
  for (i in batch) {
    p0 <- AugmentPlot(DimPlot(seurat, reduction = "umap", group.by = i, pt.size = .1) + NoLegend() + ggtitle("Before harmony"))
    p1 <- AugmentPlot(DimPlot(object = seurat, reduction = "pca", group.by = i, pt.size = .1) + NoLegend())
    seurat_corrected <- seurat %>%
      RunHarmony(i, plot_convergence = FALSE, verbose = verbose) %>%
      RunUMAP(reduction = "harmony", dims = 1:data$pca_dims, a = .5, b = 1.2, verbose = verbose) %>%
      FindNeighbors(dims = 1:2, k.param = 30, reduction = "umap", verbose = verbose)
    
    p3 <- AugmentPlot(DimPlot(object = seurat_corrected, reduction = "harmony", group.by = i, pt.size = .1) + NoLegend())
    p2 <- AugmentPlot(DimPlot(seurat_corrected, reduction = "umap", group.by = i, pt.size = .1) + NoLegend() + ggtitle("After harmony"))
    p <- (p0 | p2) / (p1 | p3)
    ggsave(plot = p, filename = paste0("temp/QC/Harmony_", i, ".png"))
    
    ## Compute the percentage of batch in cell neighbors
    neighbors <- list()
    for (j in as.factor(unique(seurat_corrected@meta.data[, i]))) {
      temp <- rownames(seurat_corrected@meta.data[seurat@meta.data[ , i] == j, ])
      neighbors[[j]] <- rowSums(as.matrix(seurat_corrected@graphs$RNA_nn[, temp]))/30
    }
    neighbors <- as.data.frame(neighbors)
    ## Compute entropy per cell
    optimum <- table(seurat@meta.data[, i]) / ncol(seurat)
    batch_harmony[[i]] <- apply(neighbors, 1, Entropy)
    batch_harmony[[i]] <- batch_harmony[[i]] / Entropy(optimum)
  }
  batch_harmony <- do.call(cbind, batch_harmony)
  ## Plot entropy over all batches
  #batch_entropy$harmony <- "Before"
  batch_entropy <- batch_entropy %>%
    as.data.frame() %>%
    mutate(harmony = "Before") %>%
    tibble::rownames_to_column("cell") %>%
    gather("batch", "entropy", -cell, -harmony)
  
  batch_harmony <- batch_harmony %>%
    as.data.frame() %>%
    mutate(harmony = "After") %>%
    tibble::rownames_to_column("cell") %>%
    gather("batch", "entropy", -cell, -harmony)
  
  batch_entropy <- rbind(batch_entropy, batch_harmony)
  batch_entropy$harmony <- factor(batch_entropy$harmony, levels = c("Before", "After"))
  
  p <- batch_entropy %>% 
    ggplot(aes(y = as.numeric(entropy), x = batch, col = harmony)) +
    geom_boxplot() +
    scale_y_continuous("Entropy") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(plot = p, filename = "temp/QC/batch_entropy.png", width = 10, height = 10)
  
} #else {
## Plot entropy over all batches
#batch_entropy <- as.data.frame(batch_entropy)
#p <- batch_entropy %>% 
#  tibble::rownames_to_column("cell") %>%
#  gather("batch", "entropy", -cell) %>%
#  ggplot(aes(y = as.numeric(entropy), x = batch)) +
#  geom_boxplot() +
#  scale_y_continuous("Entropy") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave(plot = p, filename = "temp/QC/batch_entropy.png", width = 10, height = 10)
#}

# QC
print("STEP 3: CREATE QC PLOTS")
if (length(batch) >= 1) {
  for (i in batch) {
    p2 <- DimPlot(seurat, reduction = "umap", pt.size = .1, group.by = i, label = TRUE) + NoLegend()
    p3 <- AugmentPlot(VlnPlot(seurat, features = "nFeature_RNA", pt.size = 0.1, group.by = i, log = TRUE)) + 
      NoLegend() +
      scale_y_log10("Genes", expand = c(0,0)) + 
      geom_hline(yintercept = data$QC_feature_min, color = "red") + 
      theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p4 <- AugmentPlot(VlnPlot(seurat, features = "nCount_RNA", pt.size = 0.1, group.by = i, log = TRUE)) + 
      NoLegend() + 
      scale_y_log10("Counts", expand = c(0,0)) + 
      theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    if ("percent_mt" %in% colnames(seurat@meta.data)) {
      p5 <- AugmentPlot(VlnPlot(seurat, features = "percent_mt", pt.size = 0.1, group.by = i)) + 
        NoLegend() +
        geom_hline(yintercept = data$QC_mt_max, color = "red") + 
        scale_y_continuous("Mito", expand = c(0,0)) +
        theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text())
    } else {
      seurat@meta.data$percent_mt <- 0
      p5 <- AugmentPlot(VlnPlot(seurat, features = "percent_mt", pt.size = 0.1, group.by = i)) + 
        NoLegend() +
        geom_hline(yintercept = data$QC_mt_max, color = "red") + 
        scale_y_continuous("Mito", expand = c(0,0)) +
        theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text())
    }
    if ("CD3D" %in% rownames(seurat)) {
      p6 <- AugmentPlot(FeaturePlot(seurat, features = "CD3D", pt.size = .1)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_text())
    } else if ("CD68" %in% rownames(seurat)) {
      p6 <- AugmentPlot(FeaturePlot(seurat, features = "CD68", pt.size = .1)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_text())
    } else if ("CLDN5" %in% rownames(seurat)) {
      p6 <- AugmentPlot(FeaturePlot(seurat, features = "CLDN5", pt.size = .1)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_text())
    } else {
      p6 <- AugmentPlot(DimPlot(seurat, group.by = i))
    }
    p <- (p_lbw + p2) / (p3 + p5) / (p4 + p6)
    ggsave(plot = p, filename = paste0("temp/QC/QC_", i, ".png"))
  }
} else {
  print("No batches found in metadata!")
  p1 <- AugmentPlot(VlnPlot(seurat, features = "nFeature_RNA", pt.size = 0.1, log = TRUE)) + 
    NoLegend() +
    scale_y_log10("Genes", expand = c(0,0)) + 
    geom_hline(yintercept = data$QC_feature_min, color = "red") + 
    theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p2 <- AugmentPlot(VlnPlot(seurat, features = "nCount_RNA", pt.size = 0.1, log = TRUE)) + 
    NoLegend() + 
    scale_y_log10("Counts", expand = c(0,0)) + 
    theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p3 <- AugmentPlot(VlnPlot(seurat, features = "percent_mt", pt.size = 0.1)) + 
    NoLegend() +
    geom_hline(yintercept = data$QC_mt_max, color = "red") + 
    scale_y_continuous("Mito", expand = c(0,0)) +
    theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text())
  p <- p_lbw /  p1 / p2 / p3
  ggsave(plot = p, filename = "temp/QC/QC.png")
}

## Save data.json
if (length(batch) == 1) {
  data$batch = "patient"
} else {
  data$batch = FALSE
}
data <- toJSON(data)

if (!file.exists("out/data.json")) {
  write(data, "out/data.json")
} else {
  print("data.json already exists, writing to copy")
  write(data, "out/data_copy.json")
}

