#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
seurat_obj = args[1] #path of seurat object
batch_var = args[2] #batch variable if known
verbose = FALSE

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
  library(DescTools)
  library(tidyr)
  library(tibble)
  library(jsonlite)
  library(harmony)
})

print("STEP 1: CHECKING SEURAT OBJECT")

if (is.na(seurat_obj)) {
  seurat_obj <- normalizePath(list.files(pattern = ".rds$"))
  if (length(seurat_obj) != 1) {
    stop("Specify seurat object in arguments")
  }
}
seurat_temp <- readRDS(seurat_obj)
seurat <- CreateSeuratObject(counts = seurat_temp[["RNA"]]@counts, meta.data = seurat_temp@meta.data, min.cells = 10, min.features = 200)

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
  data$annotation = c("seurat_clusters","annotation_CHETAH","annotation_major","annotation_immune","annotation_minor", colnames(seurat@meta.data)[grepl("Cluster|cluster|author|Author|Annotation|annotation|Cell_type|cell_type", colnames(seurat@meta.data))])
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
} else {
  data$norm <- TRUE
  print("Normalized counts supplied, data won't be normalized")
}

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-")
seurat[["RNA"]]@counts[1:5,1:5]
dplyr::glimpse(seurat@meta.data)

# Batch
print("STEP 2a: ESTIMATING BATCH VARIABLES")

## Sample object to max nSample specified
#if (!is.na(data$nSample) & ncol(seurat) > data$nSample) {
#  samples <- sample(colnames(seurat), data$nSample, replace = FALSE)
#  seurat_sampled <- seurat[, samples]
#} else {
#  data$nSample <- NA
  seurat_sampled <- seurat
#}

## Remove bad quality cells
seurat_sampled <- subset(seurat_sampled, subset = nFeature_RNA > data$QC_feature_min & percent.mt < data$QC_mt_max)

## Select potential batch columns from meta.data
meta <- seurat_sampled@meta.data[, sapply(seurat_sampled@meta.data, class) %in% c("character", "factor")] #Select all columns that are factor or character
meta <- meta[, sapply(sapply(meta, unique), length) != 1, drop = FALSE] #Remove all columns that have only one variable
meta <- apply(meta, 2, function(x) gsub("^$|^ $", NA, x)) #Remove all columns with NAs
if (!"metadata" %in% names(data)) {data$metadata = colnames(meta)} #save metadata columns to data.json
if (length(data$batch) > 1) {
  batch <- data$batch
} else if (is.na(data$batch)) {
  batch <- meta[, apply(meta, 2, function(x) !any(is.na(x))), drop = FALSE] #Remove all columns with NAs
  batch <- colnames(batch)
} else {
  batch <- data$batch
}

## Create nearest neighbour graph
if (data$norm == FALSE) {
  seurat_sampled <- Seurat::NormalizeData(seurat_sampled, verbose = verbose)
} else {seurat_sampled[["RNA"]]@data <- seurat_sampled[["RNA"]]@counts}

seurat_sampled <- seurat_sampled %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = data$features_var, verbose = verbose) %>% 
  ScaleData(verbose = verbose) %>% 
  RunPCA(pc.genes = seurat_sampled@var.genes, npcs = data$pca_dims+20, verbose = verbose) %>%
  RunUMAP(dims = 1:data$pca_dims, a = .5, b = 1.2, verbose = verbose) %>%
  FindNeighbors(dims = 1:2, k.param = 30, reduction = "umap", verbose = verbose)

p_lbw <- ElbowPlot(seurat_sampled, ndims = data$pca_dims+20) + geom_vline(xintercept = data$pca_dims, color = "red") + ylab("STDEV PCA") + theme(axis.title.x = element_blank())

## Compute the percentage of batch in cell neighbors
batch_entropy <- list()
for (b in batch) {
  neighbors <- list()
  for (i in unique(seurat_sampled@meta.data[, b])) {
    temp <- rownames(seurat_sampled@meta.data[seurat_sampled@meta.data[ , b] == i, ])
    neighbors[[i]] <- rowSums(as.matrix(seurat_sampled@graphs$RNA_nn[, temp]))/30
  }
  neighbors <- as.data.frame(neighbors)
  ## Compute entropy per cell
  batch_entropy[[b]] <- apply(neighbors, 1, Entropy)
}
batch_entropy <- do.call(cbind, batch_entropy)

## Save batch variables with entropy < 2
batch_var <- list()
for (i in colnames(batch_entropy)) {
  if (median(as.numeric(batch_entropy[, i])) < 2) {
    batch_var[[i]] <- median(as.numeric(batch_entropy[, i]))
  }
}
batch_var <- batch_var[!duplicated(batch_var)]
batches <- paste(names(batch_var), sep = ", ")
print(paste0("Possible batch(es): ", batches))

# Run harmony

if (length(batch_var >= 1)) {
  print("STEP 2b: RUN HARMONY")
  batch_harmony <- list()
  p0 <- AugmentPlot(DimPlot(seurat_sampled, reduction = "umap", group.by = i, pt.size = .1) + NoLegend() + ggtitle("Before harmony"))
  p1 <- AugmentPlot(DimPlot(object = seurat_sampled, reduction = "pca", group.by = i, pt.size = .1) + NoLegend())
  
  for (i in names(batch_var)) {
    seurat_sampled <- seurat_sampled %>%
      RunHarmony(i, plot_convergence = FALSE, verbose = verbose) %>%
      RunUMAP(reduction = "harmony", dims = 1:data$pca_dims, a = .5, b = 1.2, verbose = verbose) %>%
      FindNeighbors(dims = 1:2, k.param = 30, reduction = "umap", verbose = verbose)
    
    p3 <- AugmentPlot(DimPlot(object = seurat_sampled, reduction = "harmony", group.by = i, pt.size = .1) + NoLegend())
    p2 <- AugmentPlot(DimPlot(seurat_sampled, reduction = "umap", group.by = i, pt.size = .1) + NoLegend() + ggtitle("After harmony"))
    p <- (p0 | p2) / (p1 | p3)
    ggsave(plot = p, filename = paste0("temp/QC/Harmony_", i, ".png"))
    
    ## Compute the percentage of batch in cell neighbors
    neighbors <- list()
    for (j in unique(seurat_sampled@meta.data[, i])) {
      temp <- rownames(seurat_sampled@meta.data[seurat_sampled@meta.data[ , i] == j, ])
      neighbors[[j]] <- rowSums(as.matrix(seurat_sampled@graphs$RNA_nn[, temp]))/30
    }
    neighbors <- as.data.frame(neighbors)
    ## Compute entropy per cell
    batch_harmony[[i]] <- apply(neighbors, 1, Entropy)
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
} else {
  ## Plot entropy over all batches
  batch_entropy <- as.data.frame(batch_entropy)
  p <- batch_entropy %>% 
    tibble::rownames_to_column("cell") %>%
    gather("batch", "entropy", -cell) %>%
    ggplot(aes(y = as.numeric(entropy), x = batch)) +
    geom_boxplot() +
    scale_y_continuous("Entropy") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(plot = p, filename = "temp/QC/batch_entropy.png", width = 10, height = 10)
}

# QC
print("STEP 3: CREATE QC PLOTS")
if (length(batch_var) >= 1) {
  for (i in names(batch_var)) {
    p2 <- DimPlot(seurat_sampled, reduction = "umap", pt.size = .1, group.by = i, label = TRUE) + NoLegend()
    p3 <- AugmentPlot(VlnPlot(seurat, features = "nFeature_RNA", pt.size = 0.1, group.by = i, log = TRUE)) + 
      NoLegend() +
      scale_y_log10("Genes", expand = c(0,0)) + 
      geom_hline(yintercept = data$QC_feature_min, color = "red") + 
      theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p4 <- AugmentPlot(VlnPlot(seurat, features = "nCount_RNA", pt.size = 0.1, group.by = i, log = TRUE)) + 
      NoLegend() + 
      scale_y_log10("Counts", expand = c(0,0)) + 
      theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p5 <- AugmentPlot(VlnPlot(seurat, features = "percent.mt", pt.size = 0.1, group.by = i)) + 
      NoLegend() +
      geom_hline(yintercept = data$QC_mt_max, color = "red") + 
      scale_y_continuous("Mito", expand = c(0,0)) +
      theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text())
    if ("CD3D" %in% rownames(seurat)) {
      p6 <- AugmentPlot(FeaturePlot(seurat_sampled, features = "CD3D", pt.size = .1)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_text())
    } else if ("CD68" %in% rownames(seurat)) {
      p6 <- AugmentPlot(FeaturePlot(seurat_sampled, features = "CD68", pt.size = .1)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_text())
    } else if ("CLDN5" %in% rownames(seurat)) {
      p6 <- AugmentPlot(FeaturePlot(seurat_sampled, features = "CLDN5", pt.size = .1)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_text())
    } else {
      p6 <- AugmentPlot(DimPlot(seurat_sampled, group.by = i))
    }
    p <- (p_lbw + p2) / (p3 + p5) / (p4 + p6)
    ggsave(plot = p, filename = paste0("temp/QC/QC_", i, ".png"))
  }
} else {
  print("No batch effect in dataset!")
  p1 <- AugmentPlot(VlnPlot(seurat, features = "nFeature_RNA", pt.size = 0.1, log = TRUE)) + 
    NoLegend() +
    scale_y_log10("Genes", expand = c(0,0)) + 
    geom_hline(yintercept = data$QC_feature_min, color = "red") + 
    theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p2 <- AugmentPlot(VlnPlot(seurat, features = "nCount_RNA", pt.size = 0.1, log = TRUE)) + 
    NoLegend() + 
    scale_y_log10("Counts", expand = c(0,0)) + 
    theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p3 <- AugmentPlot(VlnPlot(seurat, features = "percent.mt", pt.size = 0.1)) + 
    NoLegend() +
    geom_hline(yintercept = data$QC_mt_max, color = "red") + 
    scale_y_continuous("Mito", expand = c(0,0)) +
    theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text())
  p <- p_lbw /  p1 / p2 / p3
  ggsave(plot = p, filename = "temp/QC/QC.png")
}

## Save data.json
if (length(names(batch_var)) >=1) {
  data$batch = names(batch_var)
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

