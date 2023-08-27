#!/usr/bin/env Rscript

#All plots and end object put in data derived
#out_dir <- gsub("Original", "Derived", dir)
#dir_qc <- paste0(out_dir, "/QC")
#dir_integration <- paste0(out_dir, "/Integration")

# Define directories
ressource_dir <- "./scProcessor_v2/ressources/"
table_dir <- "./scProcessor_v2/Tables/"

dir <- getwd()
out_dir <- dir
dir_qc <- paste0(out_dir, "/QC/")
dir_integration <- paste0(out_dir, "/Integration/")

if (!dir.exists(out_dir)) {dir.create(out_dir)}
if (!dir.exists(dir_qc)) {dir.create(dir_qc)}
if (!dir.exists(dir_integration)) {dir.create(dir_integration)}

# Define some variables
data_json <- paste0(out_dir, "/data.json")
data_json_copy <- paste0(out_dir, "/data_copy.json")
json_examples <- paste0(ressource_dir, "example.json")
tidy_metadata_path <- paste0(table_dir, "tidy_metadata.xlsx")

# Load packages and set environment

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
  library(scater)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(scran)
  library(cowplot)
})

cat("
#######################################################################################################################\n
##  STEP 1: CHECKING SEURAT OBJECT  ###################################################################################\n
#######################################################################################################################"
)

# Load data.json or create standard
if (file.exists(data_json)) {
  
  data <- fromJSON(data_json)
  print("The data.json found in the output dir is used")
  
  } else {
    
    # Search for seurat object
    files <- list.files(path = out_dir, pattern = ".rds$")
    files <- files[!(files %in% "raw_QC.rds")]
    seurat_obj <- normalizePath(files)
    
    # Create data.json file from standard json
    data <- fromJSON(json_examples)
    data <- toJSON(data)
    write(data, data_json)
    
    if (length(seurat_obj) != 1) {
      data$object_path <- paste0(out_dir, "/", seurat_obj)
      stop("Please check the data.json that has been created in the output folder, adjust mandatory fields
           and specify seurat object in data.json")
    } else {
      stop("Please check the data.json that has been created in the output folder and adjust mandatory fields")
    }
}

data_copy <- toJSON(data)
write(data_copy, data_json_copy)
rm(data_copy)
print("The initial data.json file has been save as data_copy.json, but will not further be used in this pipeline")

if (!exists(data$verbose)) {data$verbose = FALSE}

# Read in seurat object 
seurat_obj <- normalizePath(data$object_path)
seurat <- readRDS(seurat_obj)

if ("test" %in% names(data) & data$test) {
  print("Testing script, subsampling data")
  seurat <- seurat[, sample(colnames(seurat), data$nSample)]
}

# Validate seurat object and data.json

lib <- modules::use("~/scProcessoR/Programs/scProcessoR_v2/R/validate_object.R")
lib$validateSeurat(seurat, data)


# Check if data is normalized
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


cat("
#######################################################################################################################\n
##  STEP 2: QUALITY CONTROL  ##########################################################################################\n
#######################################################################################################################"
)

# Add mitochondrial percentages
seurat[["percent_mt"]] <- PercentageFeatureSet(seurat, pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-")
seurat[["percent_ribo"]] <- PercentageFeatureSet(seurat, pattern = "^Rp[sl]|RP[SL]")
seurat[["percent_hemo"]] <- PercentageFeatureSet(seurat, pattern = "^Hb[^(p)]|HB[^(P)]")
dplyr::glimpse(seurat@meta.data)
if (any(is.na(seurat[[data$batch_var, drop = TRUE]]))) {
  print("Removing NAs from batch column")
  seurat <- seurat[, !is.na(seurat[[data$batch_var, drop = TRUE]])]
  seurat[[data$batch_var]] <- droplevels(seurat[[data$batch_var]])
}

# Plot QC

theme_qc <- theme(
  axis.title.x = element_blank()
  , plot.title = element_blank()
  , axis.title.y = element_text()
  , legend.position = "none"
  )

if (length(data$batch_var) >= 1) {
  for (i in data$batch_var) {
    p1 <- AugmentPlot(VlnPlot(seurat, features = "nFeature_RNA", pt.size = 0, group.by = i, log = TRUE)) + 
      scale_y_log10(expand = c(0,0)) + 
      geom_hline(yintercept = data$QC_feature_min, color = "red") + 
      theme_qc
    p2 <- AugmentPlot(VlnPlot(seurat, features = "nCount_RNA", pt.size = 0, group.by = i, log = TRUE)) + 
      scale_y_log10(expand = c(0,0)) + 
      theme_qc
    if (!"percent_mt" %in% colnames(seurat@meta.data)) {
      seurat@meta.data$percent_mt <- 0
    }
      p3 <- AugmentPlot(VlnPlot(seurat, features = "percent_mt", pt.size = 0, group.by = i)) + 
        geom_hline(yintercept = data$QC_mt_max, color = "red") + 
        scale_y_continuous(expand = c(0,0)) +
        theme_qc
    if (!"percent_ribo" %in% colnames(seurat@meta.data)) {
      seurat@meta.data$percent_ribo <- 0
    }
    p4 <- AugmentPlot(VlnPlot(seurat, features = "percent_ribo", pt.size = 0, group.by = i)) + 
      scale_y_continuous(expand = c(0,0)) + 
      theme_qc
    if (!"percent_hemo" %in% colnames(seurat@meta.data)) {
      seurat@meta.data$percent_hemo <- 0
    }
    p5 <- AugmentPlot(VlnPlot(seurat, features = "percent_hemo", pt.size = 0, group.by = i)) + 
      scale_y_continuous(expand = c(0,0)) + 
      theme_qc
    p <- (p1 + p2) / (p3 + p4) / (p5)
    ggsave(plot = p, filename = paste0(dir_qc, i, ".png"), width = 10, height = 10)
  }
}

#Filter cells and genes out on minimal QC thresholds
seurat <- subset(seurat, subset = nFeature_RNA > data$QC_feature_min & percent_mt < data$QC_mt_max)
seurat <- seurat[rowSums(GetAssayData(seurat, slot = "counts")) > 3, ]

cat("
#######################################################################################################################\n
##  STEP 3: DOUBLET FINDER  ###########################################################################################\n
#######################################################################################################################"
)

# Switch seurat to SinglecCelleExperiment
sce <- as.SingleCellExperiment(seurat)
if (!data$norm) {
  sce <- logNormCounts(sce)
}
dec <- modelGeneVar(sce, block = sce[[data$batch_run]])
hvgs = getTopHVGs(dec, n = data$features_var)
sce <- runPCA(sce, subset_row = hvgs)
sce <- runUMAP(sce, pca = data$pca_dims)
sce <- scDblFinder(sce, dims = data$pca_dims, samples = data$batch_run)
table(sce$scDblFinder.class)

# Plot doublets
p_dbl <- plot_grid(plotUMAP(sce, colour_by = "scDblFinder.score"), 
                   plotUMAP(sce, colour_by = "scDblFinder.class"),
                   plotUMAP(sce, colour_by = data$batch_var), 
                   plotUMAP(sce, colour_by = data$batch_run),
                   ncol = 2)
ggsave(plot = p_dbl, filename = paste0(dir_qc, "/DoubletFinder.png"), width = 12, height = 10)

# Add ScDblFinder information to seurat meta.data
seurat <- AddMetaData(seurat, colData(sce) %>%
                        as.data.frame() %>%
                        select(scDblFinder.score, scDblFinder.class))

# Remove bad quality cells
seurat <- subset(seurat, subset = scDblFinder.class == "singlet")

# Save seurat object
saveRDS(seurat, paste0(dir, "/raw_QC.rds"))

# Batch
cat("
#######################################################################################################################\n
##  STEP 4: ESTIMATING BATCH VARIABLES  ###############################################################################\n
#######################################################################################################################"
)

## Subsample datasets larger than 50k cells
if (ncol(seurat) > data$nSample) {
  subset_size <- 0.1 #subsample to 10% of the data
  subset_id <- sample.int(n = ncol(seurat), size = floor(subset_size * ncol(seurat)), replace=FALSE)
  seurat <- seurat[, subset_id]
}

## Create nearest neighbour graph
if (data$norm == FALSE) {
  seurat <- Seurat::NormalizeData(seurat, verbose = verbose)
} else {seurat[["RNA"]]@data <- seurat[["RNA"]]@counts}

seurat <- seurat %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = data$features_var, verbose = verbose) %>% 
  ScaleData(verbose = verbose) %>% 
  RunPCA(pc.genes = seurat@var.genes, npcs = data$pca_dims+20, verbose = verbose) %>%
  RunUMAP(dims = 1:data$pca_dims, verbose = verbose) %>%
  FindNeighbors(dims = 1:2, k.param = 30, reduction = "umap", verbose = verbose)

p_lbw <- ElbowPlot(seurat, ndims = data$pca_dims+20) + geom_vline(xintercept = data$pca_dims, color = "red") + theme(axis.title.x = element_blank())
ggsave(plot = p_lbw, filename = paste0(dir_qc, "/Elbow_plot.png"), width = 5, height = 5)

## Compute the percentage of batch in cell neighbors
print("Compute batch entropy")
  neighbors <- list()
  for (i in as.factor(unique(seurat@meta.data[, data$batch_var]))) {
    temp <- rownames(seurat@meta.data[seurat@meta.data[ , data$batch_var] == i, ])
    neighbors[[i]] <- rowSums(as.matrix(seurat@graphs$RNA_nn[, temp]))/30
  }
  neighbors <- as.data.frame(neighbors)
  ## Compute entropy per cell
  optimum <- table(seurat@meta.data[, data$batch_var]) / ncol(seurat)
  batch_entropy <- apply(neighbors, 1, Entropy)
  batch_entropy <- batch_entropy / Entropy(optimum)

# Run harmony

cat("
#######################################################################################################################\n
##  STEP 5: RUN HARMONY  ##############################################################################################\n
#######################################################################################################################"
)

p0 <- AugmentPlot(DimPlot(seurat, reduction = "umap", group.by = data$batch_var, pt.size = .1) + NoLegend() + ggtitle("Before harmony"))
p1 <- AugmentPlot(DimPlot(object = seurat, reduction = "pca", group.by = data$batch_var, pt.size = .1) + NoLegend())
seurat_corrected <- seurat %>%
  RunHarmony(data$batch_var, plot_convergence = FALSE, verbose = verbose) %>%
  RunUMAP(reduction = "harmony", dims = 1:data$pca_dims, verbose = verbose) %>%
  FindNeighbors(dims = 1:2, k.param = 30, reduction = "umap", verbose = verbose)

p3 <- AugmentPlot(DimPlot(object = seurat_corrected, reduction = "harmony", group.by = data$batch_var, pt.size = .1) + NoLegend())
p2 <- AugmentPlot(DimPlot(seurat_corrected, reduction = "umap", group.by = data$batch_var, pt.size = .1) + NoLegend() + ggtitle("After harmony"))
p <- (p0 | p2) / (p1 | p3)
ggsave(plot = p, filename = paste0(dir_integration, "/Harmony_", data$batch_var, ".png"))

# Compute the percentage of batch in cell neighbors
neighbors <- list()
for (j in as.factor(unique(seurat_corrected@meta.data[, data$batch_var]))) {
  temp <- rownames(seurat_corrected@meta.data[seurat@meta.data[ , data$batch_var] == j, ])
  neighbors[[j]] <- rowSums(as.matrix(seurat_corrected@graphs$RNA_nn[, temp]))/30
}
neighbors <- as.data.frame(neighbors)

# Compute entropy per cell
optimum <- table(seurat@meta.data[, data$batch_var]) / ncol(seurat)
batch_harmony <- apply(neighbors, 1, Entropy)
batch_harmony <- batch_harmony / Entropy(optimum)

# Plot entropy over all batches
batch_entropy <- batch_entropy %>%
  as.data.frame() %>%
  mutate(harmony = "Before") %>%
  tibble::rownames_to_column("cell")
batch_entropy <- plyr::rename(batch_entropy, replace = c("." = "entropy"))

batch_harmony <- batch_harmony %>%
  as.data.frame() %>%
  mutate(harmony = "After") %>%
  tibble::rownames_to_column("cell")
batch_harmony <- plyr::rename(batch_harmony, replace = c("." = "entropy"))

batch_entropy <- rbind(batch_entropy, batch_harmony)
batch_entropy$harmony <- factor(batch_entropy$harmony, levels = c("Before", "After"))

# Boxplot entropy
p <- batch_entropy %>% 
  ggplot(aes(y = as.numeric(entropy), x = harmony, col = harmony)) +
  geom_boxplot() +
  scale_y_continuous("Entropy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = p, filename = paste0(dir_integration, "/Entropy_boxplot.png"), width = 10, height = 10)

# DimPlot harmony
seurat <- AddMetaData(seurat, batch_entropy %>%
  spread(harmony, entropy) %>%
  tibble::column_to_rownames("cell"))
p <- FeaturePlot(seurat, reduction = "umap", features =  c("Before", "After"), pt.size = .1, ncol = 1)
ggsave(plot = p, filename = paste0(dir_integration, "/Entropy_umap_", data$batch_var, ".png"), width = 5, height = 10)

# Recommend integration
if (median(seurat$Before) < .5) {
  if (median(seurat$After) > .7) {
    if (median(seurat$After) - median(seurat$Before) > .3) {
      data$integrate = TRUE
      print("Integration of the data is recomended")
    }
  }
} else {
  data$integrate = FALSE
  print("Integration of the data probably not necessary")
}

# Plot some marker genes
genes <- data$genes
genes_plot <- genes[genes %in% rownames(seurat)]
p <- FeaturePlot(seurat, features = genes_plot, pt.size = .1)
ggsave(plot = p, filename = paste0(dir_qc, "/marker_genes.png"), width = 20, height = 20)

# Save data.json
data <- toJSON(data)
write(data, data_json)
print("Your data.json has beend adjusted based on the QC. 
      You might want to check again before running scProcessor 'Annotate' ")
    
cat("
#######################################################################################################################\n
##  ALL DONE - QC  ####################################################################################################\n
#######################################################################################################################"
)