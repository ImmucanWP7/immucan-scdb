#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
seurat_obj = args[1] #path of seurat object
QC_feature_min = 250 #Minimal features threshold
QC_mt_max = 20 #Maximum mitochondrial content threshold
pca_dims = 30 #Amount of PCA dimensions to use
features_var = 2000 #Amount of variable features to select
verbose = FALSE

dir <- getwd()
setwd(dir)
ifelse(!dir.exists("temp"), dir.create("temp"), "temp/ already exists")
ifelse(!dir.exists("out"), dir.create("out"), "out/ already exists")

library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(DescTools)
library(tidyr)
library(tibble)
library(jsonlite)

print("STEP 1: CHECKING SEURAT OBJECT")

seurat_temp <- readRDS(seurat_obj)
seurat <- CreateSeuratObject(counts = seurat_temp[["RNA"]]@counts, meta.data = seurat_temp@meta.data, min.cells = 10, min.features = 200)

print(paste0("nCell = ", ncol(seurat)))
print(paste0("nGene = ", nrow(seurat)))

if (sum(colnames(seurat) == rownames(seurat@meta.data)) == ncol(seurat)) {
  print("Cell IDs linked correctly")
} else {
  stop("Cell IDs linked uncorrectly")
}

gapdh <- grepl("GAPDH|gapdh", rownames(seurat))
data <- list()
if (sum(grepl("//.", seurat[["RNA"]]@counts[gapdh, ])) == 0) {
  data$norm <- FALSE
  print("Raw counts supplied")
} else {
  data$norm <- TRUE
  print("Normalized counts supplied, data won't be normalized")
}

## Add mitochondrial fraction information
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-")

seurat[["RNA"]]@counts[1:5,1:5]
dplyr::glimpse(seurat@meta.data)
saveRDS(seurat, "temp/raw.rds")

# Batch
print("STEP 2: ESTIMATING BATCH VARIABLE")

## Sample object to max 20k cells
if (ncol(seurat) > 20000) {
  seurat_sampled <- seurat[, sample(colnames(seurat), 20000, replace = FALSE)]
} else {
  seurat_sampled <- seurat
}

## Remove bad quality cells
seurat_sampled <- subset(seurat_sampled, subset = nFeature_RNA > QC_feature_min & percent.mt < QC_mt_max)

## Select potential batch columns from meta.data
batch <- seurat_sampled@meta.data[, sapply(seurat_sampled@meta.data, class) %in% c("character", "factor")]
batch <- batch[, sapply(sapply(batch, unique), length) != 1]

## Create nearest neighbour graph
if (data$norm == FALSE) {
  seurat_sampled <- Seurat::NormalizeData(seurat_sampled, verbose = verbose)
}
seurat_sampled <- seurat_sampled %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = features_var, verbose = verbose) %>% 
  ScaleData(verbose = verbose) %>% 
  RunPCA(pc.genes = seurat_sampled@var.genes, npcs = pca_dims+20, verbose = verbose) %>%
  RunUMAP(dims = 1:pca_dims, a = .5, b = 1.2, verbose = verbose) %>%
  FindNeighbors(dims = 1:2, k.param = 30, reduction = "umap", verbose = verbose)

p <- ElbowPlot(seurat_sampled, ndims = pca_dims+20) + geom_vline(xintercept = pca_dims, color = "red") + ylab("STDEV PCA") + theme(axis.title.x = element_blank())
ggsave(plot = p, filename = "temp/Elbow.png")

## Compute the percentage of batch in cell neighbors
neighbors <- list()
batch_entropy <- list()

for (b in colnames(batch)) {
  for (i in unique(seurat_sampled@meta.data[, b])) {
    temp <- rownames(seurat_sampled@meta.data[seurat_sampled@meta.data[ , b] == i, ])
    neighbors[[i]] <- rowSums(as.matrix(seurat_sampled@graphs$RNA_nn[, temp]))/30
  }
  neighbors <- as.data.frame(neighbors)
  
  ## Compute entropy per cell
  entropy <- list()
  for (i in 1:nrow(neighbors)) {
    entropy[[rownames(neighbors)[i]]] <- Entropy(neighbors[i, ])
  }
  batch_entropy[[b]] <- as.matrix(entropy)
}

## Plot entropy over all batches
batch_entropy <- as.data.frame(batch_entropy)
p <- batch_entropy %>% 
  tibble::rownames_to_column("cell") %>%
  gather("batch", "entropy", -cell) %>%
  ggplot(aes(y = as.numeric(entropy), x = batch)) +
  geom_boxplot() +
  scale_y_continuous("Entropy")
ggsave(plot = p, filename = "temp/batch_entropy.png", width = 6, height = 6)

## Save batch variables with entropy < 2
batch_var <- list()
for (i in colnames(batch_entropy)) {
  if (median(as.numeric(batch_entropy[, i])) < 2) {
    print(paste0("Possible batch column: ", i))
    batch_var[[i]] <- median(as.numeric(batch_entropy[, i]))
  }
}

# QC
print("STEP 3: CREATE QC PLOTS")
if (length(batch_var) >= 1) {
  for (i in names(batch_var)) {
    p1 <- DimPlot(seurat_sampled, reduction = "pca", pt.size = 1, group.by = i, label = TRUE) + NoLegend()
    p2 <- DimPlot(seurat_sampled, reduction = "umap", pt.size = 1, group.by = i, label = TRUE) + NoLegend()
    p3 <- AugmentPlot(VlnPlot(seurat, features = "nFeature_RNA", pt.size = 0.1, group.by = i, log = TRUE)) + 
      NoLegend() +
      scale_y_log10("Genes", expand = c(0,0)) + 
      geom_hline(yintercept = QC_feature_min, color = "red") + 
      theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p4 <- AugmentPlot(VlnPlot(seurat, features = "nCount_RNA", pt.size = 0.1, group.by = i, log = TRUE)) + 
      NoLegend() + 
      scale_y_log10("Counts", expand = c(0,0)) + 
      theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p5 <- AugmentPlot(VlnPlot(seurat, features = "percent.mt", pt.size = 0.1, group.by = i)) + 
      NoLegend() +
      geom_hline(yintercept = QC_mt_max, color = "red") + 
      scale_y_continuous("Mito", expand = c(0,0)) +
      theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text())
    p <- p1 + p2 / p3 / p4 / p5
    ggsave(plot = p, filename = paste0("temp/QC_", i, ".png"))
  }
} else {
  print("No batch effect in dataset!")
  p1 <- AugmentPlot(VlnPlot(seurat, features = "nFeature_RNA", pt.size = 0.1, log = TRUE)) + 
    NoLegend() +
    scale_y_log10("Genes", expand = c(0,0)) + 
    geom_hline(yintercept = QC_feature_min, color = "red") + 
    theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p2 <- AugmentPlot(VlnPlot(seurat, features = "nCount_RNA", pt.size = 0.1, log = TRUE)) + 
    NoLegend() + 
    scale_y_log10("Counts", expand = c(0,0)) + 
    theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p3 <- AugmentPlot(VlnPlot(seurat, features = "percent.mt", pt.size = 0.1)) + 
    NoLegend() +
    geom_hline(yintercept = QC_mt_max, color = "red") + 
    scale_y_continuous("Mito", expand = c(0,0)) +
    theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text())
  p <- p1 / p2 / p3
  ggsave(plot = p, filename = "temp/QC.png")
}

## Save data.json
if (length(names(batch_var)) >=1) {
  data$batch = names(batch_var)
} else {
  data$batch = FALSE
}
data$QC_feature_min = QC_feature_min
data$QC_mt_max =  QC_mt_max
data$pca_dims = pca_dims
data$features_var = features_var
data$metadata = colnames(batch)
data$annotation = c("seurat_clusters","annotation_CHETAH","annotation_major","annotation_immune","annotation_minor", colnames(seurat@meta.data)[grepl("Cluster|cluster|author|Author|Annotation|annotation", colnames(seurat@meta.data))])
data$malignant = FALSE
data <- toJSON(data)
write(data, "out/data.json")
