#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
batch = args[1] #batch variable in the metadata slot
object_path = "temp/raw.rds" #_raw.rds file
data = "temp/data.rds" #If data is already normalized or not, stored by check_seurat.R
QC_feature_min = 200 #Minimal features threshold
QC_mt_max = 15 #Maximum mitochondrial content threshold
pca_dims = 30 #Amount of PCA dimensions to use
features_var = 2000 #Amount of variable features to select
#batch_columns = "patient|Patient|Sample|sample|plate|Plate"

library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)

dir <- getwd()
setwd(dir)
ifelse(!dir.exists("temp"), dir.create("temp"), FALSE)
ifelse(!dir.exists("out"), dir.create("out"), FALSE)

# QC

seurat <- readRDS(object_path)
data <- readRDS(data)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-")
batch <- tolower(batch)

for (i in colnames(seurat@meta.data)[grepl(batch, colnames(seurat@meta.data))]) {
  p1 <- AugmentPlot(VlnPlot(seurat, features = "nFeature_RNA", pt.size = 0.1, group.by = i, log = TRUE)) + 
    NoLegend() +
    scale_y_log10("Genes", expand = c(0,0)) + 
    geom_hline(yintercept = QC_feature_min, color = "red") + 
    theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p2 <- AugmentPlot(VlnPlot(seurat, features = "nCount_RNA", pt.size = 0.1, group.by = i, log = TRUE)) + 
    NoLegend() + 
    scale_y_log10("Counts", expand = c(0,0)) + 
    theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p3 <- AugmentPlot(VlnPlot(seurat, features = "percent.mt", pt.size = 0.1, group.by = i)) + 
    NoLegend() +
    geom_hline(yintercept = QC_mt_max, color = "red") + 
    scale_y_continuous("Mito", expand = c(0,0)) +
    theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text())
  p <- p1 / p2 / p3
  ggsave(plot = p, filename = paste0("temp/QC_", i, ".png"))
}

seurat <- subset(seurat, subset = nFeature_RNA > QC_feature_min & percent.mt < QC_mt_max)

# Prepare

if (data$norm == FALSE) {
  seurat <- Seurat::NormalizeData(seurat, verbose = TRUE)
}
seurat <- seurat %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = features_var, verbose = TRUE) %>% 
  ScaleData(verbose = TRUE) %>% 
  RunPCA(pc.genes = seurat@var.genes, npcs = pca_dims+20, verbose = TRUE) %>%
  RunUMAP(dims = 1:pca_dims, a = .5, b = 1.2, verbose = TRUE)

p4 <- ElbowPlot(seurat, ndims = pca_dims+20) + geom_vline(xintercept = pca_dims, color = "red") + ylab("STDEV PCA") + theme(axis.title.x = element_blank())
ggsave(plot = p4, filename = "temp/Elbow.png")

# Plotting

for (i in colnames(seurat@meta.data)[grepl(batch, colnames(seurat@meta.data))])  {
      p1 <- DimPlot(seurat, reduction = "pca", pt.size = 1, group.by = i, label = TRUE) + NoLegend()
      p2 <- DimPlot(seurat, reduction = "umap", pt.size = 1, group.by = i, label = TRUE) + NoLegend()
      p <- p1 / p2
      ggsave(plot = p, filename = paste0("temp/", i, ".png"))
}
