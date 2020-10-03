#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
batch = args[1] #batch variable in the metadata slot
normalized = args[2] #Boolean to indicate if data is already normalized
QC_feature_min = as.numeric(args[3]) #Minimal features threshold
QC_mt_max = as.numeric(args[4]) #Maximum mitochondrial content threshold
pca_dims = as.numeric(args[5]) #Amount of PCA dimensions to use

features_var = 2000 #Amount of variable features to select
cluster_resolution = c(1) #At which resolutions to cluster the data
object_path = "raw.rds" #_raw.rds file
cellMarker_file = "/gpfs01/home/glanl/scripts/IMMUcan/TME_markerGenes.xlsx"
#garnett_classifier = "/gpfs01/bhcbio/projects/research_studies/20190920_IMMUCan_Public_data/Garnett_train_datasets/NSCLC_Unbiased_Lambrechts/NSCLC_ALL_10X_garnettTrain_garnett_classifier.rds" #path of garnett classifier


# Make and set directories
dir <- getwd()
setwd(dir)
ifelse(!dir.exists("temp"), dir.create("temp"), FALSE)
ifelse(!dir.exists("out"), dir.create("out"), FALSE)


# Load packages and set environment
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)


# QC

seurat <- readRDS(object_path)
cells_before_QC <- ncol(seurat)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-")

p1 <- AugmentPlot(VlnPlot(seurat, features = "nFeature_RNA", pt.size = 0.1, group.by = batch, log = TRUE)) + 
  NoLegend() +
  scale_y_log10("Genes", expand = c(0,0)) + 
  geom_hline(yintercept = QC_feature_min, color = "red") + 
  theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p2 <- AugmentPlot(VlnPlot(seurat, features = "nCount_RNA", pt.size = 0.1, group.by = batch, log = TRUE)) + 
  NoLegend() + 
  scale_y_log10("Counts", expand = c(0,0)) + 
  theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p3 <- AugmentPlot(VlnPlot(seurat, features = "percent.mt", pt.size = 0.1, group.by = batch)) + 
  NoLegend() +
  geom_hline(yintercept = QC_mt_max, color = "red") + 
  scale_y_continuous("Mito", expand = c(0,0)) +
  theme(axis.title.x = element_blank(), plot.title = element_blank(), axis.title.y = element_text())

seurat <- subset(seurat, subset = nFeature_RNA > QC_feature_min & percent.mt < QC_mt_max)
seurat <- seurat[rowSums(as.matrix(seurat[["RNA"]]@counts == 0)), ] # remove gene with 0 total counts

# Prepare

if (normalized == FALSE) {
  seurat <- Seurat::NormalizeData(seurat, verbose = TRUE)
}
seurat <- seurat %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = features_var, verbose=TRUE) %>% 
  ScaleData(verbose = TRUE) %>% 
  RunPCA(pc.genes = seurat@var.genes, npcs = pca_dims+20, verbose = TRUE) %>%
  RunUMAP(dims = 1:pca_dims, a = .5, b = 1.2, verbose = TRUE)

p4 <- ElbowPlot(seurat, ndims = pca_dims+20) + geom_vline(xintercept = pca_dims, color = "red") + ylab("STDEV PCA") + theme(axis.title.x = element_blank())
p <- p4 / p1 / p2 / p3
ggsave(plot = p, filename = "out/QC.png")

p0 <- AugmentPlot(DimPlot(seurat, reduction = "umap", group.by = batch, pt.size = .1) + 
                    NoLegend() + 
                    ggtitle("Before harmony"))


# Harmony

p1 <- AugmentPlot(DimPlot(object = seurat, reduction = "pca", pt.size = .1, group.by = batch) + NoLegend())
p2 <- AugmentPlot(VlnPlot(object = seurat, features = "PC_1", group.by = batch, pt.size = .1) + NoLegend() + theme(plot.title = element_blank()))

seurat <- seurat %>% 
  RunHarmony(batch, plot_convergence = FALSE)

p3 <- AugmentPlot(DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by = batch) + NoLegend())
p4 <- AugmentPlot(VlnPlot(object = seurat, features = "harmony_1", group.by = batch, pt.size = .1) + NoLegend() + theme(plot.title = element_blank()))


# Dimensionality reduction and clustering

seurat <- seurat %>% 
  RunUMAP(reduction = "harmony", dims = 1:pca_dims, a = .5, b = 1.2, verbose = TRUE) %>%
  RunTSNE(reduction = "harmony", dims = 1:pca_dims, check_duplicates = FALSE)  %>%
  FindNeighbors(reduction = "harmony", dims = 1:pca_dims, verbose = TRUE) %>% 
  FindClusters(resolution = cluster_resolution, verbose = TRUE) %>% 
  identity()

p5 <- AugmentPlot(DimPlot(seurat, reduction = "umap", group.by = batch, pt.size = .1) + 
                    NoLegend() + 
                    ggtitle("After harmony"))
p <- (p0 | p5) / (p1 | p3) / (p2 | p4)
ggsave(plot = p, filename = "out/Harmony.png")


# Supervised annotation

#fdata <- data.frame("gene_short_name" = rownames(seurat[["RNA"]]@counts))
#rownames(fdata) <- fdata$gene_short_name
#pdata <- seurat@meta.data
#pdata$garnett_cluster <- pdata$seurat_clusters
#cds <- newCellDataSet(as(seurat[["RNA"]]@counts, "dgCMatrix"), phenoData = new("AnnotatedDataFrame", data = pdata), featureData = new("AnnotatedDataFrame", data = fdata)) %>%
#  estimateSizeFactors(cds)
#cds_classifier <- readRDS(garnett_classifier) #Read Garnett classifier
#cds <- classify_cells(cds, cds_classifier, db = org.Hs.eg.db, cluster_extend = TRUE, cds_gene_id_type = "SYMBOL")
#seurat$Annotation_garnett <- pData(cds)$cluster_ext_type
#seurat$Annotation_garnett <- gsub("Unknown", NA, seurat$Annotation_garnett)
#print(AugmentPlot(DimPlot(seurat, reduction = "umap", pt.size = .1, group.by = "Annotation_garnett", label = TRUE) + 
#                    ggthemes::scale_color_tableau(palette = "Tableau 20") +
#                    ggsave(paste0("annotation_garnett.png"))))
#DimPlot(seurat, group.by = "Annotation_garnett", split.by = "Annotation_garnett", ncol = 4, pt.size = .1) + 
#  ggthemes::scale_color_tableau(palette = "Tableau 20") + 
#  NoLegend()
#Idents(seurat) <- seurat$Annotation_garnett


# Plot cell makers

cell.markers <- readxl::read_excel(cellMarker_file)
p0 <- DotPlot(seurat, features = unique(cell.markers$gene), group.by = "seurat_clusters") + coord_flip() + NoLegend()
ggsave(plot = p0, filename = "temp/Dotplot_seuratClusters.png", dpi = 100, height = 12, width = 12)
p1 <- AugmentPlot(DimPlot(seurat, group.by = "seurat_clusters", label = TRUE, label.size = 12))
cell.markers <- cell.markers[cell.markers$gene %in% rownames(seurat), ]
for (type in unique(cell.markers$category)) {
  p2 <- FeaturePlot(seurat, features = cell.markers[cell.markers$category == type, ]$gene, pt.size = .1)
  p3 <- VlnPlot(seurat, features = cell.markers[cell.markers$category == type, ]$gene, pt.size = .1, group.by = "seurat_clusters")
  p <- (p1 | p2) / p3
  ggsave(plot = p, filename = paste0("temp/", type, ".png"), height = 20, width = 20, dpi = 100)
}


# Summary statistics

harmony_summary = data.frame(
  "Input_file" = object_path,
  "Batch" = batch,
  "QC_features_min" = QC_feature_min,
  "QC_mito_max" = QC_mt_max,
  "Variable_features" = features_var,
  "PCA_dimensions" = pca_dims,
  "Amount_genes" = nrow(seurat),
  "Genes_detected_per_cell" = median(seurat@meta.data$nFeature_RNA),
  "Cells_before_QC" = cells_before_QC,
  "Cells_after_QC" = ncol(seurat)
)
seurat@misc <- list(harmony_summary)
write.csv(x = harmony_summary, file = "out/harmony_summary.csv", row.names = FALSE)

#seurat@meta.data %>%
#  group_by(Annotation_garnett) %>%
#  tally() %>%
#  write.csv(temp, file = "garnett_count.csv", row.names = FALSE)

# Save RDS and convert to h5ad with searatdisk

saveRDS(seurat, paste0("temp/harmony.rds"))

