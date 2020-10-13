#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
batch = args[1] #batch variable in the metadata slot if no batch fill in empty string
QC_feature_min = as.numeric(args[2]) #Minimal features threshold
QC_mt_max = as.numeric(args[3]) #Maximum mitochondrial content threshold
pca_dims = as.numeric(args[4]) #Amount of PCA dimensions to use
integrate == TRUE
data = "temp/data.rds" #If data is already normalized or not, stored by check_seurat.R
features_var = 2000 #Amount of variable features to select
cluster_resolution = c(1) #At which resolutions to cluster the data
object_path = "temp/raw.rds" #_raw.rds file
cellMarker_path = "/gpfs01/home/glanl/scripts/IMMUcan/TME_markerGenes.xlsx"
chetahClassifier_path = "/gpfs01/home/glanl/scripts/IMMUcan/CHETAH_reference_updatedAnnotation.RData"

# Make and set directories
dir <- getwd()
setwd(dir)
ifelse(!dir.exists("temp"), dir.create("temp"), FALSE)
ifelse(!dir.exists("out"), dir.create("out"), FALSE)

# Load packages and set environment
library(Seurat)
library(SingleCellExperiment)
library(CHETAH)
library(harmony)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(WriteXLS)
library(pheatmap)
RNGkind(sample.kind = "Rounding")
set.seed(111)

# Recreate seurat object

seurat <- readRDS(object_path)
data <- readRDS(data)
if (batch == "") {
  batch = "orig.ident"
  integrate = FALSE
  }

# QC

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

# Prepare

if (data$norm == FALSE) {
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

if (integrate == TRUE) {
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
} else {
  seurat <- seurat %>% 
    RunUMAP(reduction = "pca", dims = 1:pca_dims, a = .5, b = 1.2, verbose = TRUE) %>%
    RunTSNE(reduction = "pca", dims = 1:pca_dims, check_duplicates = FALSE)  %>%
    FindNeighbors(reduction = "pca", dims = 1:pca_dims, verbose = TRUE) %>% 
    FindClusters(resolution = cluster_resolution, verbose = TRUE) %>% 
    identity()
}

# Supervised annotation

load(chetahClassifier_path)
input <- SingleCellExperiment(assays = list(counts = seurat[["RNA"]]@data),
                              reducedDims = SimpleList(TSNE = seurat@reductions$umap@cell.embeddings))
input <- CHETAHclassifier(input = input, ref_cells = reference, n_genes = 500, thresh = 0.05)

p1 <- PlotCHETAH(input, return = TRUE) 
nodes <- c("Node1" = "Immune", "Node2" = "Immune", "Node3" = "Lymphoid", "Node4" = "Lymphoid", "Node5" = "NKT", "Node6" = "T", "Node7" = "T", "Node8" = "Myeloid", "Node9" = "Macro/DC", "Node10"= "Stromal", "Node11" = "Stromal")
input$celltype_CHETAH <- plyr::revalue(input$celltype_CHETAH, replace = nodes[names(nodes) %in% input$celltype_CHETAH])
seurat@meta.data$annotation_CHETAH <- input$celltype_CHETAH
ggsave(plot = p1, filename = "out/CHETAH_classification.pdf", height = 6, width = 12)

# Split object
#Tcells <- c("T", "CD4 T cell", "CD8 T cell", "NK", "NKT", "reg. T cell")
#Myeloid <- c("Myeloid", "Macro/DC", "Macrophage", "Dendritic")
#seurat_T <- seurat[, seurat$annotation_CHETAH %in% Tcells]
#seurat_myeloid <- seurat[, seurat$annotation_CHETAH %in% Myeloid]
#seurat_T <- seurat_T %>% 
#  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose=TRUE) %>% 
#  ScaleData(verbose = TRUE) %>% 
#  RunPCA(npcs = 30, verbose = TRUE) %>%
#  RunUMAP(reduction = "harmony", dims = 1:10, a = .5, b = 1.2, verbose = TRUE) %>%
#  #RunTSNE(reduction = "harmony", dims = 1:pca_dims, check_duplicates = FALSE)  #%>%
#  FindNeighbors(reduction = "harmony", dims = 1:10, verbose = TRUE) %>% 
#  FindClusters(resolution = 0.8, verbose = TRUE) %>% 
#  identity()

# Plot cell markers

cell.markers <- readxl::read_excel(cellMarker_path)
markers <- list()
for (i in as.character(na.omit(unique(cell.markers$cell_type)))) {
    temp <- rownames(seurat)[rownames(seurat) %in% na.omit(cell.markers[cell.markers$cell_type == i, "gene"])]
    if (length(temp) > 0) {
      markers[i] <- temp
    }
}

temp <- AddModuleScore(seurat, features = markers)
p <- DotPlot(temp, features = colnames(temp@meta.data)[grepl("Cluster", colnames(temp@meta.data))], cluster.idents = TRUE) + scale_x_discrete(labels = names(markers)) + RotatedAxis()
ggsave(plot = p, filename = "temp/Dotplot_seuratClusters_geneModules.png", dpi = 100, height = 12, width = 12)
p0 <- DotPlot(seurat, features = unique(cell.markers$gene), group.by = "seurat_clusters", cluster.idents = TRUE) + coord_flip() + NoLegend()
WriteXLS(x = list("annotation" = tibble("seurat_clusters" = 0:(length(unique(seurat$seurat_clusters))-1), "abbreviation" = "Fill in")), ExcelFileName = "out/annotation.xls")
ggsave(plot = p0, filename = "temp/Dotplot_seuratClusters_genes.png", dpi = 100, height = 12, width = 12)
p1 <- AugmentPlot(DimPlot(seurat, group.by = "seurat_clusters", label = TRUE, label.size = 12))
cell.markers <- cell.markers[cell.markers$gene %in% rownames(seurat), ]
for (type in unique(cell.markers$category)) {
  p2 <- FeaturePlot(seurat, features = unique(cell.markers[cell.markers$category == type, ]$gene), pt.size = .1, ncol = 5)
  p3 <- DotPlot(seurat, features = unique(cell.markers[cell.markers$category == type, ]$gene), group.by = "seurat_clusters", cluster.idents = TRUE) + coord_flip() + NoLegend()
  layout <- "
  ACC
  BBB
  BBB
  "
  p <- p1 + p2 + p3 + plot_layout(design = layout)
  ggsave(plot = p, filename = paste0("temp/", type, ".png"), height = 20, width = 20, dpi = 100)
}

temp <- table(seurat$seurat_clusters, seurat$annotation_CHETAH)
temp <- apply(temp, 1, function(x) x / sum(x))
pheatmap::pheatmap(temp, filename = "temp/cluster_comparison.pdf")

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

# Save RDS and convert to h5ad with seuratdisk

saveRDS(seurat, paste0("temp/harmony.rds"))

