#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
batch = args[1] #batch variable in the metadata slot if no batch fill in empty string
QC_feature_min = as.numeric(args[2]) #Minimal features threshold
QC_mt_max = as.numeric(args[3]) #Maximum mitochondrial content threshold
pca_dims = as.numeric(args[4]) #Amount of PCA dimensions to use
malignant = args[5]
data = "temp/data.rds" #If data is already normalized or not, stored by check_seurat.R
features_var = 2000 #Amount of variable features to select
cluster_resolution = c(1) #At which resolutions to cluster the data
object_path = "temp/raw.rds" #_raw.rds file
cellMarker_path = "/home/jordi_camps/IMMUcan/TME_markerGenes.xlsx"
chetahClassifier_path = "/home/jordi_camps/IMMUcan/CHETAH_reference_updatedAnnotation.RData"

# Make and set directories
dir <- getwd()
setwd(dir)
ifelse(!dir.exists("temp"), dir.create("temp"), "temp/ already exists")
ifelse(!dir.exists("out"), dir.create("out"), "out/ already exists")

# Load packages and set environment
library(Seurat)
library(SingleCellExperiment)
library(CHETAH)
library(harmony)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(openxlsx)
library(pheatmap)
library(DescTools)
library(copykat)
RNGkind(sample.kind = "Rounding")
set.seed(111)

# Recreate seurat object

seurat <- readRDS(object_path)
data <- readRDS(data)
if (batch == "none") {
  print("NO BATCH SPECIFIED => NO INTEGRATION")
  batch = "orig.ident"
}

# QC

print("STEP 1: QC")
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

# Entropy

print("STEP 2: MEASURING BATCH EFFECT")
if (data$norm == FALSE) {
  seurat <- Seurat::NormalizeData(seurat, verbose = TRUE)
}
seurat <- seurat %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = features_var, verbose = TRUE) %>% 
  ScaleData(verbose = TRUE) %>% 
  RunPCA(pc.genes = seurat@var.genes, npcs = pca_dims+20, verbose = TRUE) %>%
  RunUMAP(dims = 1:pca_dims, a = .5, b = 1.2, verbose = TRUE) %>%
  FindNeighbors(dims = 1:2, k.param = 30, reduction = "umap", verbose = TRUE)

p4 <- ElbowPlot(seurat, ndims = pca_dims+20) + geom_vline(xintercept = pca_dims, color = "red") + ylab("STDEV PCA") + theme(axis.title.x = element_blank())
p <- p4 / p1 / p2 / p3
ggsave(plot = p, filename = "out/QC.png")

p0 <- AugmentPlot(DimPlot(seurat, reduction = "umap", group.by = batch, pt.size = .1) + 
                    NoLegend() + 
                    ggtitle("Before harmony"))

## Compute the percentage of batch in cell neighbors
neighbors <- list()
for (i in unique(seurat@meta.data[, batch])) {
  temp <- rownames(seurat@meta.data[seurat@meta.data[ , batch] == i, ])
  neighbors[[i]] <- rowSums(as.matrix(seurat@graphs$RNA_nn[, temp]))/30
}
neighbors <- as.data.frame(neighbors)

## Compute entropy per cell
entropy <- list()
for (i in 1:nrow(neighbors)) {
  entropy[[rownames(neighbors)[i]]] <- Entropy(neighbors[i, ])
}
entropy <- as.matrix(entropy)
median_entropy <- median(as.numeric(entropy[,1]))

entropy <- as.data.frame(entropy)
p <- ggplot(entropy, aes(y = as.numeric(V1), x = 1)) +
  geom_boxplot() +
  scale_y_continuous("Entropy")
ggsave(plot = p, filename = "out/entropy.png", width = 2, height = 4)

# Harmony

if (median_entropy < 1) {
  print("STEP 3: INTEGRATING BATCH")
  p1 <- AugmentPlot(DimPlot(object = seurat, reduction = "pca", pt.size = .1, group.by = batch) + NoLegend())
  p2 <- AugmentPlot(VlnPlot(object = seurat, features = "PC_1", group.by = batch, pt.size = .1) + NoLegend() + theme(plot.title = element_blank()))
  
  seurat <- seurat %>% 
    RunHarmony(batch, plot_convergence = FALSE)
  
  p3 <- AugmentPlot(DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by = batch) + NoLegend())
  p4 <- AugmentPlot(VlnPlot(object = seurat, features = "harmony_1", group.by = batch, pt.size = .1) + NoLegend() + theme(plot.title = element_blank()))
  
  # Dimensionality reduction and clustering
  
  print("STEP 4: CLUSTERING")
  seurat <- seurat %>% 
    RunUMAP(reduction = "harmony", dims = 1:pca_dims, a = .5, b = 1.2, verbose = FALSE) %>%
    RunTSNE(reduction = "harmony", dims = 1:pca_dims, check_duplicates = FALSE)  %>%
    FindNeighbors(reduction = "harmony", dims = 1:pca_dims, verbose = FALSE) %>% 
    FindClusters(resolution = cluster_resolution, verbose = FALSE) %>% 
    identity()
  
  p5 <- AugmentPlot(DimPlot(seurat, reduction = "umap", group.by = batch, pt.size = .1) + 
                      NoLegend() + 
                      ggtitle("After harmony"))
  p <- (p0 | p5) / (p1 | p3) / (p2 | p4)
  ggsave(plot = p, filename = "out/Harmony.png")
} else {
  print("STEP 4: CLUSTERING")
  seurat <- seurat %>% 
    RunUMAP(reduction = "pca", dims = 1:pca_dims, a = .5, b = 1.2, verbose = FALSE) %>%
    RunTSNE(reduction = "pca", dims = 1:pca_dims, check_duplicates = FALSE)  %>%
    FindNeighbors(reduction = "pca", dims = 1:pca_dims, verbose = FALSE) %>% 
    FindClusters(resolution = cluster_resolution, verbose = FALSE) %>% 
    identity()
}

# Supervised annotation

print("STEP 5: SUPERVISED ANNOTATION")
load(chetahClassifier_path)
input <- SingleCellExperiment(assays = list(counts = seurat[["RNA"]]@data),
                              reducedDims = SimpleList(TSNE = seurat@reductions$umap@cell.embeddings))
input <- CHETAHclassifier(input = input, ref_cells = reference, n_genes = 500, thresh = 0.05)

p1 <- PlotCHETAH(input, return = TRUE) 
#nodes <- c("Node1" = "Immune", "Node2" = "Immune", "Node3" = "Lymphoid", "Node4" = "Lymphoid", "Node5" = "NKT", "Node6" = "T", "Node7" = "T", "Node8" = "Myeloid", "Node9" = "Macro/DC", "Node10"= "Stromal", "Node11" = "Stromal")
#input$celltype_CHETAH <- plyr::revalue(input$celltype_CHETAH, replace = nodes[names(nodes) %in% input$celltype_CHETAH])
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

# copyKat

if (malignant == TRUE) {
  print("STEP 6: CALLING COPY NUMBER ABBERATIONS")
  counts <- as.matrix(seurat[["RNA"]]@counts)
  if (ncol(seurat) > 50000) {
    samples <- sample(ncol(seurat), 50000, replace = FALSE)
    seurat_sampled <- seurat[, samples]
    normal_cells <- rownames(seurat_sampled@meta.data[seurat_sampled$annotation_CHETAH %in% c("CD8 T cell", "Macrophage"), ])
    if (length(normal_cells) > 100) {
      copykat.test <- copykat(rawmat=counts, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.15, distance="euclidean", norm.cell.names=normal_cells, n.cores=4)
    } else {
      copykat.test <- copykat(rawmat=counts, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.15, distance="euclidean", norm.cell.names="", n.cores=4)
    }
    pred.test <- data.frame(copykat.test$prediction)
    seurat_sampled@meta.data <- merge(seurat_sampled@meta.data, pred.test[, "copykat.pred", drop = FALSE], by = "row.names", all = TRUE) %>% 
      tibble::column_to_rownames("Row.names")
    p1 <- DimPlot(seurat_sampled, group.by = "copykat.pred")
    p2 <- FeaturePlot(seurat_sampled, features = "EPCAM")
    p3 <- DimPlot(seurat_sampled, group.by = "seurat_clusters", label = TRUE) + NoLegend()
    p <- p1 + p2 + p3
    ggsave(plot = p, filename = "out/copyKat_umap.pdf", height = 5, width = 15)
  }
}

# Plot cell markers

print("STEP 7: CREATING MARKER GENE PLOTS")
cell.markers <- readxl::read_excel(cellMarker_path)
markers <- list()
for (i in as.character(na.omit(unique(cell.markers$cell_type)))) {
  temp <- rownames(seurat)[rownames(seurat) %in% na.omit(cell.markers[cell.markers$cell_type == i, "gene", drop = TRUE])]
  if (length(temp) > 0) {
    markers[[i]] <- temp
  }
}

temp <- AddModuleScore(seurat, features = markers)
p <- DotPlot(temp, features = colnames(temp@meta.data)[grepl("Cluster[[:digit:]]", colnames(temp@meta.data))], cluster.idents = TRUE) + scale_x_discrete(labels = names(markers)) + RotatedAxis()
ggsave(plot = p, filename = "temp/Dotplot_seuratClusters_geneModules.png", dpi = 100, height = 12, width = 12)
p0 <- DotPlot(seurat, features = unique(cell.markers$gene), group.by = "seurat_clusters", cluster.idents = TRUE) + coord_flip() + NoLegend()

##CHETAH recommendation
fraction_chetah <- seurat@meta.data %>%
  group_by(seurat_clusters, annotation_CHETAH) %>%
  tally() %>%
  mutate(fraction_CHETAH = n/sum(n)) %>%
  select(-n) %>%
  arrange(desc(fraction_CHETAH), .by_group = TRUE) %>%
  slice_head(n = 1)

##copykat recommendation
fraction_copykat <- seurat@meta.data %>%
  group_by(seurat_clusters, copykat.pred) %>%
  tally() %>%
  mutate(fraction_copykat = n/sum(n)) %>%
  select(-n) %>%
  arrange(desc(fraction_copykat), .by_group = TRUE) %>%
  slice_head(n = 1)

annotation <- inner_join(fraction_chetah, fraction_copykat, by = "seurat_clusters")
annotation$abbreviation <- ""

##Create annotation.xlsx
if (!file.exists("out/annotation.xlsx")) {
  write.xlsx(x = annotation, "out/annotation.xlsx")
} else {
  print("Not overwriting annotation.xlsx, saving as copy")
  write.xlsx(x = annotation, "out/annotation_copy.xlsx")
}

ggsave(plot = p0, filename = "temp/Dotplot_seuratClusters_genes.png", dpi = 100, height = 12, width = 12)
p1 <- AugmentPlot(DimPlot(seurat, group.by = "seurat_clusters", label = TRUE, label.size = 12))
cell.markers <- cell.markers[cell.markers$gene %in% rownames(seurat), ]
for (type in unique(cell.markers$category)) {
  p2 <- FeaturePlot(seurat, features = unique(cell.markers[cell.markers$category == type, ]$gene), pt.size = .1)
  p3 <- DotPlot(seurat, features = unique(cell.markers[cell.markers$category == type, ]$gene), group.by = "seurat_clusters", cluster.idents = TRUE) + coord_flip() + NoLegend()
  layout <- "
  ACC
  BBB
  BBB
  "
  p <- p1 + p2 + p3 + plot_layout(design = layout)
  ggsave(plot = p, filename = paste0("temp/", type, ".png"), height = 30, width = 20, dpi = 100)
}

temp <- table(seurat$seurat_clusters, seurat$annotation_CHETAH)
temp <- apply(temp, 1, function(x) x / sum(x))
pheatmap::pheatmap(temp, filename = "temp/cluster_comparison.pdf")

# Summary statistics

print("STEP 8: CREATING SUMMARY STATISTICS")
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
  "Cells_after_QC" = ncol(seurat),
  "Entropy" = median_entropy
)
seurat@misc <- list(harmony_summary)
write.csv(x = harmony_summary, file = "out/harmony_summary.csv", row.names = FALSE)

# Save RDS and convert to h5ad with seuratdisk

print("STEP 8: SAVING RESULTS")
saveRDS(seurat, paste0("temp/harmony.rds"))
print("ALL DONE")
