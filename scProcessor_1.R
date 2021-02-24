#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
data = "out/data.json" #If data is already normalized or not, stored by check_seurat.R
object_path = "temp/raw.rds" #_raw.rds file
cellMarker_path = "/home/jordi_camps/IMMUcan/TME_markerGenes.xlsx"
chetahClassifier_path = "/home/jordi_camps/IMMUcan/CHETAH_reference_updatedAnnotation.RData"
verbose = FALSE

# Make and set directories
dir <- getwd()
print(dir)
setwd(dir)
if (!file.exists("temp/raw.rds")) {
  stop("first run check_seurat.R")
}

# Load packages and set environment
suppressPackageStartupMessages({
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
  library(future)
  library(jsonlite)
})
suppressWarnings(RNGkind(sample.kind = "Rounding"))
set.seed(111)
options(future.globals.maxSize= 150000*1024^2)
plan("multisession", workers = 4)

# Recreate seurat object

seurat <- readRDS(object_path)
data <- fromJSON("out/data.json")
if (length(data$batch) > 1) {stop("More than one batch specified, select the correct batch")}
if (!"cluster_resolution" %in% names(data)) {data$cluster_resolution = seq(from = 0.4, to = 3, by = 0.2)}

# QC

print("STEP 1a: QC")
cells_before_QC <- ncol(seurat)
seurat <- subset(seurat, subset = nFeature_RNA > data$QC_feature_min & percent.mt < data$QC_mt_max)
if (data$norm == FALSE) {
  seurat <- Seurat::NormalizeData(seurat, verbose = verbose)
}
seurat <- suppressWarnings(seurat %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = data$features_var, verbose = verbose) %>% 
  ScaleData(verbose = verbose) %>% 
  RunPCA(pc.genes = seurat@var.genes, npcs = data$pca_dims+20, verbose = verbose) %>%
  RunUMAP(dims = 1:data$pca_dims, a = .5, b = 1.2, verbose = verbose))

# Harmony

if (data$batch != FALSE) {
  print("STEP 1b: INTEGRATING BATCH")
  p0 <- AugmentPlot(DimPlot(seurat, reduction = "umap", group.by = data$batch, pt.size = .1) + 
                      NoLegend() + 
                      ggtitle("Before harmony"))
  p1 <- AugmentPlot(DimPlot(object = seurat, reduction = "pca", pt.size = .1, group.by = data$batch) + NoLegend())
  p2 <- AugmentPlot(VlnPlot(object = seurat, features = "PC_1", group.by = data$batch, pt.size = .1) + NoLegend() + theme(plot.title = element_blank()))
  
  seurat <- seurat %>% 
    RunHarmony(data$batch, plot_convergence = FALSE, verbose = verbose)
  
  p3 <- AugmentPlot(DimPlot(object = seurat, reduction = "harmony", pt.size = .1, group.by = data$batch) + NoLegend())
  p4 <- AugmentPlot(VlnPlot(object = seurat, features = "harmony_1", group.by = data$batch, pt.size = .1) + NoLegend() + theme(plot.title = element_blank()))
  
}

# Dimensionality reduction and clustering
print("STEP 2: CLUSTERING")
  
if (data$batch != FALSE) {
  seurat <- seurat %>% 
    RunUMAP(reduction = "harmony", dims = 1:data$pca_dims, a = .5, b = 1.2, verbose = verbose) %>%
    RunTSNE(reduction = "harmony", dims = 1:data$pca_dims, check_duplicates = FALSE)  %>%
    FindNeighbors(reduction = "harmony", dims = 1:data$pca_dims, verbose = verbose) %>% 
    FindClusters(resolution = data$cluster_resolution, verbose = verbose)
  p5 <- AugmentPlot(DimPlot(seurat, reduction = "umap", group.by = data$batch, pt.size = .1) + 
                      NoLegend() + 
                      ggtitle("After harmony"))
  p <- (p0 | p5) / (p1 | p3) / (p2 | p4)
  ggsave(plot = p, filename = "out/Harmony.png")
} else {
  seurat <- seurat %>% 
    RunUMAP(reduction = "pca", dims = 1:data$pca_dims, a = .5, b = 1.2, verbose = verbose) %>%
    RunTSNE(reduction = "pca", dims = 1:data$pca_dims, check_duplicates = FALSE)  %>%
    FindNeighbors(reduction = "pca", dims = 1:data$pca_dims, verbose = verbose) %>% 
    FindClusters(resolution = data$cluster_resolution, verbose = verbose)
}

if (length(data$cluster_resolution) > 1) {
print("Defining optimal cluster resolution")
  for (i in seq_along(data$cluster_resolution)) {
    if (ncol(seurat) > 20000) {
      samples <- sample(colnames(seurat), 20000, replace = FALSE)
      seurat_sampled <- seurat[, samples]
    } else {
      seurat_sampled <- seurat
    }
    print(paste0("Checking resolution ", data$cluster_resolution[i]))
    Idents(seurat_sampled) <- seurat_sampled[[paste0("RNA_snn_res.", data$cluster_resolution[i])]]
    seurat.markers <- FindAllMarkers(seurat_sampled, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = verbose)
    seurat.markers.unique <- seurat.markers[!duplicated(seurat.markers$gene) & seurat.markers$p_val_adj < 0.05, ]
    if (length(levels(seurat.markers$cluster)) != sum(table(seurat.markers.unique$cluster) > 1)) {
      if (i == 1) {
      print("optimal cluster resolution lower as the minimum, consider 1) checking QC, 2) running a doublet identifier or 3) adapting the resolution")
      seurat$seurat_clusters <- seurat[[paste0("RNA_snn_res.", data$cluster_resolution[i])]]
      data$cluster_resolution <- data$cluster_resolution[[i]]
      } else {
        seurat$seurat_clusters <- seurat[[paste0("RNA_snn_res.", data$cluster_resolution[i-1])]]
        data$cluster_resolution <- data$cluster_resolution[[i-1]]
      }
      break
    }
  }
}
seurat@meta.data <- seurat@meta.data[, !grepl("RNA_snn_res.", colnames(seurat@meta.data))]
Idents(seurat) <- seurat$seurat_clusters #Set seurat_clusters to Idents

# Supervised annotation

print("STEP 3a: SUPERVISED ANNOTATION")
load(chetahClassifier_path)
input <- SingleCellExperiment(assays = list(counts = seurat[["RNA"]]@data),
                              reducedDims = SimpleList(TSNE = seurat@reductions$umap@cell.embeddings))
input <- CHETAHclassifier(input = input, ref_cells = reference, n_genes = 500, thresh = 0.05)
p1 <- PlotCHETAH(input, return = TRUE) 
#nodes <- c("Node1" = "Immune", "Node2" = "Immune", "Node3" = "Lymphoid", "Node4" = "Lymphoid", "Node5" = "NKT", "Node6" = "T", "Node7" = "T", "Node8" = "Myeloid", "Node9" = "Macro/DC", "Node10"= "Stromal", "Node11" = "Stromal")
#input$celltype_CHETAH <- plyr::revalue(input$celltype_CHETAH, replace = nodes[names(nodes) %in% input$celltype_CHETAH])
seurat@meta.data$annotation_CHETAH <- input$celltype_CHETAH
ggsave(plot = p1, filename = "out/CHETAH_classification.pdf", height = 6, width = 12)

##CHETAH recommendation
fraction_chetah <- seurat@meta.data %>%
  group_by(seurat_clusters, annotation_CHETAH) %>%
  tally(name = "nCells_CHETAH") %>%
  mutate(fraction_CHETAH = round(nCells_CHETAH/sum(nCells_CHETAH), digits = 2)) %>%
  select(-nCells_CHETAH) %>%
  arrange(desc(fraction_CHETAH), .by_group = TRUE) %>%
  slice_head(n = 1)

# copyKat

if (data$malignant == TRUE) {
  print("STEP 3b: CALLING COPY NUMBER ABBERATIONS")
  if (ncol(seurat) > 20000) {
    samples <- sample(colnames(seurat), 20000, replace = FALSE)
    seurat_sampled <- seurat[, samples]
  } else {
    seurat_sampled <- seurat
  }
  counts <- as.matrix(seurat_sampled[["RNA"]]@counts)
  if (is.na(data$normal_cells)) {
    normal_cells <- rownames(seurat_sampled@meta.data[seurat_sampled$annotation_CHETAH %in% c("Macrophage"), ])
    print("Running copykat with Macrophages as normal cells")
    copykat.test <- copykat(rawmat=counts, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.15, distance="euclidean", norm.cell.names=normal_cells, n.cores=4)
  } else if (data$normal_cells == FALSE) {
    print("Running copykat without normal cells")
    copykat.test <- copykat(rawmat=counts, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.15, distance="euclidean", norm.cell.names="", n.cores=4)
  } else {
    normal_cells <- rownames(seurat_sampled@meta.data[seurat_sampled$annotation_CHETAH %in% c(data$normal_cells), ])
    print(paste0("Running copykat with ", data$normal_cells, " as normal cells"))
    copykat.test <- copykat(rawmat=counts, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.15, distance="euclidean", norm.cell.names=normal_cells, n.cores=4)
  }
  pred.test <- data.frame(copykat.test$prediction)
  pred.test <- pred.test[, "copykat.pred", drop = FALSE]
  seurat@meta.data <- seurat@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    left_join(pred.test %>% 
                tibble::rownames_to_column("cell"),
              by = "cell") %>%
    tibble::column_to_rownames("cell")
  
  p1 <- DimPlot(seurat, group.by = "copykat.pred")
  p2 <- FeaturePlot(seurat, features = "EPCAM")
  p3 <- DimPlot(seurat, group.by = "seurat_clusters", label = TRUE) + NoLegend()
  p <- p1 + p2 + p3
  ggsave(plot = p, filename = "out/copyKat_umap.pdf", height = 5, width = 15)

  ##copykat recommendation
  fraction_copykat <- seurat@meta.data %>%
    group_by(seurat_clusters, copykat.pred) %>%
    tally(name = "nCells_copykat") %>%
    filter(is.na(copykat.pred) == FALSE) %>%
    mutate(fraction_copykat = round(nCells_copykat/sum(nCells_copykat), digits = 2)) %>%
    arrange(desc(fraction_copykat), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    select(-nCells_copykat, -fraction_copykat)
  
  annotation <- inner_join(fraction_chetah, fraction_copykat, by = "seurat_clusters")
  annotation$abbreviation <- as.character("")
  annotation[annotation$fraction_CHETAH >= .8, "abbreviation"] <- annotation[annotation$fraction_CHETAH >= .8, "annotation_CHETAH"]
  annotation[annotation$copykat.pred == "aneuploid", "abbreviation"] <- "mal"
} else {
  annotation <- fraction_chetah
  annotation$abbreviation <- as.character("")
  annotation[annotation$fraction_CHETAH >= .8, "abbreviation"] <- annotation[annotation$fraction_CHETAH >= .8, "annotation_CHETAH"]
}

##Create annotation.xlsx
if (!file.exists("out/annotation.xlsx")) {
  write.xlsx(x = annotation, "out/annotation.xlsx")
} else {
  print("Not overwriting annotation.xlsx, saving as copy")
  write.xlsx(x = annotation, "out/annotation_copy.xlsx")
}

# Plot cell markers

print("STEP 4: CREATING MARKER GENE PLOTS")
cell.markers <- readxl::read_excel(cellMarker_path)
markers <- list()
for (i in as.character(na.omit(unique(cell.markers$cell_type)))) {
  temp <- rownames(seurat)[rownames(seurat) %in% na.omit(cell.markers[cell.markers$cell_type == i, "gene", drop = TRUE])]
  if (length(temp) > 0) {
    markers[[i]] <- temp
  }
}

#Idents(seurat) <- seurat$seurat_clusters #set seurat_clusters as idents
temp <- AddModuleScore(seurat, features = markers)
p <- DotPlot(temp, features = colnames(temp@meta.data)[grepl("Cluster[[:digit:]]", colnames(temp@meta.data))], group.by = "seurat_clusters", cluster.idents = TRUE) + scale_x_discrete(labels = names(markers)) + RotatedAxis()
ggsave(plot = p, filename = "temp/Dotplot_seuratClusters_geneModules.png", dpi = 100, height = 12, width = 12)

p0 <- DotPlot(seurat, features = unique(cell.markers$gene), group.by = "seurat_clusters", cluster.idents = TRUE) + coord_flip()
ggsave(plot = p0, filename = "temp/Dotplot_seuratClusters_genes.png", dpi = 100, height = 12, width = 12)

p1 <- AugmentPlot(DimPlot(seurat, label = TRUE, label.size = 12))
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

print("STEP 5: CREATING SUMMARY STATISTICS")
harmony_summary = data.frame(
  "Input_file" = object_path,
  "Batch" = data$batch,
  "QC_features_min" = data$QC_feature_min,
  "QC_mito_max" = data$QC_mt_max,
  "Variable_features" = data$features_var,
  "PCA_dimensions" = data$pca_dims,
  "Amount_genes" = nrow(seurat),
  "Genes_detected_per_cell" = median(seurat@meta.data$nFeature_RNA),
  "Cells_before_QC" = cells_before_QC,
  "Cells_after_QC" = ncol(seurat),
  "Cluster_resolution" = data$cluster_resolution
)
seurat@misc <- list(harmony_summary)
write.csv(x = harmony_summary, file = "out/harmony_summary.csv", row.names = FALSE)

# Save RDS and convert to h5ad with seuratdisk

print("STEP 6: SAVING RESULTS")
saveRDS(seurat, paste0("temp/harmony.rds"))
data <- toJSON(data)
write(data, "out/data.json")
print("ALL DONE")
