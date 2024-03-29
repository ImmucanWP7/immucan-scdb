#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
object_path = "temp/harmony.rds" #harmony.rds file
if (file.exists("out/annotation.xlsx")) {
  annotationFile_path = "out/annotation.xlsx" #path to annotation file
} else {annotationFile_path = "out/annotation.xls"} #path to annotation file
cellOntology_path = "/home/jordi_camps/IMMUcan/cell_ontology.xlsx"
verbose = FALSE
if (!dir.exists("temp")) {dir.create("temp")}
if (!dir.exists("temp/plots")) {dir.create("temp/plots")}
if (!dir.exists("out")) {dir.create("out")}
if (!dir.exists("out/plots")) {dir.create("out/plots")}

dir <- getwd()
setwd(dir)
print(dir)
if (!file.exists("temp/harmony.rds")) {
  stop("first run scProcessor_1.R")
}
suppressPackageStartupMessages({
library(sceasy)
library(reticulate)
use_condaenv('sceasy')
loompy <- reticulate::import('loompy')
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(genesorteR)
library(data.table)
library(future)
library(jsonlite)
})
suppressWarnings(RNGkind(sample.kind = "Rounding"))
set.seed(111)
options(future.globals.maxSize= 150000*1024^2)
plan("multiprocess", workers = 12)
seurat <- readRDS(object_path)
data <- fromJSON("out/data.json")
if (!is.na(data$nSample) & ncol(seurat) > data$nSample) {subsamples <- sample(ncol(seurat), data$nSample, replace = FALSE)}

#makeReference, takes a Seurat Object and name of meta data column that contains the clusters. Returns a ranking of genes.
makeReference = function(seuratObj, groupBy) {
  groupBy = which(colnames(seuratObj@meta.data) == groupBy)
  gs = sortGenes(seuratObj@assays$RNA@counts, factor(seuratObj@meta.data[,groupBy]), binarizeMethod = "naive", cores = 12)
  pp = getPValues(gs, numPerm = 5, cores = 1)
  pp = apply(pp$adjpval, 1, function(x) any(x < 0.1))
  mm = getMarkers(gs)
  ref = mm$gene_shannon_index
  ref[!pp] = max(ref[[2]])
  return(sort(ref, decreasing = FALSE))
}

# Annotate

print("STEP 1: LINKING CELL ONTOLOGY")
anno_clust <- readxl::read_excel(annotationFile_path)
#anno_clust <- arrange(anno_clust, seurat_clusters)
new.cluster.ids <- tolower(anno_clust$abbreviation)
names(new.cluster.ids) <- levels(seurat)
seurat <- RenameIdents(seurat, new.cluster.ids)
seurat$abbreviation <- Idents(seurat)
cell_ont <- readxl::read_excel(cellOntology_path)
cell_ont$abbreviation <- tolower(cell_ont$abbreviation)
seurat@meta.data <- seurat@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  left_join(cell_ont, by = "abbreviation") %>%
  tibble::column_to_rownames("cell")
if (any(is.na(seurat$cell_ontology)) == TRUE) {
  print(distinct(seurat@meta.data[is.na(seurat$cell_ontology), c("abbreviation", "cell_ontology")]))
  stop("NOT ALL CELL TYPE ABBREVIATIONS FIT")
} 
Idents(seurat) <- seurat$seurat_clusters

# Remove annotations with less than 10 cells
for (i in data$annotation) {
  temp <- names(table(seurat@meta.data[[i]]))[table(seurat@meta.data[[i]]) <= 10]
  seurat@meta.data[seurat@meta.data[[i]] %in% temp, i] <- NA
}

# Plotting

temp <- colnames(seurat@meta.data)[tolower(colnames(seurat@meta.data)) %in% tolower(data$annotation)]
for (i in temp) {
  if (is.numeric(seurat@meta.data[[i]]) == TRUE) {
    p <- FeaturePlot(seurat, features = i, reduction = "umap")
    ggsave(plot = p, filename = paste0("out/plots/", i, ".png"), dpi = 300, width = 10, height = 10)
  } else if (length(unique(seurat@meta.data[[i]])) <= 20) {
    p <- DimPlot(seurat, reduction = "umap", pt.size = 1, group.by = i, label = TRUE) + ggthemes::scale_color_tableau(palette = "Tableau 20")
    ggsave(plot = p, filename = paste0("out/plots/", i, ".png"), dpi = 300, width = 10, height = 10)
  } else {
    p <- DimPlot(seurat, reduction = "umap", pt.size = 1, group.by = i, label = TRUE)
    ggsave(plot = p, filename = paste0("out/plots/", i, ".png"), dpi = 300, width = 10, height = 10)
  }
}

temp <- colnames(seurat@meta.data)[tolower(colnames(seurat@meta.data)) %in% tolower(data$metadata)]
for (i in temp) {
  p <- ggplot(seurat@meta.data, aes_string(x = "cell_ontology", fill = i)) + geom_bar(position = "fill") + RotatedAxis()
  ggsave(plot = p, filename = paste0("out/plots/", i, ".png"), dpi = 300, width = 10, height = 10)
}

# DE
print("STEP 2: CALCULATING MARKER GENES")
#if (exists("subsamples")) {
#  seurat_sampled <- seurat[, subsamples]
#} else {
  seurat_sampled <- seurat
#}
annoCounts <- list()
for (i in data$annotation) {
  Idents(seurat_sampled) <- seurat_sampled[[i, drop = TRUE]]
  seurat.markers <- FindAllMarkers(seurat_sampled, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = verbose)
  write.table(seurat.markers, paste0("out/DE_", i, ".tsv"), sep = "\t")
  temp <- table(seurat[[i]])
  annoCounts[[i]] <- data.frame(annotation = i, cell_type = temp)
  write.csv(data.table::rbindlist(annoCounts), file = "out/cell_count.csv", row.names = FALSE)
}

# Gene entropy ranking
print("STEP 3: CALCULATING GENE ENTROPY RANKING")
geneIndex <- list()
for (i in data$annotation) {
  if (length(unique(seurat_sampled@meta.data[[i]])) > 1) {
    geneIndex[[i]] <- makeReference(seuratObj = seurat_sampled, groupBy = i)
  }
}
geneIndex <- do.call(cbind, geneIndex)
write.table(geneIndex, "out/gene_index.tsv", row.names = TRUE, sep = "\t")


# Export
print("STEP 4: SAVING RESULTS")
#Seurat
saveRDS(seurat, "out/harmony.rds")
#SaveH5Seurat(seurat, filename = "out/harmony.h5Seurat", overwrite = TRUE)
#Convert("out/harmony.h5Seurat", dest = "h5ad", overwrite = TRUE)

# Export average gene expression over cluster
#write.csv(x = seurat[["RNA"]]@data, file = "out/normCounts.csv", row.names = TRUE)
Idents(seurat) <- seurat$annotation_major
temp <- AverageExpression(seurat, assays = "RNA")
write.table(x = temp$RNA, file = "out/avgExpr_major.tsv", row.names = TRUE, sep = "\t")
Idents(seurat) <- seurat$annotation_immune
temp <- AverageExpression(seurat, assays = "RNA")
write.table(x = temp$RNA, file = "out/avgExpr_immune.tsv", row.names = TRUE, sep = "\t")
Idents(seurat) <- seurat$annotation_minor
temp <- AverageExpression(seurat, assays = "RNA")
write.table(x = temp$RNA, file = "out/avgExpr_minor.tsv", row.names = TRUE, sep = "\t")
Idents(seurat) <- seurat$annotation_CHETAH
temp <- AverageExpression(seurat, assays = "RNA")
write.table(x = temp$RNA, file = "out/avgExpr_CHETAH.tsv", row.names = TRUE, sep = "\t")

Idents(seurat) <- seurat$annotation_minor
seurat@meta.data <- seurat@meta.data[, sapply(sapply(seurat@meta.data, unique), length) != 1, drop = FALSE] #Remove all columns that have only one variable
seurat@meta.data <- seurat@meta.data[, !grepl("RNA_snn_res|abbreviation|cell_id|cell.id", colnames(seurat@meta.data))]

# Convert to h5ad with sceasy for immediate use with cellxgene
sceasy::convertFormat(seurat, from="seurat", to="anndata", outFile= "out/cellxgene.h5ad")

# Export metadata with umap coordinates
write.table(x = cbind(seurat@meta.data, seurat@reductions$umap@cell.embeddings), file = "out/metadata.tsv", row.names = TRUE, sep = "\t")

#Subsample object to 10k cells
if (exists("subsamples")) {seurat <- seurat[, subsamples]}

# Export metadata with umap coordinates
write.table(x = cbind(seurat@meta.data, seurat@reductions$umap@cell.embeddings), file = "out/metadata_10k.tsv", row.names = TRUE, sep = "\t")

# Convert to h5ad with sceasy for immediate use with cellxgene
sceasy::convertFormat(seurat, from="seurat", to="anndata", outFile= "out/cellxgene_10k.h5ad")

#zip and checksum
print("STEP 5: ZIP AND CHECKSUM")
folder_name <- tail(unlist(strsplit(dir, "/")), n=1)
dir.create(folder_name)
out_files <- paste0("out/", list.files("out/"))
file.copy(out_files, folder_name, recursive = TRUE)
zip(paste0(folder_name, ".zip"), folder_name)
checksum <- tools::md5sum(paste0(folder_name, ".zip"))
file.rename(paste0(folder_name, ".zip"), paste0(folder_name, "_-_", checksum, ".zip"))
file.copy(paste0(folder_name, "_-_", checksum, ".zip"), "../")
unlink(paste0(folder_name, "_-_", checksum, ".zip"))
unlink(folder_name, recursive = TRUE)
