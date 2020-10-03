#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
object_path = "temp/harmony.rds" #harmony.rds file
annotationFile_path = "annotation.xlsx" #path to annotation file
cellOntology_path = "/gpfs01/home/glanl/scripts/IMMUcan/cell_ontology.xlsx"
meta_cols_umap = c("age", "patient", "biopsy", "tissue", "sample", "seurat_clusters", "annotation_authors", "nCount_RNA", "nFeature_RNA", "annotation_major", "annotation_immune", "annotation_minor")
meta_cols_barplot = c("patient", "biopsy", "tissue", "treatment", "treatment_prior", "treatment_response", "treatment_timepoint", "disease_stage")

library(sceasy)
library(reticulate)
use_condaenv('sceasy')
loompy <- reticulate::import('loompy')
library(Seurat)
#library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(genesorteR)
dir <- getwd()
setwd(dir)
seurat <- readRDS(object_path)

#makeReference, takes a Seurat Object and name of meta data column that contains the clusters. Returns a ranking of genes.
makeReference = function(seuratObj, groupBy) {
  
  groupBy = which(colnames(seuratObj@meta.data) == groupBy)
  gs = sortGenes(seuratObj@assays$RNA@counts, factor(seuratObj@meta.data[,groupBy]), binarizeMethod = "naive", cores = 1)
  pp = getPValues(gs, numPerm = 5, cores = 1)
  pp = apply(pp$adjpval, 1, function(x) any(x < 0.1))
  mm = getMarkers(gs)
  ref = mm$gene_shannon_index
  ref[!pp] = max(ref[[2]])
  return(sort(ref, decreasing = FALSE))
  
}

# Annotate

anno_clust <- readxl::read_excel(annotationFile_path)
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
Idents(seurat) <- seurat$cell_ontology


# Plotting

temp <- colnames(seurat@meta.data)[tolower(colnames(seurat@meta.data)) %in% tolower(meta_cols_umap)]
for (i in temp) {
  if (is.numeric(seurat@meta.data[[i]]) == TRUE) {
    p <- FeaturePlot(seurat, features = i, reduction = "umap")
    ggsave(plot = p, filename = paste0("out/", i, ".png"), dpi = 300, width = 7, height = 6)
  } else if (length(unique(seurat@meta.data[[i]])) <= 20) {
    p <- DimPlot(seurat, reduction = "umap", pt.size = 1, group.by = i, label = TRUE) + ggthemes::scale_color_tableau(palette = "Tableau 20")
    ggsave(plot = p, filename = paste0("out/", i, ".png"), dpi = 300, width = 7, height = 6)
  } else {
    p <- DimPlot(seurat, reduction = "umap", pt.size = 1, group.by = i, label = TRUE)
    ggsave(plot = p, filename = paste0("out/", i, ".png"), dpi = 300, width = 7, height = 6)
  }
}

temp <- colnames(seurat@meta.data)[tolower(colnames(seurat@meta.data)) %in% tolower(meta_cols_barplot)]
for (i in temp) {
    p <- ggplot(seurat@meta.data, aes_string(x = "cell_ontology", fill = i)) + geom_bar(position = "fill") + RotatedAxis()
    ggsave(plot = p, filename = paste0("out/", i, ".png"), dpi = 300, width = 7, height = 6)
}


# DE

seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.table(seurat.markers, paste("out/DE_genes.tsv", sep = "_"), sep = "\t")

seurat@meta.data %>%
  group_by(cell_ontology) %>%
  tally() %>%
  write.csv(file = "out/cell_count.csv", row.names = FALSE)


# Gene entropy ranking

geneIndex <- makeReference(seuratObj = seurat, groupBy = "cell_ontology")
write.table(geneIndex, "out/gene_index.tsv", row.names = TRUE, sep = "\t")


# Export

#Seurat
saveRDS(seurat, "out/harmony.rds")
#SaveH5Seurat(seurat, filename = paste(gsub("_raw.rds", "", object_path), "harmony.h5Seurat", sep = "_"), overwrite = TRUE)
#Convert(paste(gsub("_raw.rds", "", object_path), "harmony.h5Seurat", sep = "_"), dest = "h5ad", overwrite = TRUE)

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
# Export metadata with umap coordinates
write.table(x = cbind(seurat@meta.data, as.data.frame(seurat@reductions$umap@cell.embeddings)), file = "out/metadata.tsv", row.names = TRUE, sep = "\t")

Idents(seurat) <- seurat$cell_ontology
seurat@meta.data <- seurat@meta.data[, !grepl("RNA_snn_res|abbreviation|cell_id|cell.id|orig.ident|cell_ontology", colnames(seurat@meta.data))]

# Convert to h5ad with sceasy for immediate use with cellxgene
sceasy::convertFormat(seurat, from="seurat", to="anndata", outFile= "out/cellxgene.h5ad")
