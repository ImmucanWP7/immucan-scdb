suppressPackageStartupMessages({
  library(sceasy)
  library(reticulate)
  use_condaenv('sceasy')
  loompy <- reticulate::import('loompy')
  library(Seurat)
  library(tidyverse)
  library(readxl)
})

tidy_metadata_path <- "/home/jordi_camps/IMMUcan/tidy_metadata.xlsx"
seurat_obj <- "out/harmony.rds"
data <- fromJSON("out/data.json")

meta_cols <- read_excel(tidy_metadata_path)
seurat <- readRDS(seurat_obj)
if (!is.na(data$nSample) & ncol(seurat) > data$nSample) {subsamples <- sample(ncol(seurat), data$nSample, replace = FALSE)}
#colnames(seurat@meta.data) <- tolower(colnames(seurat@meta.data))
print(colnames(seurat@meta.data))

if (any(meta_cols$col_names %in% colnames(seurat@meta.data))) {
  for (i in seq_along(colnames(seurat@meta.data))) {
    hit <-  meta_cols$col_names %in% colnames(seurat@meta.data)[i]
    print(paste0("changing ", colnames(seurat@meta.data)[i], " to ", meta_cols$general[hit]))
    colnames(seurat@meta.data)[i] <- meta_cols$general[hit]
  }
} else {
    stop("No meta.data columns to tidy up")
}

Idents(seurat) <- seurat$annotation_minor
saveRDS("out/harmony.rds")
seurat@meta.data <- seurat@meta.data[, !grepl("RNA_snn_res|abbreviation|cell_id|cell.id|cell_ontology", colnames(seurat@meta.data))]

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

