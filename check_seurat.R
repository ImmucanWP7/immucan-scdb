#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
seurat_obj = args[1] #path of seurat object

library(Seurat)

seurat <- readRDS(seurat_obj)

print(paste0("CELL_ID = ", colnames(seurat)[1]))
print(paste0("CELL_NUMBER = ", ncol(seurat)))
print(paste0("GENE_ID = ", rownames(seurat)[1]))
print(paste0("GENE_NUMBER = ", nrow(seurat)))

if (sum(colnames(seurat) == rownames(seurat@meta.data)) == ncol(seurat)) {
  print("CELL ID CORRECTLY LINKED")
} else {
  ptint("CELL ID NOT CORRECTLY LINKED")
}

gapdh <- grepl("GAPDH|gapdh", rownames(seurat))
if (sum(grepl("//.", seurat[["RNA"]]@counts[gapdh, ])) == 0) {
  print("RAW COUNTS SUPPLIED")
} else {
  print("NORMALIZED COUNTS SUPPLIED")
}

seurat[["RNA"]]@counts[1:5,1:5]
dplyr::glimpse(seurat@meta.data)