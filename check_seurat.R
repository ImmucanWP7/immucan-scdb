#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
seurat_obj = args[1] #path of seurat object

dir <- getwd()
setwd(dir)
ifelse(!dir.exists("temp"), dir.create("temp"), FALSE)
ifelse(!dir.exists("out"), dir.create("out"), FALSE)

library(Seurat)

seurat_temp <- readRDS(seurat_obj)
seurat <- CreateSeuratObject(counts = seurat_temp[["RNA"]]@counts, meta.data = seurat_temp@meta.data, min.cells = 10, min.features = 200)

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
data <- list()
if (sum(grepl("//.", seurat[["RNA"]]@counts[gapdh, ])) == 0) {
  data$norm <- FALSE
  print("RAW COUNTS SUPPLIED")
} else {
  data$norm <- TRUE
  print("NORMALIZED COUNTS SUPPLIED")
}
saveRDS(data, "temp/data.rds")

seurat[["RNA"]]@counts[1:5,1:5]
dplyr::glimpse(seurat@meta.data)

saveRDS(seurat, "temp/raw.rds")