dir <- getwd()
setwd(dir)
print(dir)

suppressPackageStartupMessages({
  library(sceasy)
  library(reticulate)
  use_condaenv('sceasy')
  loompy <- reticulate::import('loompy')
  library(Seurat)
  library(tidyverse)
  library(readxl)
  library(jsonlite)
})

tidy_metadata_path <- "/home/jordi_camps/IMMUcan/tidy_metadata.xlsx"
seurat_obj <- "out/harmony.rds"
data <- fromJSON("out/data.json")

meta_cols <- read_excel(tidy_metadata_path)
seurat <- readRDS(seurat_obj)
if (!is.na(data$nSample) & ncol(seurat) > data$nSample) {subsamples <- sample(ncol(seurat), data$nSample, replace = FALSE)}
#colnames(seurat@meta.data) <- tolower(colnames(seurat@meta.data))
glimpse(seurat@meta.data)

if (any(meta_cols$col_names %in% colnames(seurat@meta.data))) {
  change_cols <- colnames(seurat@meta.data)[colnames(seurat@meta.data) %in% meta_cols$col_names]
  for (i in change_cols) {
    clean = FALSE
    hit_1 <- grepl(i, colnames(seurat@meta.data))
    hit_2 <- meta_cols$col_names %in% i
    print(paste0("changing ", i, " to ", meta_cols$general[hit_2]))
    colnames(seurat@meta.data)[hit_1] <- meta_cols$general[hit_2]
  }
} else {
  clean = TRUE
  print("No meta.data columns to tidy up")
}

print("Updating data.json")
if (isFALSE(clean)) {
  change_cols <- data$metadata[data$metadata %in% meta_cols$col_names]
  for (i in change_cols) {
    hit_1 <- grepl(i, data$metadata)
    hit_2 <- meta_cols$col_names %in% i
    data$metadata[hit_1] <- meta_cols$general[hit_2]
  }
  change_cols <- data$annotation[data$annotation %in% meta_cols$col_names]
  for (i in change_cols) {
    hit_1 <- grepl(i, data$annotation)
    hit_2 <- meta_cols$col_names %in% i
    data$annotation[hit_1] <- meta_cols$general[hit_2]
  }
  data <- toJSON(data)
  write(data, "out/data.json")
}

print("Saving objects")
Idents(seurat) <- seurat$annotation_minor
saveRDS(seurat, "out/harmony.rds")
seurat@meta.data <- seurat@meta.data[, sapply(sapply(seurat@meta.data, unique), length) != 1, drop = FALSE] #Remove all columns that have only one variable
seurat@meta.data <- seurat@meta.data[, !grepl("RNA_snn_res|abbreviation|cell_id|cell.id|cell_id", colnames(seurat@meta.data))]

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