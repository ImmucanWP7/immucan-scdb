library(Seurat)
library(tidyverse)
library(tidyr)

count_file <- ""
meta_file <- ""
patient_file <- ""

print("Read count matrix")
counts_test <- read.csv(count_file, header = TRUE, row.names = 1, sep = " ", nrows = 6) #Read try
counts_test[,1:6]
counts <- read.csv(count_file, header = TRUE, row.names = 1, sep = " ")

print("Metadata")
meta <- read.csv(meta_file, header = TRUE, sep = "\t", row.names = 1)
print(head(meta))

print("Patient info")
patient <- read.csv(patient_file, header = TRUE, sep = "\t", row.names = 1)

if (ncol(counts) != sum(colnames(counts) == rownames(meta))) {stop("colnames counts not equal to rownames meta.data")}

print("Create Seurat object")
seurat <- CreateSeuratObject(counts = counts, meta.data = meta)
saveRDS(seurat, "raw.rds")