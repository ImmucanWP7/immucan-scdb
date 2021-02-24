library(Seurat)
library(tidyverse)
library(vroom)
library(tidyr)

# Read count matrix
read.csv("", header = TRUE, row.names = 1, sep = "\t", nrows = 6)[,1:6] #Read try
counts <- read.csv("", header = TRUE, row.names = 1, sep = "\t")
meta <- data.frame("cell_id" = colnames(counts))
meta <- meta %>% 
  tidyr::separate(col = cell_id, into = c("barcode", "Patient_id"), remove = FALSE) %>%
  select(-barcode) %>%
  tibble::column_to_rownames("cell_id")

# Patient info

# Create Seurat object
seurat <- CreateSeuratObject(counts = counts, meta.data = meta)
saveRDS(seurat, "raw.rds")
