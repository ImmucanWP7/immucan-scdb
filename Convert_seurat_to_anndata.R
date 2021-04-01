#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
object_path = args[1]

suppressPackageStartupMessages({
  library(sceasy)
  library(reticulate)
  use_condaenv('sceasy')
  loompy <- reticulate::import('loompy')
  library(Seurat)
  library(SeuratDisk)
})

seurat <- readRDS(object_path)
#SaveH5Seurat(seurat, filename = "out/harmony.h5Seurat", overwrite = TRUE)
#Convert("out/harmony.h5Seurat", dest = "h5ad", overwrite = TRUE)

# Convert to h5ad with sceasy for immediate use with cellxgene
sceasy::convertFormat(seurat, from="seurat", to="anndata", outFile=gsub(".rds$", ".h5ad", object_path))
