#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
ref_clust = as.character(args[1]) #number of seurat cluster(s) to use a reference
cutoff = as.numeric(args[2]) #cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics

object_path = "temp/harmony.rds" #_raw.rds file
gene_order_file = "/gpfs01/bhcbio/projects/research_studies/20190920_IMMUCan_Public_data/gencode_v35_gene_positions.txt"

# Make and set directories
dir <- getwd()
setwd(dir)
ifelse(!dir.exists("temp"), dir.create("temp"), FALSE)
ifelse(!dir.exists("out"), dir.create("out"), FALSE)


library(Seurat)   
library(infercnv)

# Prepare data
seurat <- readRDS(object_path)
counts_matrix <- as.matrix(seurat[["RNA"]]@counts)

#type can be cell types or cell states or cluster IDs
anno <- seurat$seurat_clusters
anno <- as.matrix(anno)
rownames(anno) <- colnames(seurat)

gene_order_file <- read.table(gene_order_file, sep = "\t",header = FALSE, row.names = 1)

# Create infercnvobject
cnv_obj = CreateInfercnvObject(counts_matrix, 
                               gene_order_file = gene_order_file,
                               annotations_file = anno, 
                               ref_group_names = list(ref_clust), 
                               delim = "\t")

# Run infercnv
cnv_obj = infercnv::run(cnv_obj,
                        cutoff = cutoff,
                        out_dir = "out/", 
                        cluster_by_groups = FALSE, 
                        plot_steps = FALSE,
                        denoise = TRUE,
                        HMM = TRUE,
                        HMM_type = "i3",
                        no_plot = FALSE,
                        no_prelim_plot = FALSE,
                        k_obs_groups = NULL,
                        analysis_mode = "subclusters",
                        num_threads = 5)
