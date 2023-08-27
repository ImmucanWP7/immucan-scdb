#!/usr/bin/env Rscript


# Define directories
ressource_dir <- "./scProcessor_v2/ressources/"
table_dir <- "./scProcessor_v2/Tables/"

dir <- getwd()
out_dir <- dir
dir_anno <- paste0(out_dir, "/Annotation/")
dir_copykat <- paste0(out_dir, "/Copykat/")

if (!dir.exists(out_dir)) {dir.create(out_dir)}
if (!dir.exists(dir_anno)) {dir.create(dir_anno)}
if (!dir.exists(dir_copykat)) {dir.create(dir_copykat)}

# Define some variables
data_json <- paste0(out_dir, "/data.json")
json_examples <- paste0(ressource_dir, "example.json")
cellMarker_path <- paste0(table_dir, "TME_markerGenes.xlsx")

# Load packages and set environment

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(dplyr)
  library(openxlsx)
  library(DescTools)
  library(copykat)
  library(future)
  library(jsonlite)
  library(tidyr)
  library(SeuratDisk)
})

suppressWarnings(RNGkind(sample.kind = "Rounding"))
set.seed(111)
options(future.globals.maxSize= 150000*1024^2)
plan("multisession", workers = 4)

# Load data.json or create standard
if (file.exists(data_json)) {
  
  data <- fromJSON(data_json)
  print("The data.json found in the output dir is used")
  
} else {
  
  # Create data.json file from standard json
  data <- fromJSON(json_examples)
  data <- toJSON(data)
  write(data, data_json)
  stop("Please check the data.json that has been created in the output folder and adjust mandatory fields")

}

if (!"cluster_resolution" %in% names(data)) {data$cluster_resolution = c(.3, .5, 1)}
if (!exists(data$verbose)) {data$verbose = FALSE}

### Load seurat object ###

path_seurat_rawQC <- paste0(dir, "/raw_QC.rds")
if (file.exists(path_seurat_rawQC)) {
  
  if(data$noQC == TRUE) {
    print("scProcessor 'QC' has been run before, but you specified 'noQC = true' in the data.json file.
          scProcessor Annotate will run on the original sample provided in 'object_path' of the data.json")
    seurat <- readRDS(data$object_path)
  } else {
    print("scProcessor has been run before and the rawQC.rds found in the output folder will be used for analysis")
    seurat <- readRDS(path_seurat_rawQC)
  }
  }
  
if (!file.exists(path_seurat_rawQC)) {
  
  if(data$noQC == FALSE){
    stop("No rawQC.rds file has beend found in the output directory \n\n
         Please run scProcessor 'QC' first or specify 'noQC = true' in the data.json file to skipp QC")
  } else {
    print("scProcessor Annotate will be run on the .rds file specified in the data.json and no QC will be performed
          If you don't want this please run scProcessor 'QC' first and specify 'noQC = false' in the data.json file ")
    seurat <- readRDS(data$object_path)
  }
}


if ("test" %in% names(data) & data$test) {
  print("Testing script, subsampling data; Also only three items from 'batch_var' will be selected")
  keep_runvar <- sample(unique(seurat[[data$batch_var, drop = TRUE]]), 3)
  Idents(seurat) <- data$batch_var
  seurat <- subset(seurat, idents = keep_runvar)
  seurat <- seurat[, sample(colnames(seurat), data$nSample)]
}

# Validate seurat object and data.json

lib <- modules::use("~/scProcessoR/Programs/scProcessoR_v2/R/validate_object.R")
lib$validateSeurat(seurat, data)
  
cat("
#######################################################################################################################\n
##  STEP 1: INTEGRATING BATCH  ########################################################################################\n
#######################################################################################################################"
)

#Prepare
if (data$norm == FALSE) {
  seurat <- Seurat::NormalizeData(seurat, verbose = verbose)
  } else if (any(dim(seurat@assays$RNA@data) == 0)) {
    seurat[["RNA"]]@data <- seurat[["RNA"]]@counts
    print("'norm' has been set to false in the daata.json and no data has been found in the 'data' slot of the object \n
          It is assumed that normalised counts are instead stored in the 'counts' slot")
  }

seurat <- suppressWarnings(seurat %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = data$features_var, verbose = verbose) %>% 
  ScaleData(verbose = verbose) %>% 
  RunPCA(pc.genes = seurat@var.genes, npcs = data$pca_dims+20, verbose = verbose) %>%
  RunUMAP(dims = 1:data$pca_dims, verbose = verbose))

# Harmony
if (data$integrate) {
  reduction <- "harmony"
  seurat <- suppressWarnings(seurat %>% 
    RunHarmony(data$batch_var, plot_convergence = FALSE, verbose = verbose))
} else {
  reduction <- "pca"
  print("No integration necessary")
}

cat("
#######################################################################################################################\n
##  STEP 2: DIMENSIONALITY REDUCTION AND CLUSTERING  ###################################################################\n
#######################################################################################################################"
)

seurat <- seurat %>% 
  RunUMAP(reduction = reduction, dims = 1:data$pca_dims, verbose = verbose) %>%
  FindNeighbors(reduction = reduction, dims = 1:data$pca_dims, verbose = verbose) %>% 
  FindClusters(resolution = data$cluster_resolution, verbose = verbose)

cat("
#######################################################################################################################\n
##  STEP 3: ATLAS MAPPING  ############################################################################################\n
#######################################################################################################################"
)

lib <- modules::use("~/scProcessoR/Programs/scProcessoR_v2/R/ctmapping.R")

path_ref_list <- list(
  path_ref_l1 = "~/scProcessoR/Programs/scProcessoR_v2/reference_bundle/level1"
  , path_ref_l2_lymphoid = "~/scProcessoR/Programs/scProcessoR_v2/reference_bundle/level2/level2-lymphoid"
  )
path_ref_list

res_mapping <- lib$mapCtLevels(
  se_query = seurat,             # Query object as Seurat object
  path_ref_list = path_ref_list, # Named list for references. For details see function description; for level 1 "path_ref_l1" is sufficient
  batchvar_query = data$batch_var,     # Batch variable of the query data, regressed out with harmony 
  query_assay = "RNA",
  query_slot = "data",
  annotation_col_ref = "ct_id",  # Annotation column to map to in the references (have to be the same in all references)
  k = 5,
  save_dir = dir_anno,  # Here all results will be saved (table with predictions & plots)
  save_all = TRUE                             # if set to TRUE, all mapping outputs on every level are saved
  )

# Add cell type prediction to meta.data
rownames(res_mapping) <- res_mapping$cell
res_mapping$cell <- NULL
res_mapping$cell_type_pred <- as.character(res_mapping$cell_type_pred_l2)
seurat <- AddMetaData(seurat, res_mapping)

# Create summary of predicted cell types per seurat cluster
fraction <- seurat@meta.data %>%
  group_by(seurat_clusters, cell_type_pred) %>%
  tally(name = "cell_type_pred_count") %>%
  mutate(cell_type_pred_percentage = round(cell_type_pred_count / sum(cell_type_pred_count), 2)) %>%
  arrange(desc(cell_type_pred_percentage)) %>%
  slice_head(n = 1)

cat("
#######################################################################################################################\n
##  STEP 3: CALLING CNA  ##############################################################################################\n
#######################################################################################################################"
)

# copykat

if (data$malignant == TRUE) {
  copykat_list <- list()
  pred_list <- list()
  #Loop over patient
  setwd(dir_copykat)
  for (patient in unique(seurat[[data$batch_var, drop = T]])) {
    seurat_patient <- seurat[, seurat[[data$batch_var, drop = T]] == patient]
    #Subset object if more than 50000 cells per patient
    if (ncol(seurat_patient) > 50000) {
      subsamples <- sample(ncol(seurat_patient), 50000, replace = FALSE)
      seurat_patient <- seurat_patient[, subsamples]
    }
    counts <- as.matrix(seurat_patient[["RNA"]]@counts)
    if (is.null(data$normal_cells)) {
      normal_cells <- rownames(seurat_patient@meta.data[seurat_patient$cell_type_pred %in% c("Macrophage"), ])
      print("Running copykat with Macrophages as normal cells")
    } else if (data$normal_cells == FALSE) {
      print("Running copykat without normal cells")
      normal_cells <- ""
    } else {
      normal_cells <- rownames(seurat_patient@meta.data[seurat_patient$cell_type_pred %in% c(data$normal_cells), ])
      print(paste0("Running copykat with ", data$normal_cells, " as normal cells"))
    }
    copykat_list[[patient]] <- copykat(rawmat = counts,
                            id.type = "S", 
                            ngene.chr = 5, 
                            win.size = 25, 
                            KS.cut = 0.15, 
                            distance = "euclidean", 
                            norm.cell.names = normal_cells, 
                            n.cores = 4,
                            plot.genes = F,
                            sam.name = patient)
    pred_list[[patient]] <- data.frame(copykat_list[[patient]]$prediction)
    pred_list[[patient]] <- pred_list[[patient]][, "copykat.pred", drop = FALSE]
  }
  setwd(out_dir)

pred <- bind_rows(pred_list, .id = data$batch_var)
seurat <- AddMetaData(seurat, pred)

# Plot copykat output
p1 <- DimPlot(seurat, group.by = "copykat.pred")
p2 <- DimPlot(seurat, group.by = "seurat_clusters", label = TRUE) + NoLegend()
if ("Patient" %in% colnames(seurat@meta.data)) {
  p3 <- DimPlot(seurat, group.by = "Patient", label = TRUE) + NoLegend()
} else {p3 <- ggplot() + theme_void()}
if ("EPCAM" %in% rownames(seurat)) {
  p4 <- FeaturePlot(seurat, features = "EPCAM")
} else {p4 <- ggplot() + theme_void()}
p <- (p1 + p2) / (p3 + p4)
ggsave(plot = p, filename = paste0(dir_copykat, "copyKat_umap.png"), height = 10, width = 15)

# copykat recommendation
fraction_copykat <- seurat@meta.data %>%
  group_by(seurat_clusters, copykat.pred) %>%
  tally(name = "copykat.pred_count") %>%
  mutate(copykat.pred_percentage = round(copykat.pred_count / sum(copykat.pred_count), 2)) %>%
  arrange(desc(copykat.pred_percentage)) %>%
  slice_head(n = 1)
fraction <- left_join(fraction, fraction_copykat)
}

cat("
#######################################################################################################################\n
##  STEP 5: LINKING TO BAYER ONTOLOGY  ################################################################################\n
#######################################################################################################################"
)

# Add annotation summary
fraction$annotation_summary <- fraction$cell_type_pred
if (data$malignant) {
  fraction[fraction$cell_type_pred == "Epithelial" & 
             fraction$copykat.pred == "aneuploid" & 
             !is.na(fraction$copykat.pred), "annotation_summary"
           ] <- "Malignant"
}
metadata <- seurat@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  left_join(fraction %>%
              select(seurat_clusters, annotation_summary)) %>%
  tibble::column_to_rownames("cell_id")
seurat <- AddMetaData(seurat, metadata)

fraction$My_annotation <- "Fill in"
write.xlsx(x = fraction, paste0(dir_anno, "/annotation.xlsx"))

p <- AugmentPlot(DimPlot(seurat, label = TRUE, label.size = 8, group.by = "annotation_summary"))
ggsave(plot = p, 
       filename = paste0(dir_anno, "UMAP_annotation_summary.png"), 
       height = 10, 
       width = 10)


### BAYER ontology ###



cat("
#######################################################################################################################\n
##  STEP 6: CREATING MARKER GENE PLOTS  ###############################################################################\n
#######################################################################################################################"
)

# Plot cell markers

cell.markers <- readxl::read_excel(cellMarker_path)
markers <- list()
for (i in as.character(na.omit(unique(cell.markers$cell_type)))) {
  temp <- rownames(seurat)[rownames(seurat) %in% na.omit(cell.markers[cell.markers$cell_type == i, "gene", drop = TRUE])]
  if (length(temp) > 0) {
    markers[[i]] <- temp
  }
}

temp <- AddModuleScore(seurat, features = markers)
p <- DotPlot(temp, features = colnames(temp@meta.data)[grepl("Cluster[[:digit:]]", colnames(temp@meta.data))],
             group.by = "seurat_clusters", cluster.idents = TRUE) + 
  scale_x_discrete(labels = names(markers)) + 
  theme(axis.text.y = element_text(size = 8)) + 
  RotatedAxis()
ggsave(plot = p, filename = paste0(dir_anno, "Dotplot_seuratClusters_geneModules.png"), height = 12, width = 12)

for (group in c("seurat_clusters", data$metadata)){
  p0 <- DotPlot(seurat, features = unique(cell.markers$gene), 
                group.by = group, cluster.idents = TRUE) +  
    theme(axis.text.y = element_text(size = 8)) + 
    coord_flip()
  ggsave(plot = p0, filename = paste0(dir_anno, "Dotplot_", group, "_genes.png"),
         height = 12, width = 12)
}


p1 <- AugmentPlot(DimPlot(seurat, label = TRUE, label.size = 8))
cell.markers <- cell.markers[cell.markers$gene %in% rownames(seurat), ]
for (type in unique(cell.markers$category)) {
  markers <- unique(cell.markers[cell.markers$category == type, ]$gene)
  if (length(markers) >= 6) {p2 <- FeaturePlot(seurat, features = markers, pt.size = .1, ncol = 6)} else {p2 <- FeaturePlot(seurat, features = markers, pt.size = .1)}
  p3 <- DotPlot(seurat, features = markers, group.by = "seurat_clusters", cluster.idents = TRUE) + theme(axis.text.y = element_text(size = 8)) + coord_flip() + NoLegend()
  layout <- "
  ACC
  BBB
  BBB
  "
  p <- p1 + p2 + p3 + plot_layout(design = layout)
  ggsave(plot = p, filename = paste0(dir_anno, type, ".png"), height = 20, width = 30)
}


cat("
#######################################################################################################################\n
##  STEP 7: CREATING SUMMARY STATISTICS  ##############################################################################\n
#######################################################################################################################"
)

# Summary statistics

summary = data.frame(
  "Input_file" = data$object_path,
  "Batch" = data$batch_var,
  "Batch_run" = data$batch_run,
  "QC_features_min" = data$QC_feature_min,
  "QC_mito_max" = data$QC_mt_max,
  "Variable_features" = data$features_var,
  "PCA_dimensions" = data$pca_dims,
  "Amount_genes" = nrow(seurat),
  "Genes_detected_per_cell" = median(seurat@meta.data$nFeature_RNA),
  "Cells_after_QC" = ncol(seurat),
  "Cluster_resolution" = data$cluster_resolution
)
seurat@misc <- list(summary)

cat("
#######################################################################################################################\n
##  STEP 8: SAVING RESULTS  ###########################################################################################\n
#######################################################################################################################"
)

# Save RDS and convert to h5ad with seuratdisk

Idents(seurat) <- seurat$annotation_summary
saveRDS(seurat, paste0(out_dir, "scProcessoR_out.rds"))

SaveH5Seurat(seurat, filename = paste0(out_dir, "scProcessoR_out.h5Seurat"), overwrite = TRUE)
Convert(paste0(out_dir, "scProcessoR_out.h5Seurat"), dest = "h5ad", overwrite = TRUE)

# Save meta data as a table
write.csv(, file = paste0(out_dir, "scProcessor_meta.csv")

# Save data.json
data <- toJSON(data)
write(data, data_json)


cat("
#######################################################################################################################\n
##  ALL DONE - ANNOTATE  ##############################################################################################\n
#######################################################################################################################"
)