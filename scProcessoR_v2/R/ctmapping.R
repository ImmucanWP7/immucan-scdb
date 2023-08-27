import("stats", "setNames")
import("tidyr")
import("utils")
import("Seurat")
import("symphony")
import("ggplot2")

export("mapCtLevels")


mapQuery <- function(
  se_query, # Seurat object to predict labels for
  batchvar_query, # Batch variable of the query data, regressed out with harmony
  path_symphony_reference_dir, # this is a path to a ready to use symphony reference as .rds and complementary uwot-model
  query_assay = "RNA",# Assay to use from the query object
  query_slot = "data", # Data slot to use from the query object; has to be "counts" or "data"
  annotation_col_ref, # reference annotation column; we want to fix this to our corresponding level names
  k = 5, # this might be adjusted depending on achieved mapping confidence
  save_dir = NULL, # if we want to save the output of this function or a csv file with cell type predictions this has to be set
  plot_dir = NULL,
  plt_sufix = NULL
){
  
  set.seed(42)
  
  # Get normalised expression counts from query
  exprs_norm <- Seurat::GetAssayData(
    se_query, 
    assay = query_assay, 
    slot = query_slot
  )
  
  # Get meta data from query
  meta_query <- Seurat::FetchData(se_query, vars = c(batchvar_query))
  meta_query <- meta_query %>%
    dplyr::mutate(across(everything(), as.character))
  
  # Map to reference
  stopifnot(query_slot %in% c("counts", "data"))
  do_normalize <- ifelse(query_slot == "data", FALSE, TRUE)
  ## read reference model
  sy_ref_file <- list.files(path = path_symphony_reference_dir, pattern = ".rds", full.names = T)
  
  stopifnot(length(sy_ref_file) == 1)
  sy_ref <- readRDS(sy_ref_file)
  
  sy_ref_model <- list.files(path = path_symphony_reference_dir, pattern = ".model", full.names = T)
  stopifnot(length(sy_ref_model) == 1)
  ### Fix uwot model path for reference
  sy_ref$save_uwot_path <- sy_ref_model
  
  message("     Mapping query to reference from directory: ", path_symphony_reference_dir, "\n", 
          "          with reference: ", sy_ref_file, "\n",
          "          and uwot-model: ", sy_ref$save_uwot_path)
  
  query = symphony::mapQuery(
    exp_query = exprs_norm,        # query gene expression (genes x cells)
    metadata_query = meta_query,                  # query metadata (cells x attributes)
    ref_obj = sy_ref,                      # Symphony reference object
    do_normalize = do_normalize, # perform log(CP10k+1) normalization on query
    vars = batchvar_query
  )           
  
  message("    Predicting knn")
  ## Predict query cell types using k-NN
  query = symphony::knnPredict(
    query_obj = query, 
    ref_obj = sy_ref,
    train_labels = sy_ref$meta_data[, annotation_col_ref],
    save_as = "cell_type_pred", 
    k = k,
    confidence = TRUE, 
    seed = 0
  )
  
  ## Add UMAP embeddings and comparison to reference
  sy_ref$meta_data$ref_query <- "reference"
  sy_ref$meta_data$cell_type_combined <- sy_ref$meta_data[, annotation_col_ref]
  sy_ref$meta_data$cell_type_pred_prob <- 0
  data_ref <- cbind(sy_ref$meta_data, sy_ref$umap$embedding)
  
  meta_query <- query$meta_data
  query$meta_data$cell_type_combined <- query$meta_data$cell_type_pred
  query$meta_data$ref_query <- "query"
  data_query <- cbind(query$meta_data, query$umap)
  
  meta_data_combined <- plyr::rbind.fill(data_ref, data_query)
  query$plot_data <- meta_data_combined
  
  ## Save query data object
  if(!is.null(save_dir)){
    message("Saving csv to directory:", save_dir)
    dir.create(save_dir, recursive = T)
    #saveRDS(object = query, file = paste0(save_dir, "/query_data.rds"))
    write.csv(meta_query, file = paste0(save_dir, "/ct_pred.csv"))
  }
  
  ## Save plots
  if(!is.null(plot_dir)){
    dir.create(plot_dir, recursive = T)
    if(!(is.null(plt_sufix))){
      plt_sufix <- paste0("_", plt_sufix)
    }
    p1 <- ggplot(query$plot_data, aes(UMAP1, UMAP2, color=cell_type_combined)) +
      geom_point() +
      theme_classic() +
      facet_grid(.~ref_query)
    ggsave(filename = paste0(plot_dir, "/UMAP_ct_pred", plt_sufix, ".png"), plot = p1)
    
    p2 <- ggplot(query$plot_data, aes(UMAP1, UMAP2, color=cell_type_pred_prob)) +
      geom_point() +
      theme_classic() +
      facet_grid(.~ref_query)
    ggsave(filename = paste0(plot_dir, "/UMAP_ct_pred_score", plt_sufix, ".png"), plot = p2)
    
    data_plot <- query$plot_data
    data_plot <- data_plot[data_plot$ref_query == "query", ]
    p3 <- ggplot(data_plot, 
                 aes(cell_type_pred, cell_type_pred_prob, color=cell_type_pred)) +
      geom_point(position = "jitter") +
      geom_violin() +
      theme_classic()
    ggsave(filename = paste0(plot_dir, "/VlnPlot_ct_pred_score", plt_sufix, ".png"), plot = p3)
  }
  
  ## Return data
  message("    ** Query mapping done")
  return(query)
  
}

mapCtLevels <- function(
  se_query,    # Query object as Seurat object
  path_ref_list,    # Names list for references for details see examples below; for level 1 "path_ref_l1" is sufficient
  batchvar_query = "Patient",    # Batch variable of the query data, regressed out with harmony 
  query_assay = "RNA",
  query_slot = "data",
  annotation_col_ref = "ct_id",
  k = 5,
  save_dir = NULL,
  save_all = FALSE #if set to TRUE, all mapping outputs on every level are saved, otherwise only final table
){
  
  #######################################################################################################################
  #######################################################################################################################
  ## Function for mapping query object to multiple references in multiple layers of annotation.
  ## Query has to be a Seurat object
  ## References have to be symphony objects.
  ## "path_ref_list" is a list that has to adhere to following example structure:
  ##
  ##
  ## For the reference list all cell types at level 1 are evaluated for deeper annotation. Mapping to deeper levels can be 
  ## skipped by providing "NULL" instead of path to symphony-reference directory.
  ## Name of reference in list has to follow "path_ref_l2_[lowercase cell annotation level1]".
  ##
  ### path_ref_list <- list(
  ###  path_ref_l1 = "./reference-l1"
  ###  , path_ref_l2_epithelial = "./reference-l2-epi"
  ###  , path_ref_l2_stromal = "./reference-l2-stromal"
  ###  , path_ref_l2_endothelial = "./reference-l2-endo"
  ###  , path_ref_l2_tcell = "./reference-l2-tcell"
  ###  , path_ref_l2_myeloid = "./reference-l2-myeloid"
  ###  , path_ref_l2_cycling = NULL
  ###  , path_ref_l2_mast = NULL
  ###  , path_ref_l2_bcell = NULL
  ###  , path_ref_l2_plasma = NULL)
  #######################################################################################################################
  
  # Prepare some variables
  plot_dir <- paste0(save_dir, "/plots")
  if(!(is.null(save_dir))){
    save_dir <- paste0(save_dir, "/data")
  }
  
  # First level mapping
  message("**** First level mapping *****")
  # Run mapping for level 1
  save_level_dir <- ifelse(
    all(save_all & !(is.null(save_dir))), 
    paste0(save_dir, "/level1"), 
    NULL
  )
  
  data_l1 <- mapQuery(
    se_query,
    batchvar_query = batchvar_query,
    path_symphony_reference_dir = path_ref_list[["path_ref_l1"]],
    query_assay = query_assay, 
    query_slot = query_slot,
    annotation_col_ref = annotation_col_ref, 
    k = k,
    save_dir = save_level_dir,
    plot_dir = plot_dir,
    plt_sufix = "level1"
  )
  
  message("**** DONE - First level mapping *****")
  
  # Fix meta data
  l1_ann <- data_l1$meta_data
  names(l1_ann)[names(l1_ann) == 'cell_type_pred'] <- 'cell_type_pred_l1'
  names(l1_ann)[names(l1_ann) == 'cell_type_pred_prob'] <- 'cell_type_pred_prob_l1'
  l1_ann$cell_type_combined <- NULL
  l1_ann$ref_query <- NULL
  
  # Second level mapping
  message("**** Second level mapping *****")
  
  # Save all cells in data set
  cells_all <- rownames(data_l1$meta_data)
  
  # Split up at first level
  l1_ann$cell <- rownames(l1_ann)
  
  l1_ann$cell_type_pred_l1 <- factor(l1_ann$cell_type_pred_l1, levels = unique(l1_ann$cell_type_pred_l1))
  l1_ann <- split(l1_ann, l1_ann$cell_type_pred_l1)
  
  ct_l1 <- setNames(as.list(names(l1_ann)), names(l1_ann))
  
  data_l2 <- lapply(ct_l1, function(ct){
    
    data_tmp <- l1_ann[[ct]]
    #Next level annotation is only performed if there are more than 10 cells in the level 1 group
    
    # Get variables
    print(ct)
    stopifnot(length(ct) == 1)
    cell_use <- data_tmp$cell
    
    message("***************** \n\n", "Cell type ", ct)
    
    # Set which reference to use
    ref_use <- grep(tolower(ct), names(path_ref_list), value = T)
    
    if(identical(ref_use, character(0))){
      ref_use <- NULL
    } else {
      ref_use <- path_ref_list[[ref_use]]
    }
    
    ncells_min <- 10
    if(length(data_tmp$cell_type_pred_l1) < ncells_min){
      
      message("!!! Less than ", ncells_min, " for ", ct, " - not processed")
      
      # Carry over cell type level 1 to level 2, since there is no reference for deeper annotation
      data_tmp$cell_type_pred_l2 <- data_tmp$cell_type_pred_l1
      data_tmp$cell_type_pred_prob_l2 <- NA
      
      print("** not enough cells")
      return(data_tmp)
      
    }
    
    if(is.null(ref_use)){
      
      message("!!! No reference found for cell type ", ct)
      
      # Carry over cell type level 1 to level 2, since there is no reference for deeper annotation
      data_tmp$cell_type_pred_l2 <- data_tmp$cell_type_pred_l1
      data_tmp$cell_type_pred_prob_l2 <- NA
      
      print("** no reference found")
      return(data_tmp)
      
    }
    
    print("** run main mapping")
    # Run level 2 atlas mapping
    message("#### Processing l2 - ", ct, " with reference in ", ref_use, " #### ")
    
    # Subset query object
    se_query_sub <- subset(se_query, cells = cell_use)
    
    save_level_dir <- ifelse(
      all(save_all & !(is.null(save_dir))), 
      paste0(save_dir, "/level2"), 
      NULL
    )
    
    data_l2_sub <- mapQuery(
      se_query_sub,
      batchvar_query = batchvar_query,
      path_symphony_reference_dir = ref_use,
      query_assay = query_assay, 
      query_slot = query_slot,
      annotation_col_ref = annotation_col_ref, 
      k = k,
      save_dir = save_level_dir,
      plot_dir = plot_dir,
      plt_sufix = paste0("level2-", ct)
    )
    
    # Prepare output meta data table
    data_tmp_l2 <- data_l2_sub$meta_data
    data_tmp_l2$cell <- rownames(data_tmp_l2)
    names(data_tmp_l2)[names(data_tmp_l2) == 'cell_type_pred'] <- 'cell_type_pred_l2'
    names(data_tmp_l2)[names(data_tmp_l2) == 'cell_type_pred_prob'] <- 'cell_type_pred_prob_l2'
    data_tmp_l2$cell_type_combined <- NULL 
    data_tmp_l2$ref_query <- NULL
    
    stopifnot(all(data_tmp_l2$cell %in% data_tmp$cell))
    data_tmp <- merge(data_tmp, data_tmp_l2, by = "cell")
    
    print("** Done main mapping")
    return(data_tmp)
    
  })
  
  # Make column names homogeneous
  colnames_needed <- c("cell_type_pred_l1", "cell_type_pred_prob_l1", 
                       "cell_type_pred_l2", "cell_type_pred_prob_l2", 
                       "cell")
  
  data_l2 <- lapply(data_l2, function(x){
    stopifnot(all(colnames_needed %in% colnames(x)))
    x <- x[, colnames_needed]
    
    return(x)
    
  })
  
  ann_data <- do.call("rbind", data_l2)
  stopifnot(all(ann_data$cell %in% cells_all))
  
  return(ann_data)
  
}
