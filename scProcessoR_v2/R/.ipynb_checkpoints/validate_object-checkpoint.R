import("Seurat")

export("validateSeurat")


validateSeurat <- function(seurat, data){
    
  if (!data$batch_var %in% colnames(seurat@meta.data)) {
    stop("Batch variable not in metadata")
  }
    
  if (!data$batch_run %in% colnames(seurat@meta.data)) {
    stop("Batch run variable not in metadata")
  }

  if (length(data$batch_run) != 1) {
    stop("More than one batch run variable is not allowed")
  }
    
  # Check if Cell IDs are correctly linked
  if (sum(colnames(seurat) == rownames(seurat@meta.data)) != ncol(seurat)) {
    stop("Cell IDs linked uncorrectly")
  }
      
  

  # Check if human or mouse data
  check_genes <- toupper(rownames(seurat))
  if (!identical(check_genes, rownames(seurat)) & data$species == "human") {
    message(paste(rownames(seurat)[!identical(rownames(seurat)), check_genes], collapse = " /n/n"))
    message("Only mouse and human gene annotation are supported in this version of scProcessor /n/n")
    stop("the genes above seem to not follow the convention for human gene annotation, 
         please check or specify 'mouse' in the data.json file if you are processing mouse data")
  }

  if (!identical(check_genes, rownames(seurat)) & data$species != "mouse") {
    message(paste(rownames(seurat)[identical(rownames(seurat)), check_genes], collapse = " /n/n"))
    message("Only mouse and human gene annotation are supported in this version of scProcessor /n/n")
    stop("the genes above seem to not follow the convention for mouse gene annotation, 
         please check or specify 'human' in the data.json file if you are processing mouse data")
  }
    
  if (!any(data$species %in% c("mouse", "human"))) {
    stop("Only mouse and human gene annotation are supported in this version of scProcessor /n/n")
  }

  return(TRUE)
    
}



