# IMMUcan

Processing scripts for scRNA-seq database

## Install instructions

- Follow install instructions for sceasy (https://github.com/cellgeni/sceasy)
- Get CHETAH_reference_updatedAnnotation.RData from IMMUcan teams channel
- Install following R packages
```
install.packages(c("Seurat", "tidyverse", "readxl", "patchwork", "devtools", "data.table", "BiocManager", "remotes", "WriteXLS", "pheatmap", "plyr"))
BiocManager::install(c("CHETAH", "SingleCellExperiment"))
devtools::install_github("mahmoudibrahim/genesorteR") 
devtools::install_github("immunogenomics/harmony")
remotes::install_github("mojaveazure/seurat-disk")
```

On terminal
```
git clone https://github.com/soumelis-lab/IMMUcan.git
```

## Before starting

Change the paths to files provided in the script
- cellMarker_path = PATH to TME_markerGenes.xlsx
- chetahClassifier_path = PATH to CHETAH_reference_updatedAnnotation.RData
- cellOntology_path = PATH to cell_ontology.xlsx

## Run scProcessor
- Start from a seurat object
- scProcessor creates two folders where it stores files: temp/ and out/
- **FILL IN THE SQUARED BRACKETS**

### 1. Check Seurat object

Input:
- [SEURAT]: path to seurat object

``` 
Rscript check_seurat.R [SEURAT] 
```

Output:
- Prints cell_id, gene_id
- Checks if metadata is correctly linked
- Checks normalization: raw or normalized (stores it in data.rds)
- Prints meta.data slot

### 2. Select desired QC thresholds and batch variable

Input:
- [BATCH]: fill in part or full name of batch variable(s) e.g. "patient|sample|plate"

``` 
Rscript scProcessor_test.R [BATCH]
```

Output:
- QC plot and elbowplot in temp/

### 3. scProcessor_1

Bayer only: `bash scProcessor_1.sh` with slurm if on an HPC (adapt vars in bash file)

Input:
- [BATCH]: FULL name of the batch column in the metadata (often patient or sample), if no integration is desired use "none"
- [FEATURES]: minimum amount of detected genes per cell allowed
- [MITO]: maximum percentage of mitochondrial reads per cell allowed
- [PCA]: number of PCA dimensions to use for downstream processing

```
Rscript scProcessor_1.R [BATCH] [FEATURES] [MITO] [PCA]
```

Steps:
- QC
- Integration
- Dimensionality reduction + clustering
- Supervised classification

Output:
- QC plot in out/
- Integration plot in out/
- Supervised classification in out/
- Seurat clustering + gene modules in temp/
- Comparison seurat clusters and supervised classification in temp/
- annotation.xls in out/

### 4. Annotate clusters

- Check plots in temp/:
  - marker gene plots
  - dotplot
  - top10 diffentially expressed genes per seurat cluster
- In out/annotation.xls, fill in cell types as defined in cell_ontology.xlsx in the abbreviation column


### 5. scProcessor_2

Bayer only: `bash scProcessor_2.sh` with slurm if on an HPC (adapt vars in bash file)

```
Rscript scProcessor_2.R
```
Steps:
- Links annotation to seurat_clusters
- Includes cell ontology
- Differential expression

Output: 
- AverageExpression matrices and DE_results per annotation level in out/
- Metadata.tsv in out/
- cellCount.tsv in out/
- cellxgene.h5ad in out/
