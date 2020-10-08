# IMMUcan

Processing scripts for scRNA-seq database

## Install instructions

- Follow install instructions for sceasy (https://github.com/cellgeni/sceasy)
- Install following R packages
- Get CHETAH_reference_updatedAnnotation.RData from IMMUcan teams channel
```
install.packages(c("Seurat", "Harmony", "tidyverse", "readxl", "patchwork", "devtools"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CHETAH")
devtools::install_github("mahmoudibrahim/genesorteR") 
```

## Run scProcessor
- Start from a seurat object
- Run the following steps 
- scProcessor creates two folders: temp/ and out/
- **FILL SQUARED BRACKETS**

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
- [BATCH]: FULL name of the batch column in the metadata (often patient or sample)
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
