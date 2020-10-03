# IMMUcan

Processing scripts for scRNA-seq database

## Install instructions

- Follow install instructions for sceasy (https://github.com/cellgeni/sceasy)
- Install following R packages
```
install.packages(c("Seurat", "Harmony", "tidyverse", "readxl", "patchwork", "devtools"))
devtools::install_github("mahmoudibrahim/genesorteR") 
```

## Run scProcessor
Start from a seurat object that is named raw.rds

Run the following scProcessor steps **(FILL SQUARED BRACKETS)**

#### 1. Check Seurat object with check_seurat.R

``` 
Rscript check_seurat.R [PATH TO SEURAT OBJECT] 
```

Prints cell_id, gene_id, checks if metadata is correctly linked, checks normalization and stores it in data.rds and prints meta.data

#### 2. Test scProcessor to put desired QC thresholds and check batch variable

``` 
Rscript scProcessor_test.R
```

Saves QC and elbowplot in temp/

### 3. Run scProcessor_1

Bayer only: `bash scProcessor_1.sh` with slurm if on an HPC (adapt vars in bash file)

```
Rscript scProcessor_1.R [BATCH] [MIN FEATURES] [MAX MITO] [PCA DIMS]
```

Does QC, integration, dimensionality reduction, clustering and outputs marker gene plots in temp/

### 4. Annotate data

Check plots in temp folder
Create an excel file called annotation.xlsx with two columns: seurat clusters and abbreviation
Add the cell type abbreviation as defined in cell_ontology.xlsx in the abbreviation column to the corresponding seurat cluster

### 5. Run scProcessor_2

Bayer only: `bash scProcessor_2.sh` with slurm if on an HPC (adapt vars in bash file)

```
Rscript scProcessor_2.R
```

Links annotation to seurat_clusters, includes cell ontology, performs differential expression, creates output files
