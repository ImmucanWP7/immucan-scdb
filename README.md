# IMMUcan

Processing scripts for scRNA-seq database

Start from a seurat object

Run the following scProcessor steps (fill in squared brackets)

1. Check Seurat object with check_seurat.R
``` 
Rscript check_seurat.R [path to seurat object] 
```

2. Test scProcessor to put desired QC thresholds and batch variable (batch variable check follows in next update)
``` 
Rscript scProcessor_test.R [batch] [normalized] 
```

3. Run scProcessor_1
better done with slurm if on an HPC (adapt vars in bash file)
```
bash scProcessor_1.sh
```

on terminal
```
Rscript scProcessor_1.R [batch] [normalized] [min features] [max mito] [PCA dims]
```

4. Annotate data
Check plots in temp folder
Create an excel file called annotation.xlsx with two columns: seurat clusters and abbreviation
Add the cell type abbreviation as defined in cell_ontology.xlsx in the abbreviation column to the corresponding seurat cluster

5. Run scProcessor_2
better done with slurm if on an HPC (adapt vars in bash file)
```
bash scProcessor_2.sh
```

on terminal
```
Rscript scProcessor_2.R
```
