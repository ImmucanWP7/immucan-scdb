# scProcessoR
scProcessoR is an automated analysis pipeline for scRNA-seq data.  
Data can sequentially be processed to undergo quality control (QC) and cell types are automatically annotad by atlas mapping.


Data processing is done in two main steps including **quality controc (QC)** and **annotation**, by performing following steps:

- Quality control by thresholds for number of detected overall and mitochondiral genes detected
- Doublet detection [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)[^doubletfinder]
  
- Measure and correct batch effect through [Harmony](https://github.com/immunogenomics/harmony)[^harmony]
- Atlas mapping and annotation through [Symphony](https://github.com/immunogenomics/symphony)[^symphony]
- Identify tumor cells with CNA calling through [copyKAT](https://github.com/navinlabcode/copykat)[^copykat]
- Link everything to Bayer Cell ontology (in development)
- Differential expression and fraction analysis (in development)
  
- Universal output files ([SeuratDisk](https://github.com/mojaveazure/seurat-disk)), including:
    - .h5ad (compatible with [scanpy](https://scanpy.readthedocs.io/en/stable/)[^scanpy]
    - .h5adSeurat (compatible with [Seurat Tools](https://satijalab.org/seurat/)[^seurat])
    - .rds (compatible with [Seurat Tools](https://satijalab.org/seurat/))
    - .csv universal meta data file
  

## Installation instructions

Install [singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html).  Detailed instructions can be found on the [singularity website](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).  
Afterwards a singularity image of the workflow can be downloaded from here [sftp:xxx]. 

## Initial setup

For the analysis of scRNA-seq data a Seurat object, as well as a parameter file (data.json) needs to be provided.  
Detailed instructions on how to generate a Seurat object from various sources can be found in this [vignette]().  
The parameter file can be [downloaded](../scProcessoR_v2/ressources/data.json) or it will be automatically generated when first executing the pipeline.
All run parameters are specified in this single data.json file with some mandatory settings. Other parameters can be changed by the user or will otherwise be set automatically by the pipeline. See below for detailed description of the parameters.

| Parameter          | Mandatory | Description | Example / Preset |
| ------------------ |-----------| ----------- | ------- |
| object_path        | yes       | .rds-SeuratObject path, inside the <mounted path>       | ./seurat.rds      |
| batch_var          | yes       | Batch variable to integrate on                          | patient           |
| batch_run          | yes       | Batch variable that defines different experimental runs | sample            |
| malignant          | yes       | Does the data contain malignant cells [boolean]         | false             |
|                    |           |                                                         |                   |
| annotation         | no        | meta data columns that contain cell annotation          | ["author_annotation", "celltype"] |     
| metadata           | no        | meta data columns that will be kept in the final object | ["biopsy", "treatment"] |
| genes              | no        | genes of interest for which additional output will be generated | ["CD3D", "EPCAM"]
| norm               | no        | wether data is already normalised [boolean]. Will otherwise be determined automaticaly | false |
| QC_feature_min     | no        | Lower threshold for number of detected genes. Cells with lower gene number will be discarded | 250 | 
| QC_mt_max          | no        | Upper threshold for percentage of detected mitochondrial genes. Cells with higher percentage will be discarded | 15 |
| pca_dims           | no        | number of PCA dimensions to take for further processing | 30 |
| features_var       | no        | number of highly variable features to take for further processing | 2000 |
| cluster_resolution | no        | a sequence of different cluster resolutions, scProcessor will select the most optimal resolution | 0.3,0.5,1 |
| integrate          | no        | wether data should be integrated. This overwrites the pipelines suggestion | true |
|                    |           |  |   |
| normal_cells       | no        | cell type taken as normal cells to increase confidence of malingant cell prediction. null (standard Macrophages are taken), false (no normal cells taken), ["fibroblasts", "macrophages"] | null |
| nSample            | no        | number of cells to take for subsasmpling during test runs or intense computing steps | 50000 | 
| verbose            | no        | if true pipeline increases verbosity | false |
| test               | no        | if true data is subsampled to increase processing speed | false | 
    

## Run scProcessoR

The Seurat Object and the data.json file should be placed together in one folder.  
The pipeline can be run by executing the **QC** and **Annotate** apps sequentially. For a detailed walkthrough of the execution and generated outputs see this [vignette]().

### 1. Run **QC**:
```
singularity run --bind <YOUR-OUTPUT-DIR>:/mnt/outdir --app QC scProcessoR.sif
```

Afterwards check the data.json file, which might have been updated with the pipelines suggestions (A copy of your original file is kept in the output directory). Also inspect the QC plots generated in <YOUR-OUTPUT-DIR>/QC/.  

 ### 2. Run **Annotate**:
```
singularity run --bind <YOUR-OUTPUT-DIR>:/mnt/outdir --app Annotate scProcessoR.sif
```

Further help for individual apps can be accessed via ```singularity run-help --app <AppName> scProcessoR.sif```




[^seurat]: [Seurat Tools](https://satijalab.org/seurat/) ; DOI: [10.1038/nbt.3192](https://doi.org/10.1038/nbt.3192)
[^doubletfinder]: [DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors](https://doi.org/10.1016/j.cels.2019.03.003)
[^harmony]: [Fast, sensitive and accurate integration of single-cell data with Harmony](https://www.nature.com/articles/s41592-019-0619-0)
[^symphony]: [Efficient and precise single-cell reference atlas mapping with Symphony](https://www.nature.com/articles/s41467-021-25957-x)
[^copykat]: [Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes](https://doi.org/10.1038/s41587-020-00795-2)
[^scanpy]: [SCANPY: large-scale single-cell gene expression data analysis](https://doi.org/10.1186/s13059-017-1382-0)
[^NCBIID]: [NCBI ID](https://www.ncbi.nlm.nih.gov/taxonomy)
[^ncit]: [NCIT ontology tissue ID](https://www.ebi.ac.uk/ols/index)
[^mesh]: [MeSH disease IDs](https://www.ncbi.nlm.nih.gov/mesh?Db=mesh&Cmd=DetailsSearch&Term=%22Disease%22%5BMeSH+Terms%5D) 

