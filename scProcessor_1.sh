#!/bin/bash
SCRIPT='/gpfs01/home/glanl/scripts/IMMUcan/scProcessor_1.R'
BATCH='none' #Fill in
GENES='250' #Adapt
MITO='15' #Adapt
PCA='30' #Adapt

ml R/3.6.3-foss-2016b

srun --mem=100G --time=4:00:00 --cpus-per-task=5 --partition=bigmem Rscript ${SCRIPT} ${BATCH} ${GENES} ${MITO} ${PCA}
