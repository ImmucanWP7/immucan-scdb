#!/bin/bash
SCRIPT='/gpfs01/home/glanl/scripts/scProcessor_1.R'
BATCH='sample'
NORM='FALSE'
GENES='250'
MITO='10'
PCA='20'

ml R/3.6.3-foss-2016b

srun --mem=100G --time=4:00:00 --cpus-per-task=5 --partition=bigmem Rscript ${SCRIPT} ${BATCH} ${NORM} ${GENES} ${MITO} ${PCA}
