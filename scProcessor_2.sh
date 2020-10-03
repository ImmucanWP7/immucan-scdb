#!/bin/bash
SCRIPT='/gpfs01/home/glanl/scripts/IMMUcan/scProcessor_2.R'

module use /gpfs01/sw/eb-2019/modules/all /gpfs01/sw/eb-rh7/modules/all /gpfs01$
ml Anaconda3/2020.02.lua
source activate	sceasy
ml R/3.6.3-foss-2016b

srun --mem=100G --time=6:00:00 --cpus-per-task=5 --partition=bigmem Rscript ${SCRIPT}

