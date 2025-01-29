#!/bin/bash
#SBATCH -A p31621
#SBATCH -p genhimem
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --export=NONE
#SBATCH -J data-split

module purge all

module load R/4.1.0

Rscript data-split.R