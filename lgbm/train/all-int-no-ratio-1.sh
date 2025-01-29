#!/bin/bash
#SBATCH -A p31621
#SBATCH -p normal
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --export=NONE
#SBATCH -J train-all-int-no-ratio-1
#SBATCH --output=train-all-int-no-ratio-1.out
#SBATCH --error=train-all-int-no-ratio-1.err
#SBATCH --mem=32G

module purge all

module load R/4.1.0

Rscript lgbm/train/all-int-no-ratio-1.R