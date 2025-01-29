#!/bin/bash
#SBATCH -A p31621
#SBATCH -p normal
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --export=NONE
#SBATCH -J train-ps1-int-no-ratio-1
#SBATCH --output=train-ps1-int-no-ratio-1.out
#SBATCH --error=train-ps1-int-no-ratio-1.err
#SBATCH --mem=32G

module purge all

module load R/4.1.0

Rscript lgbm/train/ps1-int-no-ratio-1.R