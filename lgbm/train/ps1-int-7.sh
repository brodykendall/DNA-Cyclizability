#!/bin/bash
#SBATCH -A p31621
#SBATCH -p normal
#SBATCH -t 16:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --export=NONE
#SBATCH -J train-ps1-int-7
#SBATCH --output=train-ps1-int-7.out
#SBATCH --error=train-ps1-int-7.err
#SBATCH --mem=128G

module purge all

module load R/4.1.0

Rscript lgbm/train/ps1-int-7.R