#!/bin/bash
#SBATCH -A p31621
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --export=NONE
#SBATCH -J train-ps1-int-4
#SBATCH --output=train-ps1-int-4.out
#SBATCH --error=train-ps1-int-4.err
#SBATCH --mem=32G

module purge all

module load R/4.1.0

Rscript lgbm/train/ps1-int-4.R