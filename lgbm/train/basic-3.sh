#!/bin/bash
#SBATCH -A p31621
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --export=NONE
#SBATCH -J train-basic-3
#SBATCH --output=train-basic-3.out
#SBATCH --error=train-basic-3.err
#SBATCH --mem=16G

module purge all

module load R/4.1.0

Rscript lgbm/train/basic-3.R