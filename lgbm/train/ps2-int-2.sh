#!/bin/bash
#SBATCH -A p31621
#SBATCH -p normal
#SBATCH -t 16:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --export=NONE
#SBATCH -J train-ps2-int-2
#SBATCH --output=train-ps2-int-2.out
#SBATCH --error=train-ps2-int-2.err
#SBATCH --mem=64G

module purge all

module load R/4.1.0

Rscript lgbm/train/ps2-int-2.R