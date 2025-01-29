#!/bin/bash
#SBATCH -A p31621
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --export=NONE
#SBATCH -J test-basic-1
#SBATCH --output=test-basic-1.out
#SBATCH --error=test-basic-1.err
#SBATCH --mem=8G

module purge all

module load R/4.1.0

Rscript lgbm/test/basic-1.R