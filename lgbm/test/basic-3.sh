#!/bin/bash
#SBATCH -A p31621
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --export=NONE
#SBATCH -J test-basic-3
#SBATCH --output=test-basic-3.out
#SBATCH --error=test-basic-3.err
#SBATCH --mem=8G

module purge all

module load R/4.1.0

Rscript lgbm/test/basic-3.R