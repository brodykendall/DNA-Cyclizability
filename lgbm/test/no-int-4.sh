#!/bin/bash
#SBATCH -A p31621
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --export=NONE
#SBATCH -J test-no-int-4
#SBATCH --output=test-no-int-4.out
#SBATCH --error=test-no-int-4.err
#SBATCH --mem=8G

module purge all

module load R/4.1.0

Rscript lgbm/test/no-int-4.R