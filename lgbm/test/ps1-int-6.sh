#!/bin/bash
#SBATCH -A p31621
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --export=NONE
#SBATCH -J test-ps1-int-6
#SBATCH --output=test-ps1-int-6.out
#SBATCH --error=test-ps1-int-6.err
#SBATCH --mem=8G

module purge all

module load R/4.1.0

Rscript lgbm/test/ps1-int-6.R