#!/bin/bash
#SBATCH -A p31621
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --export=NONE
#SBATCH -J model2-network-int

module purge all

module load R/4.1.0

Rscript 10-lgbm-model2-network-int.R