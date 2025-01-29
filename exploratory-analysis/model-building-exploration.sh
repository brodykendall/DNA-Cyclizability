#!/bin/bash
#SBATCH -A p31621
#SBATCH -p short
#SBATCH -t 0:30:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --export=NONE
#SBATCH -J model-building-exploration

module purge all

module load R/4.1.0

Rscript Model-building-exploration.R