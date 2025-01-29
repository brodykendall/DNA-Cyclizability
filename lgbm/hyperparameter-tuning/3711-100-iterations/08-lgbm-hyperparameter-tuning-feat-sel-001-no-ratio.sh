#!/bin/bash
#SBATCH -A p31621
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --array=1-20
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --export=NONE
#SBATCH -J hyperparameter-tuning-no-ratio

module purge all

module load R/4.1.0

Rscript lgbm-hyperparameter-tuning-feat-sel-001-no-ratio.R ${SLURM_ARRAY_TASK_ID}

