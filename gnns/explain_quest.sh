#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguest
#SBATCH -N 1
#SBATCH --mem=180G
#SBATCH -t 42:00:00
#SBATCH --job-name="Explain Tiling"

module purge all

cd /home/cbk6686

python Documents/DNA_Cyclizability/benchmarks/gnns/explain_quest.py
