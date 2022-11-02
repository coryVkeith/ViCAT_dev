#!/bin/sh

# Your job will use 1 node, 28 cores, and 168gb of memory total.
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 1:00:00

cd ${OUT_DIR}/relative_abundance

module load R
Rscript $SCRIPT_DIR/metadata.R
