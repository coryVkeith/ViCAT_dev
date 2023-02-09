#!/bin/sh

#Your job will use 1 node, 28 cores, and 64gb of memory total.

#SBATCH -N 1
#SBATCH -n 28
#SBATCH --mem=64gb
#SBATCH -t 1:00:00

# Activates anaconda environment for this code to run.
source activate bio

#Calls the correct sample ID from the profile.txt
export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $PROFILE | tail -n 1`

# These are the variable files for the python tool, gather_contigs.py
if [[ -f "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta" ]]; then
    export F1="$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta"
    else
    export F1="$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_contigs.fasta"
fi
export F2="$OUT_DIR/$SMPLE/blast/${SMPLE}_unique_contigs.tsv"
export F3="$OUT_DIR/virus_contigs/${SMPLE}_contigs_blasthits.fasta"

#python code for the gather_contigs.py tool.

if [[ -f "$F1" ]]; then
    python $SCRIPT_DIR/gather_contigs.py $F1 $F2 $F3
    else
    echo " Error in run_gather_contigs.sh. No contigs or scaffolds found." >> $ERRORS/${SMPLE}_error.log
fi

