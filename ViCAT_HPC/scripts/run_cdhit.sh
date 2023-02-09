#!/bin/sh

# Your job will use 1 node, 28 cores, and 168gb of memory total.
#SBATCH -N 1
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=6G
#SBATCH -t 18:00:00

source activate cdhit
## calls the profile variable to pull sample names from a list iteratively
export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $PROFILE | tail -n 1`

#_______________________________________________________________________________________________________________________|
###                                                                                                                     |
###CD-HIT remove sequences with defined similarity of 88.5% from the list of blast sequences for bowtie2 index building |
###                            Change -c to the prefered percent similarity if desired                                  |
#_______________________________________________________________________________________________________________________|

if [[ -f "${OUT_DIR}/${SMPLE}/blast_ref_fasta/"${SMPLE}_blasthits.fasta"" ]]; then
    cd-hit-est -i ${OUT_DIR}/${SMPLE}/blast_ref_fasta/"${SMPLE}_blasthits.fasta" -o ${OUT_DIR}/${SMPLE}/blast_ref_fasta/"${SMPLE}_nr_ref.fasta" -c 0.88.5
    else
    echo "Error in run_cdhit. No blast file found. Did not cluster sequences." >> $ERRORS/${SMPLE}_error.log
fi
