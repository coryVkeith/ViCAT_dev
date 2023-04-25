#!/bin/sh

# Your job will use 1 node, 28 cores, and 168gb of memory total.
#SBATCH -N 1
#SBATCH -n 28
#SBATCH --mem=168G
#SBATCH -t 48:00:00

#Load ViCAT environment
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate ViCAT

### sample names from a list correlates to raw data
export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $PROFILE | tail -n 1`

# Copies blast hit fastas to relative_abundance dir
if [ -f "${OUT_DIR}/${SMPLE}/blast_ref_fasta/${SMPLE}_blasthits.fasta" ]; then
    cp "${OUT_DIR}/${SMPLE}/blast_ref_fasta/${SMPLE}_blasthits.fasta" "${OUT_DIR}/relative_abundance"
fi

cd ${OUT_DIR}/relative_abundance/

# Grabs sequence description from the line, turns first space to a comma, removes the ">", and turns the first comma to a tab for the resulting tsv file.
if [ -f "${SMPLE}_blasthits.fasta" ]; then
    for f in ${SMPLE}*
        do
            grep ">" $f | sed 's/\ /,/' | sed 's/^.//' >> seqdescription.tsv
    done
    else 
    echo " No virus references in directory ${OUT_DIR}/relative_abundance for ${SMPLE}." >> $ERRORS/${SMPLE}_error.log
fi
# Deduplicates the files and writes taxa matrix
sort -u seqdescription.tsv -o taxmat.csv

# Adds headers to the taxa matrix
sed -i 1i"RefID,Species" taxmat.csv 

if [ -f "${SMPLE}_blasthits.fasta" ]; then
    rm *blasthits.fasta 
fi
