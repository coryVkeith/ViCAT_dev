#!/bin/sh

# Your job will use 1 node, 28 cores, and 5gb of memory total.
#SBATCH -N 1
#SBATCH -n 28
#SBATCH --mem=5G
#SBATCH -t 1:00:00

## this code writes accesion numbers to a list

export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $PROFILE | tail -n 1`

if [[ -f "$OUT_DIR/$SMPLE/blast/${SMPLE}_blastn_results.tsv" ]]; then
    f="$OUT_DIR/$SMPLE/blast/${SMPLE}_blastn_results.tsv"; else
        echo "Error in run_collect_contigs_accessions.sh. No blast results to write accessions. Stopping at writing accessions."
fi

#__________________________________________________________________________________________________________________________________|
###                                                                                                                                |
### filters hits that only have a match length of greater than 150 and less than 4500 and writes the accession number to a list.   |
###                                  Change the awk code in the first step of the pipe for different filters.                      |
#__________________________________________________________________________________________________________________________________|

echo $f >> $ERRORS/${SMPLE}_error.log
awk '$5>150 && $5<4500 {print $2}' $f | sed '$!N; /^\(.*\)\n\1$/!P; D'  >> "$OUT_DIR/$SMPLE/blast/${SMPLE}_unique_accessions.tsv"
# Takes just the Accession number from the line and writes to a new list
for l in $OUT_DIR/$SMPLE/blast/${SMPLE}_unique_accessions.tsv
    do
    while read line; do
        pre="${line%.*}"
        suff="${pre#*|}"
        printf "$suff\n" >> $OUT_DIR/$SMPLE/blast/trimmed_${SMPLE}_unique_accessions.tsv
    done <$l
done

#__________________________________________________________________________________________________________________________________________________________|
###                                                                                                                                                        |
# filters hits that only have a match length of greater than 300 and less than 4500 and then removes duplicates and then writes the contig name to a file. |
###                                  Change the awk code in the first step of the pipe for different filters.                                              |
#__________________________________________________________________________________________________________________________________________________________|

awk '$5>250 && $5<4500 {print $1}' $f | sed '$!N; /^\(.*\)\n\1$/!P; D'  >> "$OUT_DIR/$SMPLE/blast/${SMPLE}_unique_contigs.tsv"
# This code writes the fasta file from the unique accessions file from the blast hits step to a new fasta file for use in building bowtie indexes

if [[ -f "$OUT_DIR/$SMPLE/blast/trimmed_${SMPLE}_unique_accessions.tsv" ]]; then
    while read line; do
        grep ">$line" -A 1 $BLAST1line | cat >> $OUT_DIR/$SMPLE/blast_ref_fasta/${SMPLE}_blasthits.fasta
    done <$OUT_DIR/$SMPLE/blast/trimmed_${SMPLE}_unique_accessions.tsv
    else
    echo "Error in run_collect_contigs_accessions.sh. No accessions found in file. Stopping at collecting accessions." >> $ERRORS/${SMPLE}_error.log
fi
