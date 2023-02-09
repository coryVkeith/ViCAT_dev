#!/bin/sh

#Your job will use 1 node, 28 cores, and 168gb of memory total.

#SBATCH -N 1
#SBATCH -n 28
#SBATCH --mem=168gb
#SBATCH -t 48:00:00

module load blast


## This code checks if there are scaffolds or contigs and then sets the variables for use in the blast code below. If no files are present, an error log is printed.
export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $PROFILE | tail -n 1`
if [[ -f "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta" ]]; then
    export f="$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta"; elif [[ -f "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_contigs.fasta" ]]; then
        export f="$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_contigs.fasta"; else
            echo "Error in run_blastn.sh. There are not filtered contigs or scaffolds." >> $ERRORS/${SMPLE}_error.log
fi
export OUTF="$OUT_DIR/$SMPLE/blast/${SMPLE}_blastn_results.tsv"

echo $BLASTDB

#___________________________________________________________________________________________|
###                                                                                         |
### This is the blastn algorithm. To allow more dissimilar sequences: '-task dc-megablast'. | 
###          To change the max number of hits per query: '-max_target_seqs [#]'             |
###                                                                                         |
#___________________________________________________________________________________________|

if [[ -f "$f" ]]; then
    blastn -db $BLASTDB -query $f -out $OUTF -evalue 0.01 -max_target_seqs 2 -outfmt "6 qseqid sseqid ssciname pident length mismatch gapopen qstart qend sstart send evalue bitscore"; else
        echo "Error in run_blastn.sh. No contigs or scaffolds found to blast." >> $ERRORS/${SMPLE}_error.log
fi
# Copies blast files to a common folder
if [[ -f "$OUTF" ]]; then
    cp "$OUTF" "$ALLBLAST"
fi
