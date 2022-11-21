#!/bin/sh

# Your job will use 1 node, 28 cores, and 168gb of memory total.
#SBATCH -N 1
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=6G

#SBATCH -t 48:00:00


SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $PROFILE | tail -n 1`

#___________________________________________________________________________|
###                                                                         |
###      This changes the file suffix to match the data from the $PROFILE.  |
###                   Change the suffix to match the reads.                 |
#___________________________________________________________________________|

F1="${SMPLE}_R1.fastq"
R1="${SMPLE}_R2.fastq"

echo $F1 $R1 "$RAW/$F1" "$RAW/$R1"
$SPADES/rnaviralspades.py -t 28 -o "${OUT_DIR}/${SMPLE}/spades" -1 ${RAW}/"$F1" -2 ${RAW}/"$R1"
#-----------------------------------------------------------------------------------------------------------|
###                                                                                                         |
### This code filters the spades contigs by 30x coverage and 250bp - 5000bp lengths for faster blast times. |
###                        To change replade the numbers at the awk step in the pipe.                       |
###                                                                                                         |
#-----------------------------------------------------------------------------------------------------------|
c="$OUT_DIR/$SMPLE/spades/contigs.fasta"
s="$OUT_DIR/$SMPLE/spades/scaffolds.fasta"

if [[ -f "$s" ]]; then
    $SEQKIT fx2tab "$s" | $CSVTK mutate -H -t -f 1 -p "cov_(.+)" | $CSVTK mutate -H -t -f 1 -p "length_([0-9]+)" | awk -F "\t" '$4>=30 && $5>=250 && $5<10000' | $SEQKIT tab2fx >> "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta";elif [[ -f "$c" ]]; then
        $SEQKIT fx2tab "$c" | $CSVTK mutate -H -t -f 1 -p "cov_(.+)" | $CSVTK mutate -H -t -f 1 -p "length_([0-9]+)" | awk -F "\t" '$4>=30 && $5>=250 && $5<10000' | $SEQKIT tab2fx >> "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_contigs.fasta"; else
            echo "Error in run_spades.sh. No contigs were assembled from the sample. Check the sample file in the raw data." >> $ERRORS/${SMPLE}_error.log
fi

# Copies the filtered spades files to a common folder
if [[ -f "$s" ]]; then
    cp "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta" "$ALLSPADES"; elif [[ -f "$c" ]]; then
        cp "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta" "$ALLSPADES"; else
            echo "No files."
fi
