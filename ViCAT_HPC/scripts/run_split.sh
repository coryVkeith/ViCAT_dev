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

cd ${OUT_DIR}/${SMPLE}/bowtie2/alignments/
echo ${OUT_DIR}/${SMPLE}/blast_ref_fasta/"${SMPLE}_blasthits.fasta"
$BAMTOOLS split -in ${OUT_DIR}/${SMPLE}/bowtie2/alignments/"${SMPLE}_vircom_nr_sorted.bam" -reference

#_________________________________________________________________________________|
###                                                                               |
###               This code provides the stats for the read mapping.              |
###    Change the file name to match the regular expression of the BLASTdb used.  |
###                    Also change the Suff and OUT variable.                     |
##________________________________________________________________________________|
 
for f in ${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr_sorted.REF_NC*.bam
    do
    Suff=${f%%.bam}
    OUT=${Suff##*/}
    samtools stats $f >> ${OUT_DIR}/${SMPLE}/bowtie2/alignments/stats/"${OUT}_stats.txt"
done

#####This code takes just the number of reads that map from the stats file generated above, and the name of the sample.
for s in ${OUT_DIR}/${SMPLE}/bowtie2/alignments/stats/${SMPLE}_*stats.txt
    do
    Suff=${s##*/}
    HIT="${Suff%%_stats.txt}"
    awk 'NR==14' $s | sed 's/ \+ /\t/g' | awk '{ print $4 }' | xargs printf "$HIT \t %s\n" >> ${OUT_DIR}/${SMPLE}/bowtie2/alignments/stats/"${SMPLE}_readscount.tsv"
done

##### This code combines the individual reads counts into the combined read counts
for f in ${OUT_DIR}/${SMPLE}/bowtie2/alignments/stats/*_readscount.tsv
    do
    while read line
        do
        printf "$line\n" >> ${OUT_DIR}/relative_abundance/combined_readscount.tsv
    done <$f
done
