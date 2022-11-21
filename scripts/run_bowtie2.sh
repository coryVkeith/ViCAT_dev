#!/bin/sh

# Your job will use 1 node, 28 cores, and 168gb of memory total.
#SBATCH -N 1
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=6G
#SBATCH -t 18:00:00

## calls the profile variable to pull sample names from a list iteratively
export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $PROFILE | tail -n 1`

#___________________________________________________________________________|
###                                                                         |
###      This changes the file suffix to match the data from the $PROFILE.  |
###                   Change the suffix to match the reads.                 |
#___________________________________________________________________________|

F1="${SMPLE}_R1.fastq"
R1="${SMPLE}_R2.fastq"

###bowtie2 build index of filtered references above
if [[ -f "$OUT_DIR/$SMPLE/blast_ref_fasta/${SMPLE}_nr_ref.fasta" ]]; then
    $BOWTIE/bowtie2-build -f $OUT_DIR/$SMPLE/blast_ref_fasta/${SMPLE}_nr_ref.fasta $OUT_DIR/$SMPLE/blast_ref_fasta/bowtie2index/${SMPLE}_nr_idx
    else
    echo "Error in run_bowtie2buildnmap.sh, bowtie2. No CD-HIT clustered reference fasta file found to build index. Stopping at build." >> $ERRORS/${SMPLE}_error.log
fi
##bowtie2 reference mapping to indexes
if [[ -f "$OUT_DIR/$SMPLE/blast_ref_fasta/${SMPLE}_nr_ref.fasta" ]]; then    
    $BOWTIE/bowtie2 --un-conc ${OUT_DIR}/${SMPLE}/bowtie2/unused_reads/"${SMPLE}_vircom_rm.fastq" -x ${OUT_DIR}/${SMPLE}/blast_ref_fasta/bowtie2index/"${SMPLE}_nr_idx" -q -1 ${RAW}/$F1 -q -2 ${RAW}/$R1 -S ${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr.SAM
    else
    echo "Error in run_bowtie2buildnmap.sh, bowtie2. No reference fasta file found to build index. Stopping at mapping" >> $ERRORS/${SMPLE}_error.log
fi
#Samtools convert to sorted bam, extract consensus, and remove sam
if [[ -f "${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr.SAM" ]]; then
    samtools view -bo ${OUT_DIR}/${SMPLE}/bowtie2/alignments/"${SMPLE}_vircom_nr.bam" ${OUT_DIR}/${SMPLE}/bowtie2/alignments/"${SMPLE}_vircom_nr.SAM"
    samtools sort ${OUT_DIR}/${SMPLE}/bowtie2/alignments/"${SMPLE}_vircom_nr.bam"  -o ${OUT_DIR}/${SMPLE}/bowtie2/alignments/"${SMPLE}_vircom_nr_sorted.bam"
    samtools mpileup --max-depth 0 -uf ${OUT_DIR}/${SMPLE}/blast_ref_fasta/${SMPLE}_nr_ref.fasta ${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${OUT_DIR}/${SMPLE}/bowtie2/consensus/${SMPLE}_vircom_nr_consensus.fasta

#_____________________________________________________________________________________________________________________________________|    
###    This code will copy the consensus sequences to a common folder and remove the .sam and unsorted .bam files to help save space. |
###                                     Comment this code if you want to keep intermediate files.                                     |
###                                                                                                                                   |
#_____________________________________________________________________________________________________________________________________|    
    cp ${OUT_DIR}/${SMPLE}/bowtie2/consensus/${SMPLE}_vircom_nr_consensus.fasta "$ALLCONSENSUS"
    rm ${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr.SAM
    rm ${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr.bam
    else
    echo "Error in run_bowtie2buildnmap.sh, samtools. No sam alignment file found. Stopping at bam file conversion." >> $ERRORS/${SMPLE}_error.log
fi

