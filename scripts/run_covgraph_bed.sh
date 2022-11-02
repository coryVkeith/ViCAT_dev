#!/bin/sh

# Your job will use 1 node, 28 cores, and 168gb of memory total.
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 1:00:00

SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $PROFILE | tail -n 1`

##rm -r ./coverage
##mkdir coverage

###Writes coverage file for each split .bam file.

cd ${OUT_DIR}/${SMPLE}/bowtie2/alignments/
#rm *.coverage
#echo "${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}*_sorted.REF*" >> $ERRORS/${SMPLE}_error.log

if [[ -f "${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr_sorted.bam" ]]; then
    for f in ${SMPLE}*_sorted.REF*
        do
#        samtools depth $f > ${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/${f}.coverage
        module load bedtools2
        bedtools genomecov -bg -ibam $f > ${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/${f}.coverage
        sed -i -r 's/(\s+)?\S+//2' ${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/${f}.coverage
    done
    touch ${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/coverage.txt
else
    echo "Error in calculating coverage. No reference sorted bam files found." >> $ERRORS/${SMPLE}_error.log
fi

###Removes any unmapped files that may have been written, so that downstream code does not break.
if [[ -f "${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr_sorted.REF_unmapped.bam" ]]; then
    rm ${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr_sorted.REF_unmapped.bam; else
    echo "DID NOT FIND UNMAPPED" >> $ERRORS/${SMPLE}_error.log
fi

#________________________________________________________________________________________|
###                                                                                      |
### This changes the names of the files to remove uninformative text from the filenames. |
###                                 EDIT FOR EACH RUN                                    |
###                       '#':Remove prefix       '%':Remove suffix                      |
#----------------------------------------------------------------------------------------|
cd ${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/

if [[ -f "${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/coverage.txt" ]]; then
    for f in RAP*_NC*
        do
        echo $f
        samp_pref="${f#*146201_*}"
        samp_suff="${samp_pref%_i5*}"
        echo $samp_suff
        acc_pre="${f#*sorted.*}"
        acc_suff="${acc_pre%.bam*}"
        echo $acc_suff
        PLOT_NAME="${samp_suff}_${acc_suff}"
        echo $PLOT_NAME
        mv $f "${PLOT_NAME}.bam.coverage"
    done
##    rm ${SMPLE}_altspec_88.5_sorted.REF_unmapped.bam.coverage
else
    echo "Error in renaming files. No coverage files found." >> $ERRORS/${SMPLE}_error.log
fi

cd ${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/

if [[ -f "${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/${SMPLE}_vircom_nr_sorted.REF_unmapped.bam.coverage" ]]; then
    rm *REF_unmapped.bam.coverage
else
    echo "No unmapped coverage" >> $ERRORS/${SMPLE}_error.log
fi

#This calls the Rscript to plot coverage graphs

if [[ -f "${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/${PLOT_NAME}.bam.coverage" ]]; then
    module load R
    Rscript $SCRIPT_DIR/cov_for_graph.R $OUT_DIR $SMPLE --save
    touch svg.txt
else
    echo "Error in generating coverage graphs. No coverage information used to build graph" >> $ERRORS/${SMPLE}_error.log
fi


if [[ -f "${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/svg.txt" ]]; then
    cp ${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/*.svg ${OUT_DIR}/all_cov
else
    echo "Error in copying coverage graphs to common folder. No coverage graph." >> $ERRORS/${SMPLE}_error.log
fi 
