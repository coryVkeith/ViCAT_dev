#!/bin/sh
export CWD=$PWD
#___________________________________________________________________________|
###                                                                         |
###      This script builds a bowtie2 index from the blast hit reference    |
###      fasta, maps the reads to the reference index, converts the sam     |
###      to a sorted bam, splits the bam on the number of references in     |
###      the index, calls the statistics of reads mapped along the genome   |
###      of the reference, and plots a coverage graph for each reference    |
###      in the sample.                                                     |
#___________________________________________________________________________|

# Reading parameters

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "Read mapping to blast hit reference index and plots coverage depth."
   echo
   echo "Usage: ViCAT_readmap.sh -s Sample1 -o usr/ViCAT/out -f usr/data/Sample1_R1.fq.gz -r usr/data/Sample1_R2.fq.gz"
   echo "options:"
   echo "-h     Print this Help."
   echo "-s     Required: Name of sample as it relates to reads. A directory will be created for the sample name if one does not already exist."
   echo "-f     Required: path to forward reads."
   echo "-r     Required: path to reverse reads."
   echo "-o     Required: output directory to write samples. If not a current directory, one will be written."
   echo "-x     Optional: path and name to bowtie2 index"
   echo
}


## Check number of parameter passed
## $# is literal number of arguments given
## -lt is less than
if [[ $# -lt 1 ]]
then
  echo "Missing function parameter"
  Help
  exit 0
fi
## Optional Argument variable
IDX="$OUT_DIR/$SMPLE/blast_ref_fasta/bowtie2index/${SMPLE}_nr_idx"
## colon denotes required
while getopts "hs:f:r:o:x:" OPTION
do
   case $OPTION in
       h)
         # if -h, print help function and exit
         Help
         exit 0
         ;;
       s)
         SMPLE=$OPTARG
         ;;
       f)
         F1=$OPTARG
         ;;
       r)
         R1=$OPTARG
         ;;
       o)
         OUT_DIR=$OPTARG
         ;;
       x)
         IDX=$OPTARG
         ;;
       ?)
         echo "ERROR: invalid argument."
         Help
         exit -1
         ;;
     esac
done

#Initialize the directory structure

SCRIPT_DIR="$PWD/scripts"
SAMPLE_DIR="$OUT_DIR/$SMPLE"
ERRORS="$OUT_DIR/errors"
LOG="$OUT_DIR/log"
BLAST_OUT="$OUT_DIR/$SMPLE/blast"
BLAST_REF="$OUT_DIR/$SMPLE/blast_ref_fasta/bowtie2index"
STATS="${OUT_DIR}/${SMPLE}/bowtie2/alignments/stats/"
COVER="${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/"
UNUSED="${OUT_DIR}/${SMPLE}/bowtie2/unused_reads/"
CONSENSUS="${OUT_DIR}/${SMPLE}/bowtie2/consensus/"
RA="${OUT_DIR}/relative_abundance/"
ALLCON="$OUT_DIR/all_consensus/"
ALLCOV="$OUT_DIR/all_coverage/"

mkdir -p "$SAMPLE_DIR"
mkdir -p "$ERRORS"
mkdir -p "$LOG"
mkdir -p "$BLAST_OUT"
mkdir -p "$BLAST_REF"
mkdir -p "$STATS"
mkdir -p "$UNUSED"
mkdir -p "$CONSENSUS"
mkdir -p "$RA"
mkdir -p "$ALLCON"
mkdir -p "$ALLCOV"
mkdir -p "$COVER"

#
#####DEVVVV
IDX="$OUT_DIR/$SMPLE/blast_ref_fasta/bowtie2index/${SMPLE}_nr_idx"
module load bowtie2

##bowtie2 build index of filtered references above
if [[ -f "$OUT_DIR/$SMPLE/blast_ref_fasta/${SMPLE}_nr_ref.fasta" ]]; then
    bowtie2-build -f $OUT_DIR/$SMPLE/blast_ref_fasta/${SMPLE}_nr_ref.fasta $OUT_DIR/$SMPLE/blast_ref_fasta/bowtie2index/${SMPLE}_nr_idx
    else
    echo "Error in run_bowtie2buildnmap.sh, bowtie2. No CD-HIT clustered reference fasta file found to build index. Stopping at build." >> $ERRORS/${SMPLE}_error.log
fi
##bowtie2 reference mapping to indexes
if [[ -f "$OUT_DIR/$SMPLE/blast_ref_fasta/${SMPLE}_nr_ref.fasta" ]]; then    
    bowtie2 --un-conc ${OUT_DIR}/${SMPLE}/bowtie2/unused_reads/"${SMPLE}_vircom_rm.fastq" -x $IDX -q -1 $F1 -q -2 $R1 -S ${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr.SAM
    else
    echo "Error in run_bowtie2buildnmap.sh, bowtie2. No reference fasta file found to build index. Stopping at mapping" >> $ERRORS/${SMPLE}_error.log
fi
##Samtools convert to sorted bam, extract consensus, and remove sam
if [[ -f "${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr.SAM" ]]; then
    samtools view -bo ${OUT_DIR}/${SMPLE}/bowtie2/alignments/"${SMPLE}_vircom_nr.bam" ${OUT_DIR}/${SMPLE}/bowtie2/alignments/"${SMPLE}_vircom_nr.SAM"
    samtools sort ${OUT_DIR}/${SMPLE}/bowtie2/alignments/"${SMPLE}_vircom_nr.bam"  -o ${OUT_DIR}/${SMPLE}/bowtie2/alignments/"${SMPLE}_vircom_nr_sorted.bam"
    samtools mpileup --max-depth 0 -uf ${OUT_DIR}/${SMPLE}/blast_ref_fasta/${SMPLE}_nr_ref.fasta ${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${OUT_DIR}/${SMPLE}/bowtie2/consensus/${SMPLE}_vircom_nr_consensus.fasta
#
##_____________________________________________________________________________________________________________________________________|    
####    This code will copy the consensus sequences to a common folder and remove the .sam and unsorted .bam files to help save space. |
####                                     Comment this code if you want to keep intermediate files.                                     |
####                                                                                                                                   |
##_____________________________________________________________________________________________________________________________________|    
    cp ${OUT_DIR}/${SMPLE}/bowtie2/consensus/${SMPLE}_vircom_nr_consensus.fasta "$ALLCON"
    rm ${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr.SAM
    rm ${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr.bam
    else
    echo "Error in run_bowtie2buildnmap.sh, samtools. No sam alignment file found. Stopping at bam file conversion." >> $ERRORS/${SMPLE}_error.log
fi
#
### Splits the sorted bam files
#### This is not calling the correct directory so I changed the file path to have it call the CWD. I need to test this.
##echo ${OUT_DIR}/${SMPLE}/blast_ref_fasta/"${SMPLE}_blasthits.fasta"
##echo "bamtools split -in $OUT_DIR/$SMPLE/bowtie2/alignments/${SMPLE}_vircom_nr_sorted.bam -reference"
bamtools split -in $CWD/$OUT_DIR/$SMPLE/bowtie2/alignments/${SMPLE}_vircom_nr_sorted.bam -reference
#
###_________________________________________________________________________________|
#####                                                                               |
#####               This code provides the stats for the read mapping.              |
#####    Change the file name to match the regular expression of the BLASTdb used.  |
#####                    Also change the Suff and OUT variable.                     |
####________________________________________________________________________________|
##
for f in ${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr_sorted.REF_NC*.bam
    do
    Suff=${f%%.bam}
    OUT=${Suff##*/}
#    echo $OUT
    samtools stats $f >> ${OUT_DIR}/${SMPLE}/bowtie2/alignments/stats/${OUT}_stats.txt
done
##
#######This code takes just the number of reads that map from the stats file generated above, and the name of the sample.
for s in ${OUT_DIR}/${SMPLE}/bowtie2/alignments/stats/${SMPLE}_*stats.txt
    do
    Suff=${s##*/}
    HIT="${Suff%%_stats.txt}"
    awk 'NR==14' $s | sed 's/ \+ /\t/g' | awk '{ print $4 }' | xargs printf "$HIT \t %s\n" >> ${OUT_DIR}/${SMPLE}/bowtie2/alignments/stats/"${SMPLE}_readscount.tsv"
done
##
####### This code combines the individual reads counts into the combined read counts
for f in ${OUT_DIR}/${SMPLE}/bowtie2/alignments/stats/*_readscount.tsv
    do
    while read line
        do
        printf "$line\n" >> ${OUT_DIR}/relative_abundance/combined_readscount.tsv
    done <$f
done
##
####Writes coverage file for each split .bam file.
#
#
##rm *.coverage
##echo "${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}*_sorted.REF*" >> $ERRORS/${SMPLE}_error.log
cd ${OUT_DIR}/${SMPLE}/bowtie2/alignments/
if [[ -f "${SMPLE}_vircom_nr_sorted.bam" ]]; then
##    pwd
    for f in ${SMPLE}*_sorted.REF*
        do
        bedtools genomecov -bg -ibam $f > ./coverage/${f}.coverage 
        sed -i -r 's/(\s+)?\S+//2' ./coverage/${f}.coverage
    done
    touch ./coverage/coverage.txt
else
    echo "Error in calculating coverage. No reference sorted bam files found." >> $CWD/$ERRORS/${SMPLE}_error.log
fi
cd $CWD
##pwd
#####Removes any unmapped files that may have been written, so that downstream code does not break.
if [[ -f "${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr_sorted.REF_unmapped.bam" ]]; then
    rm ${OUT_DIR}/${SMPLE}/bowtie2/alignments/${SMPLE}_vircom_nr_sorted.REF_unmapped.bam; else
    echo "DID NOT FIND UNMAPPED" >> $ERRORS/${SMPLE}_error.log
fi

### Allows conda environments to be called from within a script.
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda deactivate
#This calls the R environemnt to plot coverage graphs
conda activate R

#
##________________________________________________________________________________________|
####                                                                                      |
#### This changes the names of the files to remove uninformative text from the filenames. |
####                                 EDIT FOR EACH RUN                                    |
####                       '#':Remove prefix       '%':Remove suffix                      |
##----------------------------------------------------------------------------------------|
#
cd $CWD/${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/
if [[ -f "coverage.txt" ]]; then
    for f in *vircom*
        do
        echo $f
       # samp_pref="${f#*146201_*}"
        samp_suff="${f%_vircom_nr_*}"
        #echo $samp_suff
        acc_pre="${f#*sorted.*}"
        acc_suff="${acc_pre%.bam*}"
        #echo $acc_suff
        PLOT_NAME="${samp_suff}_${acc_suff}"
        #echo $PLOT_NAME
        mv $f "${PLOT_NAME}.bam.coverage"
    done
else
    echo "Error in renaming files. No coverage files found." >> $CWD/$ERRORS/${SMPLE}_error.log
fi

if [[ -f "${SMPLE}_REF_unmapped.bam.coverage" ]]; then
    rm *REF_unmapped.bam.coverage
else
    echo "No unmapped coverage" >> $CWD/$ERRORS/${SMPLE}_error.log
fi

COV_R_DIR=$CWD/$OUT_DIR

if [[ -f "coverage.txt" ]]; then
    Rscript $SCRIPT_DIR/cov_for_graph.R $COV_R_DIR $SMPLE --save &> /xdisk/jkbrown/corykeith/RNAviral_test/ViCAT_shell/R_error.txt
    touch $CWD/${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/svg.txt
else
    echo "Error in generating coverage graphs. No coverage information used to build graph" >> $CWD/$ERRORS/${SMPLE}_error.log
fi

cd $CWD

if [[ -f "${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/svg.txt" ]]; then
    cp ${OUT_DIR}/${SMPLE}/bowtie2/alignments/coverage/*.svg ${OUT_DIR}/all_coverage
else
    echo "Error in copying coverage graphs to common folder. No coverage graph." >> $ERRORS/${SMPLE}_error.log
fi
