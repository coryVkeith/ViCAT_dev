#/bin/sh

#___________________________________________________________________________|
###                                                                         |
###      This script assemble the reads and filter the contigs according    |
###                   their length.                                         |
#___________________________________________________________________________|

# Reading parameters

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "De novo assembles reads using SPAdes, and then filters contigs on size and k-mer coverage."
   echo
   echo "Usage: ViCAT_assembly.sh -s Sample1 -o usr/ViCAT/out -f usr/data/Sample1_R1.fq.gz -r usr/data/Sample1_R2.fq.gz"
   echo "options:"
   echo "-h     Print this Help."
   echo "-s     Required: Name of sample as it relates to reads. A directory will be created for the sample name."
   echo "-f     Required: path to forward reads."
   echo "-r     Required: path to reverse reads."
   echo "-o     Required: output directory to write samples. If not a current directory, one will be written."
   echo "-c     Optional: k-mer coverage minimum threshold to filter. Default: 30x"
   echo "-M     Optional: maximum length of contigs to filter. Default: 20,000nt"
   echo "-m     Optional: minimum length of contigs to filter. Default: 150nt"
   echo "-t     Optional: threads to use for SPAdes de novo assembly. Default: 16"
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

## Default parameters for optional arguments
COV=30
MAX_LEN=20000
MIN_LEN=150
THREAD=16
## colon denotes required
while getopts "hs:f:r:o:c:M:m:t:" OPTION
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
       c)
	 COV=$OPTARG
	 ;;
       M)
         MAX_LEN=$OPTARG
         ;;
       m)
         MIN_LEN=$OPTARG
         ;;
       t)
         THREAD=$OPTARG
         ;;
       ?)
         echo "ERROR: invalid argument."
         Help
         exit -1
         ;;
     esac
done


# initialize the directory structure

SAMPLE_DIR="$OUT_DIR/$SMPLE"
SPADESOUT="$SAMPLE_DIR/spades"
CONTIGFILT="$SAMPLE_DIR/spades_filtered"
ALLSPADES="$OUT_DIR/all_spades"
ERRORS="$OUT_DIR/errors"
LOG="$OUT_DIR/log"

mkdir -p "$SAMPLE_DIR"
mkdir -p "$SPADESOUT"
mkdir -p "$CONTIGFILT"
mkdir -p "$ALLSPADES"
mkdir -p "$ERRORS"
mkdir -p "$LOG"

OUT_FILE="/xdisk/uschuch/corykeith/BLAST/NCBI_virus/test_1line.fsa"
echo -e "ViCAT_assembly: The parameters for the sample: rnaviralspades.py -t $THREAD -o "${OUT_DIR}/${SMPLE}/spades" -1 $F1 -2 $R1.\nCoverage depth filtering was set to: $COV.\nThe Minimum contig length was set to: $MIN_LEN.\nThe Maximum contig length was set to: $MAX_LEN" > ${LOG}/${SMPLE}_log.txt


# Run spades on the sample

rnaviralspades.py -t $THREAD -o "${OUT_DIR}/${SMPLE}/spades" -1 $F1 -2 $R1
#-----------------------------------------------------------------------------------------------------------|
###                                                                                                         |
### This code filters the spades contigs by 30x coverage and 250bp - 5000bp lengths for faster blast times. |
###                        To change replade the numbers at the awk step in the pipe.                       |
###                                                                                                         |
#-----------------------------------------------------------------------------------------------------------|

### Set variables to be used in if statements and column numbers for awk variables.
contig="$OUT_DIR/$SMPLE/spades/contigs.fasta"
scaffold="$OUT_DIR/$SMPLE/spades/scaffolds.fasta"
COL1=4
COL2=5

#seqkit fx2tab $scaffold | csvtk mutate -H -t -f 1 -p "cov_(.+)" | csvtk mutate -H -t -f 1 -p "length_([0-9]+)" | awk -v c1="$COL1" -v c2="$COL2" -v cov="$COV" -v min="$MIN_LEN" -v max="$MAX_LEN" -F "\t" '$c1>=cov && $c2>=min && $c2<max' | seqkit tab2fx >> "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta" 

### Conditional statment to use scaffolds if SPAdes assembles scaffolds, and contigs if not. Prints an error message if neither scaffolds or contigs were written, as to not break the pipe.
if [[ -f "$scaffold" ]]; then
     seqkit fx2tab $scaffold | csvtk mutate -H -t -f 1 -p "cov_(.+)" | csvtk mutate -H -t -f 1 -p "length_([0-9]+)" | awk -v c1="$COL1" -v c2="$COL2" -v cov="$COV" -v min="$MIN_LEN" -v max="$MAX_LEN" -F "\t" '$c1>=cov && $c2>=min && $c2<max' | seqkit tab2fx >> "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta"; elif [[ -f "$contig" ]]; then
     seqkit fx2tab $contig | csvtk mutate -H -t -f 1 -p "cov_(.+)" | csvtk mutate -H -t -f 1 -p "length_([0-9]+)" | awk -v c1="$COL1" -v c2="$COL2" -v cov="$COV" -v min="$MIN_LEN" -v max="$MAX_LEN" -F "\t" '$c1>=cov && $c2>=min && $c2<max' | seqkit tab2fx >> "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_contigs.fasta"; else
            echo "Error in ViCAT_assembly.sh. No contigs were assembled from the sample. Check the sample file in the raw data." >> $ERRORS/${SMPLE}_error.log
fi

# Copies the filtered spades files to a common folder
if [[ -f "$scaffold" ]]; then
    cp "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta" "$ALLSPADES"; elif [[ -f "$contig" ]]; then
        cp "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta" "$ALLSPADES"; else
            echo "No contig files found, could not copy." >> $ERRORS/${SMPLE}_error.log
fi

# Job done
echo "------------------------------------------"
echo "---------- Assembly successful -----------"
echo "------------------------------------------"

