#!/bin/sh
export CWD=$PWD
#___________________________________________________________________________|
###                                                                         |
###      This script runs blastn on the assembled contigs for a sample,     |
###      gathers the information on which contigs have blast hits, writes   |
###      the unique GenBank accession numbers to a list, gathers the        |
###      fasta file of the blast hit subject to a new file, and removes     |
###      any redundancy in the reference fasta by percent nucleotide        |
###      identity.                                                          |
#___________________________________________________________________________|

# Reading parameters

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "Blastn of filtered SPAdes de novo assemblies, and builds reference fasta for read mapping."
   echo
   echo "Usage: ViCAT_assembly.sh -s Sample1 -d usr/ViCAT/blastdb/db.fna -f usr/ViCAT/blastdb/1linedb.fa -o usr/ViCAT/out"
   echo "options:"
   echo "-h     Print this Help."
   echo "-s     Required: Name of sample as it relates to reads. A directory will be created for the sample name, unless one has already been created."
   echo "-d     Required: Path to blast database."
   echo "-f     Required: Path to 1 line fasta file of blast reference database."
   echo "-o     Required: Path to output directory to write samples. If not a current directory, one will be written."
   echo "-t     Optional: Number of blast hits to show. Default: 2"
   echo "-M     Optional: Maximum length of blast hit query to subject. Default: 4500"
   echo "-m     Optional: Minimum length of blast hit query to subject. Default: 150"
   echo "-r     Optional: Redundancy threshold. Percent pairwise nt identity to remove fasta sequences from blast reference fasta file."
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
TARGET=2
MIN_MATCH=150
MAX_MATCH=4500
REDUND=0.885

## colon denotes required
while getopts "hs:d:f:o:t:M:m:r:" OPTION
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
       d)
         BLASTDB=$OPTARG
         ;;
       f)
         BLAST1line=$OPTARG
         ;;
       o)
         OUT_DIR=$OPTARG
         ;;
       t)
         TARGET=$OPTARG
         ;;
       M)
         MAX_MATCH=$OPTARG
         ;;
       m)
         MIN_MATCH=$OPTARG
         ;;
       r)
         REDUND=$OPTARG
         ;;
       ?)
         echo "ERROR: invalid argument."
         Help
         exit -1
         ;;
     esac
done



# initialize the directory structure
SCRIPT_DIR="$PWD/scripts"
SAMPLE_DIR="$OUT_DIR/$SMPLE"
SPADESOUT="$SAMPLE_DIR/spades"
CONTIGFILT="$SAMPLE_DIR/spades_filtered"
ALLSPADES="$OUT_DIR/all_spades"
ERRORS="$OUT_DIR/errors"
LOG="$OUT_DIR/log"
BLAST_OUT="$OUT_DIR/$SMPLE/blast"
ALLBLAST="$OUT_DIR/all_blast"
BLAST_REF="$OUT_DIR/$SMPLE/blast_ref_fasta/"
VIRUS="$OUT_DIR/virus_contigs"
RA="${OUT_DIR}/relative_abundance"

mkdir -p "$SAMPLE_DIR"
mkdir -p "$SPADESOUT"
mkdir -p "$CONTIGFILT"
mkdir -p "$ALLSPADES"
mkdir -p "$ERRORS"
mkdir -p "$LOG"
mkdir -p "$BLAST_OUT"
mkdir -p "$ALLBLAST"
mkdir -p "$BLAST_REF"
mkdir -p "$VIRUS"
mkdir -p "$RA"

## This code checks if there are scaffolds or contigs and then sets the variables for use in the blast code below. If no files are present, an error log is printed.

if [[ -f "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta" ]]; then
    export f="$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta"; elif [[ -f "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_contigs.fasta" ]]; then
        export f="$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_contigs.fasta"; else
            echo "Error in run_blastn.sh. There are not filtered contigs or scaffolds." >> $ERRORS/${SMPLE}_error.log
fi
export OUTF="$OUT_DIR/$SMPLE/blast/${SMPLE}_blastn_results.tsv"

#___________________________________________________________________________________________|
###                                                                                         |
### This is the blastn algorithm. To allow more dissimilar sequences: '-task dc-megablast'. | 
###          To change the max number of hits per query: '-max_target_seqs [#]'             |
###                                                                                         |
#___________________________________________________________________________________________|

#cd $OUT_DIR/$SMPLE/spades_filtered/
#f="${SMPLE}_filtered_scaffolds.fasta"

## Blast algorithm
if [[ -f "$f" ]]; then    
    #echo "blastn -db $BLASTDB -query $f -out $OUTF -evalue 0.01 -max_target_seqs $TARGET -outfmt "6 qseqid sseqid ssciname pident length mismatch gapopen qstart qend sstart send evalue bitscore"" 
    blastn -db $BLASTDB -query $f -out $OUTF -evalue 0.01 -max_target_seqs $TARGET -outfmt "6 qseqid sseqid ssciname pident length mismatch gapopen qstart qend sstart send evalue bitscore"; else
        echo "Error in run_blastn.sh. No contigs or scaffolds found to blast." >> $ERRORS/${SMPLE}_error.log
fi
# Copies blast files to a common folder
if [[ -f "$OUTF" ]]; then
    cp "$OUTF" "$ALLBLAST"    
fi

# File Check
if [[ -f "$OUT_DIR/$SMPLE/blast/${SMPLE}_blastn_results.tsv" ]]; then
    f="$OUT_DIR/$SMPLE/blast/${SMPLE}_blastn_results.tsv"; else
        echo "Error in run_collect_contigs_accessions.sh. No blast results to write accessions. Stopping at writing accessions." >> $ERRORS/${SMPLE}_error.log
fi

#__________________________________________________________________________________________________________________________________|
###                                                                                                                                |
### filters hits that only have a match length of greater than 150 and less than 4500 and writes the accession number to a list.   |
###                                  Change the awk code in the first step of the pipe for different filters.                      |
#__________________________________________________________________________________________________________________________________|
COL1=5
COL2=2
echo $f >> $ERRORS/${SMPLE}_error.log
awk -v col1="$COL1" -v col2="$COL2" -v min="$MIN_MATCH" -v max="$MAX_MATCH" '$col1>min && $col1<max {print $col2}' $f | sed '$!N; /^\(.*\)\n\1$/!P; D'  >> "$OUT_DIR/$SMPLE/blast/${SMPLE}_unique_accessions.tsv"
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
COL3=1
awk -v col1="$COL1" -v col2="$COL2" -v min="$MIN_MATCH" -v max="$MAX_MATCH" -v col3="$COL3" '$col1>min && $col1<max {print $col3}' $f | sed '$!N; /^\(.*\)\n\1$/!P; D'  >> "$OUT_DIR/$SMPLE/blast/${SMPLE}_unique_contigs.tsv"
# This code writes the fasta file from the unique accessions file from the blast hits step to a new fasta file for use in building bowtie indexes

if [[ -f "$OUT_DIR/$SMPLE/blast/trimmed_${SMPLE}_unique_accessions.tsv" ]]; then
    while read line; do
        grep ">$line" -A 1 $BLAST1line | cat >> $OUT_DIR/$SMPLE/blast_ref_fasta/${SMPLE}_blasthits.fasta
    done <$OUT_DIR/$SMPLE/blast/trimmed_${SMPLE}_unique_accessions.tsv
    else
    echo "Error in run_collect_contigs_accessions.sh. No accessions found in file. Stopping at collecting accessions." >> $ERRORS/${SMPLE}_error.log
fi

if [[ -f "$OUT_DIR/$SMPLE/blast_ref_fasta/${SMPLE}_blasthits.fasta" ]]; then
    cp "$OUT_DIR/$SMPLE/blast_ref_fasta/${SMPLE}_blasthits.fasta" "$OUT_DIR/relative_abundance"
fi

# These are the variable files for the python tool, gather_contigs.py
if [[ -f "$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta" ]]; then
    export F1="$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_scaffolds.fasta"
    else
    export F1="$OUT_DIR/$SMPLE/spades_filtered/${SMPLE}_filtered_contigs.fasta"
fi
export F2="$OUT_DIR/$SMPLE/blast/${SMPLE}_unique_contigs.tsv"
export F3="$OUT_DIR/virus_contigs/${SMPLE}_contigs_blasthits.fasta"

#python code for the gather_contigs.py tool.

if [[ -f "$F1" ]]; then
    python $SCRIPT_DIR/gather_contigs.py $F1 $F2 $F3
    echo $F1 $F2 $F3
    else
    echo " Error in run_gather_contigs.sh. No contigs or scaffolds found." >> $ERRORS/${SMPLE}_error.log
fi

### Removes redundancy of reference index fasta
if [[ -f "${OUT_DIR}/${SMPLE}/blast_ref_fasta/"${SMPLE}_blasthits.fasta"" ]]; then
    cd-hit-est -i ${OUT_DIR}/${SMPLE}/blast_ref_fasta/"${SMPLE}_blasthits.fasta" -o ${OUT_DIR}/${SMPLE}/blast_ref_fasta/"${SMPLE}_nr_ref.fasta" -c $REDUND
    else
    echo "Error in run_cdhit. No blast file found. Did not cluster sequences." >> $ERRORS/${SMPLE}_error.log
fi

# Change working directory
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
