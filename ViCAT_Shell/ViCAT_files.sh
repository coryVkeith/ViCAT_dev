#!/bin/sh
export CWD=$PWD
#___________________________________________________________________________|
###                                                                         |
###      This script writes 4 files total; 3 files are in the format used   |
###      by phyloseq of which only the metadata needs editing. The other    |
###      file tabulates the contigs of all of the samples in the batch.     |
###      This should only be run, after all other samples in the batch      |
###      have been run through the ViCAT piepline.                          |
#___________________________________________________________________________|

# Reading parameters

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "Writes files for community analysis."
   echo
   echo "Usage: ViCAT_files.sh -s Sample1 -o usr/ViCAT/out"
   echo "options:"
   echo "-h     Print this Help." 
   echo "-o     Required: output directory to write samples. If not a current directory, one will be written."
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

## colon denotes required
while getopts "hs:o:" OPTION
do
   case $OPTION in
       h)
         # if -h, print help function and exit
         Help
         exit 0
         ;;
       o)
         OUT_DIR=$OPTARG
         ;;
       ?)
         echo "ERROR: invalid argument."
         Help
         exit -1
         ;;
     esac
done

#Initialize the directory structure

SAMPLE_DIR="$OUT_DIR/$SMPLE"
ERRORS="$OUT_DIR/errors"
SCRIPT_DIR="$PWD/scripts"
RA="${OUT_DIR}/relative_abundance/"

mkdir -p "$SAMPLE_DIR"
mkdir -p "$ERRORS"
mkdir -p "$RA"

cd $RA

# Deduplicates the files and writes taxa matrix
sort -u seqdescription.tsv -o taxmat.csv

# Adds headers to the taxa matrix
sed -i 1i"RefID,Species" taxmat.csv 

if [ -f *_blasthits.fasta ]; then
    rm *blasthits.fasta 
fi

cd $CWD/$OUT_DIR/virus_contigs/
#________________________________________________________________________________________|
###                                                                                      |
### This changes the names of the files to remove uninformative text from the filenames. |
###                                 EDIT FOR EACH RUN                                    |
###                 '#':Remove before match       '%':Remove after match                 |
#----------------------------------------------------------------------------------------|
for f in *contigs_blasthits.fasta
    do
    #samp_pref="${f#*146201_*}"
    samp_name="${f%_contigs*}"
    grep ">" $f >> ${samp_name}_contigs.txt
        while read l;
            do
            name="${l%_length*}"
            inter="${l#>NODE_*_}"
            length="${inter%_cov*}"
            coverage="${inter#length_*_}"
            node_num="${name#>NODE_}"
            length_num="${length#length_}"
            cov_num="${coverage#cov_}"
            echo "$samp_name,$node_num,$length_num,$cov_num" >> virus_contig_table.csv
         done < ${samp_name}_contigs.txt
done

sed -i 1i"SampleID,Node_Number,Contig_Length,k-mer_coverage" virus_contig_table.csv
rm *contigs.txt

### This code calls the R code to build the metadata table.
cd $CWD/${OUT_DIR}/relative_abundance

### Allows conda environments to be called from within script.
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate R
#Calls Rscript to write files for phyloseq
Rscript $SCRIPT_DIR/metadata.R

