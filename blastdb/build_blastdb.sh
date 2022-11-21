#!/bin/sh


## sets variables
ALL_FASTA="viraldb.fsa"
ONE_LINE="viraldb_1line.fsa"
DB_NAME="viraldb"

## This code concatenates all of the genbank refseq files into one .fsa file.
cat *.fna > $ALL_FASTA
 
## This code converts multiline fasta files to one line fasta files for using awk.
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $ALL_FASTA > $ONE_LINE

## This code builds the blastdb from the files created above
module load blast
makeblastdb -in $ALL_FASTA -parse_seqids -blastdb_version 5 -title "$DB_NAME" -dbtype nucl
