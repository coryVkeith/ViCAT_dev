#!/bin/sh

# Your job will use 1 node, 28 cores, and 168gb of memory total.
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 1:00:00

#Load ViCAT environment
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate ViCAT

cd $OUT_DIR/virus_contigs/
#________________________________________________________________________________________|
###                                                                                      |
### This changes the names of the files to remove uninformative text from the filenames. |
###                                 EDIT FOR EACH RUN                                    |
###                 '#':Remove before match       '%':Remove after match                 |
#----------------------------------------------------------------------------------------|
for f in *contigs_blasthits.fasta
    do
    samp_pref="${f#*146201_*}"
    samp_name="${samp_pref%_i5*}"
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
cd ${OUT_DIR}/relative_abundance

module load R
Rscript $SCRIPT_DIR/metadata.R
