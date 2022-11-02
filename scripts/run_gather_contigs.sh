#!/bin/sh

#Your job will use 1 node, 28 cores, and 64gb of memory total.

#SBATCH -N 1
#SBATCH -n 28
#SBATCH --mem=64gb
#SBATCH -t 1:00:00

# Activates anaconda environment for this code to run.
source activate bio

#Calls the correct sample ID from the profile.txt
export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $PROFILE | tail -n 1`

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
    else
    echo " Error in run_gather_contigs.sh. No contigs or scaffolds found." >> $ERRORS/${SMPLE}_error.log
fi

### The following code builds a table of all of the virus contig information in the virus contigs directory


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
