export CWD=$PWD

###### #Parameters 

# Location of raw data 
export RAW="/CWD/test_data/sim_resist/mock_reads"
export OUT_DIR="/CWD/out"
export PROFILE="/CWD/profile.txt"

# Location of tools and indices if installed outside of anaconda environemnt, add the path before the commands. If installed in an Anaconda environment, options for running tools can change. Otherwise these are the tools arguments.
export BOWTIE="bowtie2"
export BOWTIEBUILD="bowtie2-build"
export SPADES="rnaviralspades.py"
export SEQKIT="seqkit"
export CSVTK="csvtk"
export BAMTOOLS="bamtools"
export BEDTOOLS="bedtools"

# Location of blast database and 1 line fasta
export BLASTDB="CWD/blastdb/[database]"
export BLAST1line="CWD/blastdb/[database]_1line.fa"

#Place to store scripts
export SCRIPT_DIR="$PWD/scripts"
export GATHER="$PWD/scripts"

# User information
export QUEUE="standard"
export GROUP="[group_name]"
export MAIL_USER="[user_email@email.com]"
export MAIL_TYPE="END"


# --------------------------------------------------
# removes directories with name and makes new directory based on variable.
function init_dir {
    for dir in $*; do
        if [ -d "$dir" ]; then
            rm -rf $dir/*
        else
            mkdir -p "$dir"
        fi
    done
}
