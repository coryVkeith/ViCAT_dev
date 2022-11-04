# ViCAT_dev
<a href="https://github.com/cvk1988/CLCuD_pop_pipe/graphs/contributors">
<img src="https://contrib.rocks/image?repo=cvk1988/CLCuD_pop_pipe" />
</a>

<img src="https://komarev.com/ghpvc/?username=cvk1988"/> ![Hits](https://hitcounter.pythonanywhere.com/count/tag.svg?url=https://github.com/cvk1988/CLCuD_pop_pipe/) 

--------------------------------------------------------------------------------------------------------------------------------------------------------

# **ViCAT**

## Virus Community Assembly Tool
  A tool that assembles and characterizes virus communities from target enrichment high-throughput sequncing (TE-HTS)  data within and between samples. To be used on an HPC with slurm scheduler. 

--------------------------------------------------------------------------------------------------------------------------------------------------------

## Work Flow

![plot](https://github.com/cvk1988/ViCAT_dev/blob/main/ViCAT_workflow1.png)

## Directory Structure
![plot](https://github.com/cvk1988/ViCAT_dev/blob/main/ViCAT_directory.png)

--------------------------------------------------------------------------------------------------------------------------------------------------------

### Dependencies
- [SPAdes](https://github.com/ablab/spades): at least version 3.14.1
- [SeqKit](https://github.com/shenwei356/seqkit)
- [csvtk](https://github.com/shenwei356/csvtk)
- [Local Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [CD-HIT](http://bioinformatics.org/cd-hit/)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://www.htslib.org/)
- [bamtools](https://github.com/pezmaster31/bamtools)
- R packages: ggplot2; gridExtra; plyr
- Python (3.6 or higher) packages: Bipython; sys

--------------------------------------------------------------------------------------------------------------------------------------------------------

### Installation

If tools are not installed on your HPC download the tools and install. Change the paths to the tools in `config.sh` to include the full path of the newly installed tools.

`SPADES="path/to/spades"`


Alternatively, create an Anaconda environment and download all of the tools with conda.

`conda create -n ViCAT python=3.7`

`conda activate ViCAT`

`pip install [packages]`

`conda install -c bioconda [tool]`

--------------------------------------------------------------------------------------------------------------------------------------------------------

### Quick Usage
*OPTIONAL lines of code are identified by comment boxes in the respective scripts.*

1. Edit the `config.sh` file in the home directory to include paths to:
    - RAW: the path to the raw data for the run.
    - OUT_DIR: the path to the desired location for the outputs of the tool. This will be re-written each time the tool is used, so be careful to backup the results of previous runs.
    - BLASTDB: the location of the blast database to use.
    - BLAST1line: 1 line fasta to use for building reference index in bowtie2 mapping.
    - PROFILE: the path to the file that inlcudes the sample names as they relate to the filenames. An example profile in located in the home directory.
    - [TOOLS]: any paths to tools that were not installed via Anaconda.
    - [USER_INFO]: information for the HPC scheduler.
2. Edit the `run_spades.sh` lines 13 and 14 that set the file suffixes. The file can be found in `CWD/scripts`.

    `F1="${SMPLE}_R1.fastq"`
    
    `R1="${SMPLE}_R2.fastq"`
    
3. OPTIONAL but ***suggested***: edit the `run_covgraph_for.sh` lines 44-61 to remove portions of the sample name which may be uninformative and cause difficulties in reading coverage graph titles in the output files. The file can be found in `CWD/scripts`.
4. OPTIONAL: edit the `run_blastn.sh` file on line 32 to change blast options.
5. OPTIONAL: edit the `run_collect_contigs_accessions.sh` file on lines 25 and 42 to change filtering criteria for output files from Step 3. The file can be found in `CWD/scripts`. The file can be found in `CWD/scripts`.
6. OPTIONAL: edit the `run_cdhit.sh` file on line 20 to change percent similarity of blast hits for clustering. The file can be found in `CWD/scripts`.


To run the tool, go to the tools directory and execute the `run.sh` file.

`./run.sh`

--------------------------------------------------------------------------------------------------------------------------------------------------------

### Explanation of files for usage

- `run.sh` contains the instructions for the workflow and builds the directory structure shown above.
- `config.sh` includes the variables for the paths to the raw data, desired output directory location, and `profile.txt`, which all must be changed for each run. Also includes variables for paths to tools and information for job scheduling using SLURM. Paths to tools only need to be set after initial instalations.
- `profile.txt` contains the information for the samples as they are labeled by the reads. An example for the test data is included in `CWD`.
- `scripts/[scripts.sh]` Directory with all of the instructions for running each individual step shown in the workflow. See Usage for options on what to change.
    - `run_spades.sh`: Instructions for SPAdes *de novo* assembly of reads. Filters assembled contigs by *k*-mer coverage and length and copies the filtered contigs to `OUT_DIR/all_spades`. 
    - `run_blastn.sh`: Instructions for BLAST used for taxonomic assignment. Writes BLAST files and copies to `OUT_DIR/all_blast`
    - `run_collect_contigs_accessions.sh`: Writes contig names that have BLAST hits to a file. Writes accession numbers of BLAST hits to file and then collects RefSeq fasta file for use in building bowtie2 index for read mapping. List of contigs and accessions found in `OUT_DIR/Sample(n)/blast`. RefSeq fasta file found in `OUT_DIR/Sample(n)/blast_ref_fasta`.
    - `run_gather_contigs`: Gathers contigs that have BLAST hits to deposit in `OUT_DIR/virus_contigs` folder.
        - `gather_contigs.py`: python code to execute above task.
    - `run_cdhit`: Removes redundancy of fasta file of virus RefSeq accessions for use in building bowtie2 index. Clustering output and non-redundant fasta files are found in `OUT_DIR/Sample(n)/blast_ref_fasta`.
    - `run_bowtie2buildnmap.sh`:This file makes a bowtie2 index and places it in `OUT_DIR/Sample(n)/blast_ref_fasta/bowtie2index`. The QC reads from the first step of the workflow are then used to map against this reference alignment files can be found in `OUT_DIR/Sample(n)/bowtie2/alignment`. Each sample can have a different index, as the BLAST results indicate the taxonomic assignments used for building the bowtie2 index. This code also converts the `[alignment].sam` to an `[alignment].sorted.bam`, removes the `[alignment].sam` and intermediate `[alignment].bam` files, and extracts the consensus sequences in the sample which is deposited in `OUT_DIR/Sample(n)/bowtie2/consensus`. Copies consensus fastq files to `OUT_DIR/all_consensus`.
    - `run_split`: Splits the `[alignment].sorted.bam` into seperate files based on the RefSeq accession number. Calls the statistics on each split `[alignment].sorted.bam` and writes a hit table of total reads mapped to the respective RefSeq. Statistics and hit tables are found in `OUT_DIR/Sample(n)/bowtie2/alignment/stats`. Writes each viral taxa with reads mapped from each sample to a combined hit table of all samples in the run which is located in `OUT_DIR/relative_abundance`, labeled `combined_readscount.tsv`.
    - `run_covgraph_for.sh`: Numbers reads mapped to viral loci for writing to coverage graph and writes files to `OUT_DIR/Sample(n)/bowtie2/alignment/coverage`. Copies coverage graphs to `OUT_DIR/all_cov`.
        - `cov_for_graph.R`: R code to plot coverage of from files written in above step.   
    - `run_taxmat.sh`: Writes the taxanomic matrix of the submitted samples to `OUT_DIR/relative_abundance`. The file `taxmat.tsv` is a tab separated matrix containing the GenBank RefSeq accession number and sequence description, and is in the format used by Phyloseq in downstream community analysis.
    - `run_metadata.sh`: This converts the `combined_readscount.tsv` in `OUT_DIR/relative_abundance` to the format used by Phyloseq in downstream community analysis and writes to the new file `ViCAT_otutable.csv`. Using the `ViCAT_otutable.csv`, a metadata template file `metadata.csv` is written that retains the sample order and sample ID of the original samples as well as columns of potentially relevant useful metadata for downstream analysis. These values must be changed according to the sample information. These files are found in `OUT_DIR/relative_abundance`. This code also writes a descriptive table of the contigs which have blast hits to viral taxa in the `OUT_DIR/virus_contigs` folder.

### Explanation of directory structure




