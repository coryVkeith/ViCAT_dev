# ViCAT_dev
<a href="https://github.com/cvk1988/CLCuD_pop_pipe/graphs/contributors">
<img src="https://contrib.rocks/image?repo=cvk1988/CLCuD_pop_pipe" />
</a>

<img src="https://komarev.com/ghpvc/?username=cvk1988"/> ![Hits](https://hitcounter.pythonanywhere.com/count/tag.svg?url=https://github.com/cvk1988/CLCuD_pop_pipe/) 

--------------------------------------------------------------------------------------------------------------------------------------------------------

# **ViCAT**

## Virus Community Assembly Tool
  A tool that assembles and characterizes virus communities from target enrichment high-throughput sequncing (TE-HTS) data within and between samples. To be used on an HPC with slurm scheduler. 

--------------------------------------------------------------------------------------------------------------------------------------------------------

## Work Flow

![plot](https://github.com/coryVkeith/ViCAT_dev/blob/main/figs/ViCAT_workflow.png)

## Directory Structure
![plot](https://github.com/coryVkeith/ViCAT_dev/blob/main/figs/ViCAT_directory.png)

--------------------------------------------------------------------------------------------------------------------------------------------------------




## Installation

This pipeline has been tested on CentOS linux; should work in all linux distributions.

### Option 1 (conda version)

Conda is the easiest way to install ViCAT. If you do not have conda installed, it can be installed following this link.

```
# Clone repository
git clone git@github.com:coryVkeith/ViCAT_dev.git

cd ViCAT_dev/scripts

# Create Conda environment
conda env create -f ViCAT_v0.1.yml
conda activate ViCAT
```

### Option 2 (Manual installation of each dependancies and tools)

ViCAT relies on several external tools that can be individually installed:
- [SPAdes](https://github.com/ablab/spades): at least version 3.15.1
- [SeqKit](https://github.com/shenwei356/seqkit)
- [csvtk](https://github.com/shenwei356/csvtk)
- [Local Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [CD-HIT](http://bioinformatics.org/cd-hit/)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://www.htslib.org/)
- [bamtools](https://github.com/pezmaster31/bamtools)
- R packages: ggplot2; gridExtra; plyr, deseq2, phyloseq
- Python (3.6 or higher) packages: Bipython; sys

## Download database and dependencies

Before running ViCAT, users must download the RefSeq reference database, this task only needs to be done once. Alternatively, users can generate a custom blast database of reference.

```
cd ViCAT_dev/blastdb

bash build_blastdb.sh

# the database should take few minutes to get created
```

## Run ViCAT

The complete ViCAT pipeline contains 4 successive steps (1 optional step). Each step is contained in a separate function described below. 

This tutorial uses a simulated viral mock community that has differential relative abundance in response to host background have been included in the folder test_data. Simulated host resistance was based on the assumption that viruses for which a host is resistant will have a lower relative abundance of reads mapped.

###  1. Assembly and contig filtering

The first step of the ViCAT pipeline is the assembly of the reads and contig length filtering. ViCAT leverages SPAdes as an assembler. If the users have generated their own assembly, they may skip this step.


**VICAT_assembly.sh usage:**

```
VICAT_assembly.sh -s ${SAMPLE NAME} -f ${FORWARD READ} -r ${REVERSE READ} -o ${OUTPUT DIRECTORY}
```

As an example to run on one sample :

```
cd ViCAT_dev/scripts
TEST_SET="../test_data/sim_resist/mock_reads"

bash ViCAT_assembly.sh -s begomovirus_P10 -f $TEST_SET/begomovirus_P10_R1.fastq -r $TEST_SET/begomovirus_P10_R2.fastq -o ../test_ViCAT
```
 
**Optional arguments:**

    -c: k-mer coverage minimum threshold to filter. Default is 30x.
    -M: Maximum length of contig to filter. Default is 20,000nt.
    -m: Minimum length of contig to filter. Default is 150nt.
    -t: Number of threads to use for SPAdes de novo assembly. Default is 16. 

### 2. Taxonomic assignment

The second step of the ViCAT pipeline assigns taxonomy to the ***de novo*** assembled contigs from the step above. This is done by a blastn to the GenBank Viral RefSeq database that the user previously installed. The reference fasta sequence of the queried blast hits for each contig is written to a file for each contig in the sample, and then the RefSeq fasta file is removed of sequences that share a higher nucleotide percent identity than a user defined threshold.

**ViCAT_taxassign.sh usage:**

```
ViCAT_taxassign.sh -s $[SAMPLE NAME} -d ${BLAST DB PATH} -f ${REFSEQ 1LINE FASTA} -o ${OUTPUT DIRECTORY}
```

As an example to run on one sample :

```
cd ViCAT_dev/scripts
BLASTDB="../blastdb/viraldb.fsa"
FASTA="../blastdb//viraldb_1line.fsa"

./ViCAT_taxassign.sh -s begomovirus_P10 -d $BLASTDB -f $FASTA -o ../test_ViCAT
```

**Optional arguments:**

    -t: Number of blast hits to show. Default is two. More than two hits is discouraged, as some viruses undergo extensive recombination and this may lead to incomplete coverage graphs and read binning in subsequent steps.
    -M: Maximum length of blast hit query to subject to select query fasta files for downstream read mapping. Default is 4,500 as is set for ~1,500 nt longer than begomvirus full-length molecules to account for any concatemeric regions or misassemblies.
    -m: Minimum length of blost hit query to subject. Default is 150.
    -r: Redundancy threshold for building fasta file of RefSeq sequences to use in bowtie2 indexes in subsequent steps. Default is 88.5% PNI.
    
### 3. Read binning

The third step of ViCAT builds a bowtie2 index from the RefSeq fasta files generated in the taxonomic assignment step and then maps the sample reads to this index. It then sorts the generated alignment files, extracts consensus sequences, writes the number of reads mapped to files,  calls statistics, and plots coverage graphs of reads mapped to individual virus taxa within a sample.

**ViCAT_readmap.sh usage:**

```
ViCAT_readmap.sh -s ${SAMPLE} -f ${FORWARD READS} -r ${REVERSE READS} -o ${OUTPUT DIRECTORY}
```

As an example to run on one sample :

```
cd ViCAT_dev/scripts
TEST_SET="../test_data/sim_resist/mock_reads"

ViCAT_readmap.sh -s begomovirus_P10 -f $TEST_SET/begomovirus_P10_R1.fastq -r $TEST_SET/begomovirus_P10_R2.fastq -o ../test_ViCAT
```

**Optional argument:**

    -x: Path and name to bowtie2 index.

### 4. Write files for community analysis
The fourth step of ViCAT writes files for use in the optional community characterization and analysis step of ViCAT. This step should be run after all samples in the set have been submitted to ViCAT. This step formats all of the data to be fed into the R-script, 'Phyloseq_Community_Characterization.R', found in ViCAT_dev/Community_Assembly. This script utilizes the [phyloseq](https://joey711.github.io/phyloseq/) package in R.

**ViCAT_files.sh usage:**

```
ViCAT_files.sh
```






# THIS IS THE OLD USAGE!!!

### Dependencies
- [SPAdes](https://github.com/ablab/spades): at least version 3.15.1
- [SeqKit](https://github.com/shenwei356/seqkit)
- [csvtk](https://github.com/shenwei356/csvtk)
- [Local Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [CD-HIT](http://bioinformatics.org/cd-hit/)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://www.htslib.org/)
- [bamtools](https://github.com/pezmaster31/bamtools)
- R packages: ggplot2; gridExtra; plyr, deseq2, phyloseq
- Python (3.6 or higher) packages: Bipython; sys

--------------------------------------------------------------------------------------------------------------------------------------------------------
## Database

Before ViCAT can be run, a database needs to be created for blast. The instructions and a code for building custom databases are included in `CWD/blastdb`. This only needs to be done any time a different database is required to run.

## Installation

### Anaconda environment yaml
***Recommended***

`ViCAT.yml` is included in the home directory. If this is the installation, do not change the lines in the `config.sh` file that include the paths to the tools.

To install an environment containing all of the tools:
    `conda env create -f ViCAT.yml`

### Anaconda
Create an Anaconda environment and download all of the tools with conda or pip. If this is the installation, do not change the lines in the `config.sh` file that include the paths to the tools.

`conda create -n ViCAT python=3.7`

`conda activate ViCAT`

`pip install [packages]`
OR
`conda install -c bioconda [tool]`

### Manual installion of tools

If tools are not installed on your HPC download the tools and install. Change the paths to the tools in `config.sh` to include the full path of the newly installed tools.

    - Example: `SPADES="path/to/spades"`

--------------------------------------------------------------------------------------------------------------------------------------------------------
## Tutorial
Example data to generate all of the figures for a simulated viral mock community that has differential relative abundance in response to host background have been included in `CWD/test_data/`. Simulated host resistance was based on the assumption that viruses for which a host is resistant will have a lower relative abundance of reads mapped.

1. Navigate to `CWD/test_data/sim_resist/mock_reads/` to look at the simulated illumina data. Create a `profile.txt` that includes the sample name only, excluding the suffixes for the illumina read pairs. One is already provided with the appropriate sample names in `CWD/`
2. There are two scripts that must be edited to match the suffix of the files to include for the illumina read pair information. This will be different for each run of ViCAT.
    - `CWD/scripts/run_spades.sh` must be changed at lines 19 and 20 to match he suffix of the reads. `_R1.fastq` and `_R2.fastq` for each respective variable.
    - `CWD/scripts/run_bowtie2.sh` must be changed at lines 18 and 19 to match the suffix of the reads. `_R1.fastq` and `_R2.fastq` for each respective variable.
3. Navigate back to the `CWD/` and edit the `config.sh` file to include the path to the raw data: `CWD/test_data/sim_resist/mock_reads/`. Also edit the path to the `profile.txt` file: `CWD/profile.txt`, and the path to the blast database.
4. The paths to tools should not need to be edited, if the `ViCAT.yml` was used to install the dependencies. If not, change the pahts to the tools.
5. OPTIONAL but ***suggested***: edit the `CWD/scripts/run_covgraph.sh` line 53-60 to remove portions of the sample name which may be uninformative and cause difficulties in reading coverage graph titles in the output files. The file can be found in `CWD/scripts`.
    - The `$samp_pref` variable should be changed to exclude uniformative information in the file name to only include the unique identifiers of the samples. This will allow for easier readability in the coverage graph titles. For this set, "begomovirus" is common to all of the samples so it should be removed at line 53. 
    - The `$samp_suff` variable should be changed in a similar manner at line 54. 
6. Run the `run.sh` file: `./run.sh`.
7. Once the run has finished, look in the `CWD/OUT_DIR/relative_abundance/` directory for output files to copy to `CWD/Community_Assembly/`. 
8. Change the `metadata.csv` file to include information about plant host resistance. The relavant data is found in `CWD/Community_Assembly/plant_genotype.csv`. It is worth noting that the order of the plant IDs may be different in the `metadata.csv` file than the **a priori** `plant_genotype.csv` file.
9. Run `Phyloseq_Community_Characterization.R` on terminal, or in R-Studio.
10. Files and graphs will be generated in the `CWD/Community_Assembly/` directory.
11. Ordination plots should discriminate plant samples by genotype of plant if run correctly.


## Usage
* Lines of code that can be changed are identified by obvious comment boxes in the respective scripts.*

1. Edit the `config.sh` file in the home directory to include paths to:
    - RAW: the path to the raw data for the run.
    - OUT_DIR: the path to the desired location for the outputs of the tool. This will be re-written each time the tool is used, so be careful to backup the results of previous runs.
    - BLASTDB: the location of the blast database to use.
    - BLAST1line: 1 line fasta to use for building reference index in bowtie2 mapping.
    - PROFILE: the path to the file that inlcudes the sample names as they relate to the filenames. An example profile is located in the home directory.
    - [TOOLS]: any paths to tools that were not installed via Anaconda.
    - [USER_INFO]: information for the HPC scheduler.
2. Edit the `run_spades.sh` lines 19 and 20 to set the file suffixes. The file can be found in `CWD/scripts`.

    `F1="${SMPLE}_R1.fastq"`
    
    `R1="${SMPLE}_R2.fastq"`
3. Edit the `run_bowtie2.sh` lines 18 and 19 to set the file suffixes. The file can be found in `CWD/scripts`.
    
4. OPTIONAL but ***suggested***: edit the `run_covgraph.sh` lines 53-54 to remove portions of the sample name which may be uninformative and cause difficulties in reading coverage graph titles in the output files. The file can be found in `CWD/scripts`.
5. OPTIONAL: edit the `run_blastn.sh` file on line 32 to change blast options.
6. OPTIONAL: edit the `run_write_files.sh` file on lines 25 and 42 to change filtering criteria for output files from Step 3. The file can be found in `CWD/scripts`. The file can be found in `CWD/scripts`.
7. OPTIONAL: edit the `run_cdhit.sh` file on line 20 to change percent similarity of blast hits for clustering. The file can be found in `CWD/scripts`.


To run the tool, go to the tools directory and execute the `run.sh` file.

`./run.sh`

#### Coverage Graphs
Coverage graphs for each sample can be found in `CWD/OUT_DIR/all_cov_graphs` or in the sample directory.

#### Virus Contigs
Virus contigs can be found in `CWD/OUT_DIR/virus_contigs` and a table for contig number, length, and coverage can be found here too.

#### Figure Generation

1. Once the tool has completed, output files for generating community assembly figures can be found in `CWD/OUT_DIR/relative_abundance/`.
    - `ViCAT_otutable.csv`: reads mapped to each taxa within each sample, for each sample.
    - `taxmat.csv`: Taxonomic matrix of identified species found in all samples, redundancy removed.
    - `metadata.csv`: Metadata matrix of sample information for the run. **THIS MUST BE EDITED WITH THE METADATA**
2. Copy the files to the `CWD/Community_Assembly/` directory.
3. If you have not already edited the metadata file, do so now.
4. Install the R packages found in `Phyloseq_Community_Characterization.R`, if you have not already installed them.
5. Change lines 201-205 of `Phyloseq_Community_Characterization.R` to the desired fields.
6. Run code as a block, or in R-Studio (v.2022.07.01 or higher).
7. Output files should be generated in `CWD/OUT_DIR/relative_abundance/`


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
    - `out/`: Directory of all standard outputs from the scripts. This will contain the outputs of the individual tools in the pipelines in directorys that have the names of the jobs in the `run.sh` file.
    - `err/`: Directory of all standard errors from the scripts. This will contain the errors of the individual tools in the pipelines in directorys that have the names of the jobs in the `run.sh` file.

### Explanation of directory structure

- `OUT_DIR/all_blast/`: Contains the copied blast results for all samples in the run.
- `OUT_DIR/all_consensus/`:Contains the extracted consensus sequences for all viral taxa in a sample, for all samples in the run.
- `OUT_DIR/all_cov/`: Contains all of the reads mapped to viral RefSeq genomes within a sample, for all samples in the run.
- `OUT_DIR/all_spades/`: Contains all of the *k*-mer and contig length filtered contigs for all samples in the run.
- `OUT_DIR/errors/`: Contains the error logs for each sample in the run.
- `OUT_DIR/relative_abundance/`: Contains the `combined_readscount.tsv`, `ViCAT_otutable.csv`, `taxmat.tsv`, and `metadata.csv` files for use in downstream community analysis in Phyloseq R package.
- `OUT_DIR/virus_contigs/`: Contains all of the contigs which have blast hits to viral RefSeq accessions. Also contains `virus_contig_table.csv`.
- `OUT_DIR/[sample(n)]/`: Directory for each sample in the run.
  - `blast/`: Contains the blast results, `[sample]_blastn_results.tsv`. Also contains intermediate files containing information on unique accessions and unique contigs in the sample for use in building bowtie2 reference indices and depositing virus contigs.
  - `blast_ref_fasta/`: Contains the RefSeq fasta genomes of the blast hits of the contigs of the sample, `[sample]_blasthits.fasta`. Also contains the fasta file with similar sequences removed by percent similarity, `[sample]_nr_ref.fasta`. Other files include the information used for clustering of sequences by similarity.
    - `bowtie2index/`: Contains the bowtie2 reference index for read mapping.
  - `bowtie2/`: Directory for files generated from read mapping to the RefSeq bowtie2 index.
    - `consensus/`: Contains the consensus fastq for the reads mapped to RefSeq bowtie2 index.
    - `unused_reads/`: Contains the reads not used in the read mapping to the RefSeq bowtie2 index.
    - `alignments/`: Contains the sorted binary alignment of reads to the reference index, `[sample]_vircom_nr_sorted.bam`. Also contains the reference split binary alignments for generating statistics and coverage graphs, `[sample]_vircom_nr_sorted.REF_NC_######.#.bam`.
      - `coverage/`: Contains coverage statistics of reads mapped to individual bases and combined coverage graphs, `[sample]_combined.svg`.
      - `stats/`: Contains statistics for the bowtie2 read mapping and files used in determining reads mapped to individual RefSeq accessions.


