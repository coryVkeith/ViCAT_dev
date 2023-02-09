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
ViCAT_files.sh -o ../test_ViCAT
```


--------------------------------------------------------------------------------------------------------------------------------------------------------

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


