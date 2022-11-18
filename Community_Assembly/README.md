### Explanation of directory
Export `taxmat.tsv`, `metadata.csv`, and `ViCAT_otutable.csv` from `OUT_DIR/relative_abundance/` to this directory and edit `Pyhloseq_Community_Characterization.R` to generate figures.

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
