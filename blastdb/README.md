# Create blast database from fasta file

ViCAT can handle any custom database. This example uses the NCBI viral RefSeq sequences. 

## RefSeq Viral Database

# Download fasta files from GenBank

1. Go to viral [refseq](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/)
2. Download all `viral.#.#.genomic.fna.gz` files.
3. Upload to this directory, or the directory of all of your databases.
4. Unzip the files. 

# Concatenate files into one file and convert to 1 line fasta

1. Navigate to the directory with the database `.fna` 
2. Open `build_blastdb.sh`
3. Change the variables on lines 5-7 to the desired database names. For this example, leave as is.
4. Run `build_blastdb` in the directory with the `.fna` files from GenBank
5. Check that your database is accurate. There should be 9 files with different suffixes added after `.fsa`

## Custom database

# fasta files

1. Upload the appropriate fasta file, or files into the database directory.
2. Open `build_blastdb.sh`
3. Change the variables on lines 5-7 to the desired database names.
4. If there are multiple fasta files desired to be in the new database, change line 10 to have the file suffix matching the custom fasta files.
  - Example: `cat *.fasta > $ALL_FASTA` or `cat *.fa > $ALL_FASTA`
5. If there is only one fasta file, comment out line 10.
6. Run `build_blastdb` in the directory with the custom fasta database sequences.
7. Check that your database is accurate. There should be 9 files with different suffixes added after `.fsa`


