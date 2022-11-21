# Create blast database from fasta file

ViCAT can handle any custom database. This example uses the NCBI viral RefSeq sequences. 

## NCBI Viral RefSeq Database

### Download fasta files from GenBank
Run `build_blastdb.sh` to download the viral refseq genomic databases, and create a blast database for the blast step in ViCAT. There should be 9 blast database files with different suffixes added after `.fsa`. Example: `viraldb.fsa.ndb`

## Custom database

1. Upload the appropriate fasta file, or files into the database directory.
2. Open `build_blastdb.sh`
3. Comment out lines 4 and 7, to prevent the NCBI viral RefSeq download.
4. Change the variables on lines 10-12 to the desired database names.
5. If there are multiple fasta files desired to be in the new database, change line 15 to have the file suffix matching the custom fasta files.
    - Example: `cat *.fasta > $ALL_FASTA` or `cat *.fa > $ALL_FASTA`
6. If there is only one fasta file, comment out line 10.
7. Run `build_blastdb` in the directory with the custom fasta database sequences.
8. Check that your database is accurate. There should be 9 files with different suffixes added after `.fsa`


