# raw_fastq (Mars)

`.fastq` files can be directly downloaded from the [ENA database](https://www.ebi.ac.uk/ena/browser/view/PRJEB37924), following these steps:
1. Downloading the table with the list of links to .fastq files: click on "TSV" next to "Download report";
2. Moving the downloaded `filereport_read_run_PRJEB37924_tsv.txt` file to this current directory;
3. Subsetting this table to the 72 samples of interest in the [00_Metadata-Mars.Rmd](../../../../scripts/analysis-individual/Mars-2020/00_Metadata-Mars.Rmd) script;
4. Exporting the list of `.fastq.gz` links to the `list_files_mars.txt`;
5. Executing the `download_mars.sh` bash script (which uses `wget` to download the `.fastq.gz` files).