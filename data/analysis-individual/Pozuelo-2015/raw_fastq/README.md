# raw_fastq (Pozuelo)

We downloaded the `.fastq` files directly from the [ENA database](https://www.ebi.ac.uk/ena/browser/view/PRJNA268708), following these steps:
1. Downloading the table with the list of links to .fastq files: clicked on "TSV" next to "Download report";
2. Moving the downloaded `filereport_read_run_PRJNA268708_tsv.txt` file to this current directory;
3. Subsetting this table to the 273 samples of interest in the [00_Metadata-Pozuelo.R](../../../../scripts/analysis-individual/Pozuelo-2015/00_Metadata-Pozuelo.R) script;
4. Exporting the list of `.fastq.gz` links to the `list_files_pozuelo.txt`;
5. Executing the `download_fastq_pozuelo.sh` bash script (which uses `wget`).