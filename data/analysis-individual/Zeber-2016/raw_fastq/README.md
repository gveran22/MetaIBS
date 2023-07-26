# raw_fastq (Zeber-Lubecka)

We downloaded the `.fastq` files directly from the [ENA database](https://www.ebi.ac.uk/ena/browser/view/PRJEB11252), following these steps:
1. Downloading the table with the list of links to .fastq files: clicked on "TSV" next to "Download report";
2. Moving the downloaded `filereport_read_run_PRJEB11252_tsv.txt` file to this current directory;
3. Subsetting this table to the 90 samples of interest in the [00_Metadata-Zeber.R](../../../../scripts/analysis-individual/Zeber-2016/00_Metadata-Zeber.R) script;
4. Exporting the list of `.fastq.gz` links to the `list_files_zeber.txt`;
5. Executing the `download_fastq_zeber.sh` bash script (which uses `wget`).