# raw_fastq (Liu)

We downloaded the `.fastq` files directly from the [ENA database](https://www.ebi.ac.uk/ena/browser/view/PRJNA544721), following these steps:
1. Downloading the table with the list of links to .fastq files: clicked on "TSV" next to "Download report";
2. Moving the downloaded `filereport_read_run_PRJNA544721_tsv.txt` file to this current directory;
3. Subsetting this table to the 128 samples of interest in the [00_Metadata-Liu.R](../../../../scripts/analysis-individual/Liu-2020/00_Metadata-Liu.R) script;
4. Exporting the list of `.fastq.gz` links to the `list_files_liu.txt`;
5. Executing the `download_fastq_liu.sh` bash script (which uses `wget`).