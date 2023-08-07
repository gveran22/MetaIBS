# raw_fastq (AGP)

We downloaded the `.fastq` files directly from the [ENA database](https://www.ebi.ac.uk/ena/browser/view/PRJEB11419), following these steps:
1. Downloading the table with the list of links to .fastq files: clicked on "TSV" next to "Download report";
2. Moving the downloaded `filereport_read_run_PRJEB11419_tsv.txt` file to this current directory;
3. Subsetting this table to the 1290 samples of interest in the [00_Metadata-AGP.R](../../../../scripts/analysis-individual/AGP/00_Metadata-AGP.R) script;
4. Exporting the list of `.fastq.gz` links to the `list_files_agp.txt`;
5. Executing the `download_fastq_agp.sh` bash script (which uses `wget`).