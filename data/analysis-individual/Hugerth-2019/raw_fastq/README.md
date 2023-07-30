# raw_fastq (Hugerth)

We downloaded the `.fastq` files directly from the [ENA database](https://www.ebi.ac.uk/ena/browser/view/PRJEB31817), following these steps:
1. Downloading the table with the list of links to .fastq files: clicked on "TSV" next to "Download report";
2. Moving the downloaded `filereport_read_run_PRJEB31817_tsv.txt` file to this current directory;
3. Subsetting this table to the 607 samples of interest in the [00_Metadata-Hugerth.Rmd](../../../../scripts/analysis-individual/Hugerth-2019/00_Metadata-Hugerth.Rmd) script;
4. Exporting the list of `.fastq.gz` links to the `list_files_hugerth.txt`;
5. Executing the `download_fastq_hugerth.sh` bash script (which uses `wget`).