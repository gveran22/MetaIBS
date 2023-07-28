# raw_fastq (Zhu)

We downloaded the `.fastq` files directly from the [ENA database](https://www.ebi.ac.uk/ena/browser/view/PRJNA566284), following these steps:
1. Downloading the table with the list of links to .fastq files: clicked on "TSV" next to "Download report";
2. Moving the downloaded `filereport_read_run_PRJNA566284_tsv.txt` file to this current directory;
3. Subsetting this table to the 29 samples of interest in the [00_Metadata-Zhu.R](../../../../scripts/analysis-individual/Zhu-2019/00_Metadata-Zhu.R) script;
4. Exporting the list of `.fastq.gz` links to the `list_files_zhu.txt`;
5. Executing the `download_fastq_zhu.sh` bash script (which uses `wget`).