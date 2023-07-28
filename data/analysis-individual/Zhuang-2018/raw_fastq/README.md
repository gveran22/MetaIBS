# raw_fastq (Zhuang)

We downloaded the `.fastq` files directly from the [ENA database](https://www.ebi.ac.uk/ena/browser/view/PRJNA475187), following these steps:
1. Downloading the table with the list of links to .fastq files: clicked on "TSV" next to "Download report";
2. Moving the downloaded `filereport_read_run_PRJNA475187_tsv.txt` file to this current directory;
3. Subsetting this table to the 29 samples of interest in the [00_Metadata-Zhuang.R](../../../../scripts/analysis-individual/Zhuang-2018/00_Metadata-Zhuang.R) script;
4. Exporting the list of `.fastq.gz` links to the `list_files_zhuang.txt`;
5. Executing the `download_fastq_zhuang.sh` bash script (which uses `wget`).