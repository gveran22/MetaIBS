# raw_fastq (LoPresti)

Samples from the LoPresti dataset on the SRA have weird quality profile (and are 800bp long...). We thus downloaded the `.fastq` files directly from the [ENA database](https://www.ebi.ac.uk/ena/browser/view/PRJNA391149) instead, and followed these steps:
1. Downloading the table with the list of links to .fastq files: clicked on "TSV" next to "Download report";
2. Moving the downloaded `filereport_read_run_PRJNA391149_tsv.txt` file to this current directory;
3. Subsetting this table to the 163 samples of interest in the [00_Metadata-LoPresti.R](../../../../scripts/analysis-individual/LoPresti-2019/00_Metadata-LoPresti.R) script;
4. Exporting the list of `.fastq.gz` links to the `list_files_lopresti.txt`;
5. Executing the `download_fastq_lopresti.sh` bash script (which uses `wget`).