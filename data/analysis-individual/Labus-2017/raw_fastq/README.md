# raw_fastq (Labus)

Samples from the Labus dataset on the SRA have weird quality profile (and are 1200bp long...). We thus downloaded the `.fastq` files directly from the [ENA database](https://www.ebi.ac.uk/ena/browser/view/PRJNA373876) instead, and followed these steps:
1. Downloading the fastq files: clicked on "Download All" above "Generated FASTQ files:FTP";
2. Moved the `ena_files.zip` to the current directory
3. Unzipped `ena_files.zip` from our terminal by executing this `unzip_ena_fastq.sh` bash script:
```
chmod u+x unzip_ena_fastq.sh 
./unzip_ena_fastq.sh
```