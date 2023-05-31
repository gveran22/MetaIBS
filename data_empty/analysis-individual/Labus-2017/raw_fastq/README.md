# raw_fastq (Labus)

Samples from the Labus dataset on the SRA have weird quality profile (and are 1200bp long...). We recommend downloading the `.fastq` files directly from the [ENA database](https://www.ebi.ac.uk/ena/browser/view/PRJNA373876) instead, and following these steps:
1. Download the fastq files: click on "Download All" above "Generated FASTQ files:FTP";
2. Move the `ena_files.zip` to the current directory
3. Unzip `ena_files.zip` from your terminal by executing this `unzip_ena_fastq.sh` bash script:
```
chmod u+x unzip_ena_fastq.sh 
./unzip_ena_fastq.sh
```

<br/>

As an (easier) alternative, you can also get the fastq files directly from the `data/` folder from our Zenodo (ADD LINK).