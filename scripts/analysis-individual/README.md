# scripts/analysis-individual/

This repository contains `scripts` to preprocess raw fastq files into ASV and taxonomic tables. Our meta-analysis included 13 datasets, so there are 13 sub-directories. Even though we used a standardized pipeline to preprocess the data for all datasets, we still have separate directories & scripts as there are subtle differences between datasets (e.g. cleaning up metadata tables, different primer sequences, ...). However, we structured all these sub-directories in a similar way, and if you compare the `01_Dada2-NameDataset.Rmd` files across datasets, you will find few differences (e.g. different file paths, primer sequences, etc.).

<br/>

## Structure of each sub-directory

Within each sub-directory, you will find:
1. A `00_Metadata-NameDataset.R` script, which cleans up the SRA metadata dataframe to obtain a final metadata table containing only relevant covariates, and also exports a `.txt` file with the list of Runs to download.
2. A `download-NameDataset-samples/` folder, which was used to download raw fastq files from the SRA. This folder contains (1) a .txt file with the list of samples to download; and (2) a bash script to execute, that uses the SRA-toolkit to download the list of samples in the .txt file. Instructions to download the SRA-toolkit for any system can be found [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit), or we also put instructions below for mac users. Otherwise, we are also providing the raw `.fastq` files in our `data/` directory you downloaded from our Zenodo (ADD LINK).
3. A `01_Dada2-NameDataset.Rmd`, which preprocesses raw fastq files into ASV and taxonomic tables. It takes as input the raw fastq files downloaded from the SRA, and outputs a phyloseq object containing ASV, taxonomic, and metadata tables. These phyloseq objects can be found in the [phyloseq-without-phylotree](../../data/phyloseq-objects/phyloseq-without-phylotree/) directory, within the "data" directory. Also, there is a `01_Dada2-NameDataset_Notes.md` file with some notes on the preprocessing pipeline for each dataset (# samples, parameters used, anomalies encountered, # ASVs infered, # chimeras, etc.).
4. A `02_PhyloTree-NameDataset.R` script, which infers a phylogenetic tree and saves it in the phyloseq object (taken from the [Bioconductor Workflow](https://f1000research.com/articles/5-1492/v2)). We advise to run these scripts on a server if possible (would take a long time to run on local computer, especially for big datasets). They take as input the phyloseq objects with ASV+taxonomic+metadata tables (from step 2.), and output a phyloseq object that also contains a phylogenetic tree (can be found in the [phyloseq-objects](../../data/phyloseq-objects/) directory, within the "data" directory))
5. A `03_EDA-NameDataset.Rmd`, which performs standard exploratory data analyses (firmicutes/bacteroidota ratio, &beta;-diversity, etc.). These scripts were not used for figures in the paper.
6. HTML outputs of the `01_Dada2-NameDataset.Rmd` and `03_EDA-NameDataset.Rmd` R notebooks can be found in a `html_outputs` subdirectory, to compare your output with ours.

<br/>

## Requirements

### SRA toolkit
If you wish to download yourself the fastq files of each dataset from the SRA database (instead of using the raw fastq files we provide in our `data` directory), you can use the SRA toolkit. Instructions to download the SRA-toolkit for any system can be found [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit). For mac users, you can:
1. Go to your home directory
```cd ~```

2. Download the SRA toolkit with
```curl --output sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz```

3. Unzip the tar file:
```tar -vxzf sratoolkit.tar.gz```

4. Rename the SRA toolkit folder
```mv sratoolkit.3.0.1-mac64/ sratoolkit/```

5. Test that the SRA toolkit is functional: 
 ```fastq-dump --stdout -X 2 SRR390728```
You should see as output after a few seconds:
```
Read 2 spots for SRR390728
Written 2 spots for SRR390728
@SRR390728.1 1 length=72
CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCGATGGAGAATGACTTTGACAAGCTGAGAGAAGNTNC
+SRR390728.1 1 length=72
;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142;;;;;;;;;;;;;;;;;;;;;;;;;;;;;96&&&&(
@SRR390728.2 2 length=72
AAGTAGGTCTCGTCTGTGTTTTCTACGAGCTTGTGTTCCAGCTGACCCACTCCCTGGGTGGGGGGACTGGGT
+SRR390728.2 2 length=72
;;;;;;;;;;;;;;;;;4;;;;3;393.1+4&&5&&;;;;;;;;;;;;;;;;;;;;;<9;<;;;;;464262
```


<br/>

### Silva taxonomic reference
**Make sure to download Silva's reference fastas** and put it in the [silva-taxonomic-ref directory](../../data_empty/analysis-individual/CLUSTER/taxonomy/silva-taxonomic-ref/) to be able to do taxonomic alignment of infered ASVs:
- to reproduce results of the paper, download the `data/` folder from our Zenodo (ADD LINK), we saved the Silva v138 fastas in `data/analysis-individual/CLUSTER/taxonomy/silva-taxonomic-ref/`;
- to have the latest Silva version, download the reference fastas DADA2-formatted from [the DADA2 website](https://benjjneb.github.io/dada2/training.html) and put them in `data_empty/analysis-individual/CLUSTER/taxonomy/silva-taxonomic-ref/`.
