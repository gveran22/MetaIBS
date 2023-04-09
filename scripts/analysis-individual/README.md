# scripts/analysis-individual/

This repository contains `scripts` to preprocess raw fastq files into ASV and taxonomic tables. Our meta-analysis included 13 datasets, so there are 13 sub-directories. Even though we used a standardized pipeline to preprocess the data for all datasets, we still have separate directories & scripts as there are subtle differences between datasets (e.g. cleaning up metadata tables, different primer sequences, ...). However, we structured all these sub-directories in a similar way, and if you compare the `01_Dada2-NameDataset.Rmd` files across datasets, you will find few differences (e.g. different file paths, primer sequences, etc.).


# Structure of each sub-directory

Within each sub-directory, you will find:
1. A `download-NameDataset-samples/` folder, which was used to download raw fastq files from the SRA. This folder contains (1) a .txt file with the list of samples to download; and (2) a bash script to execute, that uses the SRA-toolkit to download the list of samples in the .txt file. Instructions to download the SRA-toolkit can be found [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit). Otherwise, you can also download the fastq files from our Zenodo repo (ADD LINK).
2. A `01_Dada2-NameDataset.Rmd`, which preprocesses raw fastq files into ASV and taxonomic tables. It takes as input the raw fastq files downloaded from the SRA, and outputs a phyloseq object containing ASV, taxonomic, and metadata tables. These phyloseq objects can be found in the [phyloseq-without-phylotree](../../data/phyloseq-objects/phyloseq-without-phylotree/) directory, within the "data" directory. Also, there is a `01_Dada2-NameDataset_Notes.md` file with some notes on the preprocessing pipeline for each dataset (# samples, parameters used, anomalies encountered, # ASVs infered, # chimeras, etc.).
3. A `02_PhyloTree-NameDataset.R` script, which infers a phylogenetic tree and saves it in the phyloseq object (taken from the [Bioconductor Workflow](https://f1000research.com/articles/5-1492/v2)). We advise to run these scripts on a server if possible (would take a long time to run on local computer, especially for big datasets). They take as input the phyloseq objects with ASV+taxonomic+metadata tables (from step 2.), and output a phyloseq object that also contains a phylogenetic tree (can be found in the [phyloseq-objects](../../data/phyloseq-objects/) directory, within the "data" directory))
4. A `03_EDA-NameDataset.Rmd`, which performs standard exploratory data analyses (firmicutes/bacteroidota ratio, &beta;-diversity, etc.). These scripts were not used for figures in the paper.
5. HTML outputs of the `01_Dada2-NameDataset.Rmd` and `03_EDA-NameDataset.Rmd` R notebooks can be found in a `html_outputs` subdirectory, to compare your output with ours.


# Important links
**Make sure to download Silva's reference fastas** and put it in the [silva-taxonomic-ref directory](../../data/preprocessing/silva-taxonomic-ref/) to be able to do taxonomic alignment of infered ASVs:
- to reproduce results of the paper, download the fasta files from our zenodo link (ADD LINK), where we used Silva v138
- to have the latest Silva version, download the reference fastas DADA2-formatted from [the DADA2 website](https://benjjneb.github.io/dada2/training.html)