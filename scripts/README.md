# Intro
This repository contains `scripts` to (1) preprocess raw fastq files into ASV and taxonomic tables; and (2) perform statistical analyses on combined datasets. All the scripts to reproduce figures from the paper will be found within this directory.

<br/>

# Structure of the repository
- [analysis-combined](analysis-combined/): `code` directory containing R scripts for any analysis combining several datasets. These scripts were used to make the figures of the paper;
- [analysis-individual](analysis-individual/): `code` directory containing R scripts for all analysis on individual datasets (downloading the raw data, preprocessing raw data, obtaining phyloseq objects, performing standard exploratory data analysis).

More details are provided in REAME.md files within each directory.