# Meta-analysis of microbiome composition measured by 16S rRNA seq in IBS case-control studies

This is the repository for the study MetaIBS - large-scale amplicon-based meta analysis of irritable bowel syndrome (ADD PAPER LINK).
This repository contains all the data and code used for the meta-analysis, and can also be used as a template to extend this meta-analysis with your own data!

![gif](https://www.dana-farber.org/uploadedImages/Newsroom/Features/Gut_Instincts/microbiome-animated.gif)

<br/>

# Structure of the repository
- [data/](data/): `data` directory containing preprocessed data (phyloseq objects, count tables aggregated at different taxonomic levels, ASV sequences imputed from the raw fastq files)
- [scripts/](scripts/): `code` directory containing R scripts for (1) preprocessing raw fastq files into phyloseq objects; and (2) performing analyses on combined data from different datasets


<br/>

# How to use this repository
It depends what analysis you want to do (obviously). :thinking:

## Reproduce figures from the paper
In that case, you can:
1. Start by downloading & preprocessing the raw data for each dataset: everything you need to know will be in the [scripts/analysis-individual/](scripts/analysis-individual/) directory;
2. Repeat the analyses performed on several datasets: everything is in the [scripts/analysis-combined](scripts/analysis-combined/) directory.


## Extend the meta-analysis by adding your own datasets
In that case, you can:
1. Take a quick look at how raw fastq files were preprocessed for the datasets already included in this meta-analysis. You will want to look at the `01_Dada2-NameDataset.Rmd` files (or their HTML outputs) that are located in each subdirectory of [scripts/analysis-individual/](scripts/analysis-individual/). Even though there is one preprocessing script for each dataset, they are all mostly identical (little changes for the primer sequences, etc.), so you don't need to look at all the 13 scripts (for the 13 datasets);
2. Create a new directory within the [scripts/analysis-individual/](scripts/analysis-individual/) directory, to preprocess your own data. Save your phyloseq object into the [data/phyloseq-objects/](data/phyloseq-objects/) directory;
3. Go into the [scripts/analysis-combined/](scripts/analysis-combined/) directory, and you can re-run whichever analysis, adding your phyloseq object(s) at the beginning of the script to be included. Or you can create your own scripts of course to run additional analyses! :upside_down_face: