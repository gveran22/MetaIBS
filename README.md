# Meta-analysis of microbiome composition measured by 16S rRNA seq in IBS case-control studies

This is the repository for the study MetaIBS - large-scale amplicon-based meta analysis of irritable bowel syndrome (add paper link when available).
This repository contains all the data and code used for the meta-analysis, and can also be used as a template to extend this meta-analysis with your own data!

<br/>

# Structure of the repository
- [aggregated-tables](aggregated-tables/): `data` directory containing count tables aggregated at all taxonomic levels (Phylum, Class, Order, etc.) for each study
- [analysis-combined](analysis-combined/): `code` directory containing R scripts for any analysis combining several datasets
- [analysis-individual](analysis-individual/): `code` directory containing R scripts for all analysis on individual datasets (downloading the raw data, preprocessing raw data, obtaining phyloseq objects, performing standard exploratory data analysis)
- [asv-sequences](asv-sequences/): `data` directory containing all ASV sequences imputed from the [DADA2 algorithm](https://benjjneb.github.io/dada2/) for each dataset
- [phyloseq-objects](phyloseq-objects/): `data` directory containing the phyloseq objects for each dataset, ready for analysis!

<br/>

# How to use this repository
It depends what analysis you want to do (obviously).

## 

## Reproduce figures from the paper
In that case, I woul