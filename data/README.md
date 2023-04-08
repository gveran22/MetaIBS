# Intro
This repository contains `data` obtained from preprocessing raw fastq files into ASV and taxonomic tables.

<br/>

# Structure of the repository
- [aggregated-tables](aggregated-tables/): `data` directory containing count tables aggregated at all taxonomic levels (Phylum, Class, Order, etc.) for each study;
- [asv-sequences](asv-sequences/): `data` directory containing all ASV sequences imputed from the [DADA2 algorithm](https://benjjneb.github.io/dada2/) for each dataset;
- [phyloseq-objects](phyloseq-objects/): `data` directory containing the phyloseq objects for each dataset, ready for analysis!

More details on how this data was obtained are provided in REAME.md files within each directory (hint: scripts to generate this data are in the [scripts/](../scripts/) directory).