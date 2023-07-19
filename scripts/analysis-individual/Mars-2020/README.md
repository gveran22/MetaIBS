# Mars dataset
Mars et al. (_Cell_, 2020) - [Longitudinal Multi-omics Reveals Subset-Specific Mechanisms Underlying Irritable Bowel Syndrome][1]

[1]: https://www.sciencedirect.com/science/article/pii/S0092867420309983


## Download samples

The raw fastq files are accessible on the SRA or ENA with the PRJEB37924 accession number. To download samples, you can either:
- from your terminal, go into the [data_empty/analysis-individual/Mars-2020/raw_fastq/](data_empty/analysis-individual/Mars-2020/raw_fastq/) directory and follow the instructions of the `README.md` file to download samples from the ENA;
- obtain directly the fastq files from the `data/` folder we have deposited in the on Zenodo (ADD LINK)


## Preprocess fastq files into ASV and taxonomic tables

Preprocessing from raw fastq files, to ASV and taxonomic tables, is done with the [01_Dada2-Mars.Rmd](01_Dada2-Mars.Rmd) script. The output of that script is a phyloseq object, that contains an ASV table, taxonomic table, and metadata (can be found in the [phyloseq-without-phylotree](../../../data/phyloseq-objects/phyloseq-without-phylotree/) directory, within the "data" directory). Some general notes on the preprocessing are reported in the [01_Dada2-Mars_Notes.md](01_Dada2-Mars_Notes.md) file. An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.

To infer a phylogenetic tree, the [02_PhyloTree-Mars.R](02_PhyloTree-Mars.R) can be run on a server (or local computer, but will take a while). The input is the phyloseq object obtained from [01_Dada2-Mars.Rmd](01_Dada2-Mars.Rmd), and the output is a phyloseq object containing ASV+taxonomic+metadata tables and a phylogenetic tree (can be found in the [phyloseq-objects](../../../data/phyloseq-objects/) directory, within the "data" directory)).


## Quick Exploratory Data Analysis (EDA)

For a quick EDA, you can take a look at [03_EDA-Mars.Rmd](03_EDA-Mars.Rmd). An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.