# Zhu dataset
Zhu (_Frontiers in Cellular and Infection Microbiology_, 2019) - [Identification of Gut Microbiota and Metabolites Signature in Patients With Irritable Bowel Syndrome][1]

[1]: https://www.frontiersin.org/articles/10.3389/fcimb.2019.00346/full 


## Download samples

The raw fastq files are accessible on the SRA or ENA with the PRJNA566284 accession number. To download samples, you can either:
- from your terminal, go into the [download-Zhu-samples](download-Zhu-samples/) directory and execute the [download_fastq_zhu.sh](download-Zhu-samples/download_fastq_zhu.sh) file;
- download directly the fastq files we have deposited on Zenodo (ADD LINK)


## Preprocess fastq files into ASV and taxonomic tables

Preprocessing from raw fastq files, to ASV and taxonomic tables, is done with the [01_Dada2-Zhu.Rmd](01_Dada2-Zhu.Rmd) script. The output of that script is a phyloseq object, that contains an ASV table, taxonomic table, and metadata (can be found in the [phyloseq-without-phylotree](../../../data/phyloseq-objects/phyloseq-without-phylotree/) directory, within the "data" directory). Some general notes on the preprocessing are reported in the [01_Dada2-Zhu_Notes.md](01_Dada2-Zhu_Notes.md) file. An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.

To infer a phylogenetic tree, the [02_PhyloTree-Zhu.R](02_PhyloTree-Zhu.R) can be run on a server (or local computer, but will take a while). The input is the phyloseq object obtained from [01_Dada2-Zhu.Rmd](01_Dada2-Zhu.Rmd), and the output is a phyloseq object containing ASV+taxonomic+metadata tables and a phylogenetic tree (can be found in the [phyloseq-objects](../../../data/phyloseq-objects/) directory, within the "data" directory)).


## Quick Exploratory Data Analysis (EDA)

For a quick EDA, you can take a look at [03_EDA-Zhu.Rmd](03_EDA-Zhu.Rmd). An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.