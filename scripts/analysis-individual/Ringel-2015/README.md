# Ringel-Kulka dataset
Ringel-Kulka et al. (_Am J Physiol Gastrointest Liver Physiol_, 2016) - [Molecular characterization of the intestinal microbiota in patients with and without abdominal bloating][1]

[1]: https://journals.physiology.org/doi/full/10.1152/ajpgi.00044.2015


## Download samples

The raw fastq files are accessible on the SRA or ENA with the SRP066323 accession number. To download samples, you can either:
- from your terminal, go into the [download-Ringel-samples](download-Ringel-samples/) directory and execute the [download_fastq_ringel.sh](download-Ringel-samples/download_fastq_ringel.sh) file;
- download directly the fastq files we have deposited on Zenodo (ADD LINK)


## Preprocess fastq files into ASV and taxonomic tables

Preprocessing from raw fastq files, to ASV and taxonomic tables, is done with the [01_Dada2-Ringel.Rmd](01_Dada2-Ringel.Rmd) script. The output of that script is a phyloseq object, that contains an ASV table, taxonomic table, and metadata (can be found in the [phyloseq-without-phylotree](../../../data/phyloseq-objects/phyloseq-without-phylotree/) directory, within the "data" directory). Some general notes on the preprocessing are reported in the [01_Dada2-Ringel_Notes.md](01_Dada2-Ringel_Notes.md) file. An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.

To infer a phylogenetic tree, the [02_PhyloTree-Ringel.R](02_PhyloTree-Ringel.R) can be run on a server (or local computer, but will take a while). The input is the phyloseq object obtained from [01_Dada2-Ringel.Rmd](01_Dada2-Ringel.Rmd), and the output is a phyloseq object containing ASV+taxonomic+metadata tables and a phylogenetic tree (can be found in the [phyloseq-objects](../../../data/phyloseq-objects/) directory, within the "data" directory)).


## Quick Exploratory Data Analysis (EDA)

For a quick EDA, you can take a look at [03_EDA-Ringel.Rmd](03_EDA-Ringel.Rmd). An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.