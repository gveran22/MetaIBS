# Nagel dataset
Nagel et al. (_Microbiome_, 2016) - [Comparison of faecal microbiota in Blastocystis-positive and Blastocystisnegative irritable bowel syndrome patients][1]

[1]: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0191-0


## Download samples

The raw fastq files were obtained after private inquiry to the authors. We demultiplexed samples from the terminal using [sabre](https://github.com/najoshi/sabre). The raw (multiplexed) fastq files and the code to demultiplex them are reported in [data/.../Nagel-2016/raw_fastq_multiplexed/](../../../data/analysis-individual/Nagel-2016/raw_fastq_multiplexed/).

To download samples, you will need to download the `data/` directory we have deposited on Zenodo (ADD LINK). Within that directory, you will find:
- in `data/analysis-individual/Nagel-2016/raw_fastq_multiplexed/` the 2 multiplexed fastq files that authors provided us, that you can demultiplex with a bashscript we have provided within that same folder;
- in `data/analysis-individual/Nagel-2016/raw_fastq/` the 30 demultiplexed fastq files you can directly use in the [01_Dada2-Nagel.Rmd](01_Dada2-Nagel.Rmd) script.


## Preprocess fastq files into ASV and taxonomic tables

Preprocessing from raw fastq files, to ASV and taxonomic tables, is done with the [01_Dada2-Nagel.Rmd](01_Dada2-Nagel.Rmd) script. The output of that script is a phyloseq object, that contains an ASV table, taxonomic table, and metadata (can be found in the [phyloseq-without-phylotree](../../../data/phyloseq-objects/phyloseq-without-phylotree/) directory, within the "data" directory). Some general notes on the preprocessing are reported in the [01_Dada2-Nagel_Notes.md](01_Dada2-Nagel_Notes.md) file. An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.

To infer a phylogenetic tree, the [02_PhyloTree-Nagel.R](02_PhyloTree-Nagel.R) can be run on a server (or local computer, but will take a while). The input is the phyloseq object obtained from [01_Dada2-Nagel.Rmd](01_Dada2-Nagel.Rmd), and the output is a phyloseq object containing ASV+taxonomic+metadata tables and a phylogenetic tree (can be found in the [phyloseq-objects](../../../data/phyloseq-objects/) directory, within the "data" directory))


## Quick Exploratory Data Analysis (EDA)

For a quick EDA, you can take a look at [03_EDA-Nagel.Rmd](03_EDA-Nagel.Rmd). An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.