# Hugerth dataset
Hugerth (_Gut Microbiota_, 2019) - [No distinct microbiome signature of irritable bowel syndrome found in a Swedish random population][1]

[1]: https://gut.bmj.com/content/69/6/1076


## Choose samples of interest to download

The Hugerth dataset sequenced not only healthy/IBS samples, but also samples from individuals with other gut diseases. We thus first investigated the different metadata tables provided in [00_Metadata-Hugerth.Rmd](00_Metadata-Hugerth.Rmd). You can look at the html version of that script [here](./html_outputs/00_Metadata-Hugerth.html). The outputs of this R script are:
- a [list of runs](./download-Hugerth-samples/list_files_hugerth.txt) (e.g. ERR3548557, ERR358558, ...) to download with the SRA toolkit;
- a metadata dataframe exported as a .csv file (in the [00_Metadata-Hugerth](00_Metadata-Hugerth/) directory, as "Metadata-Hugerth.csv")


## Download samples

The raw fastq files are accessible on the SRA or ENA with the PRJEB31817 accession number. To download samples, you can either:
- from your terminal, go into the [download-Hugerth-samples](download-Hugerth-samples/) directory and execute the [download_fastq_hugerth.sh](download-Hugerth-samples/download_fastq_hugerth.sh) file;
- download directly the fastq files we have deposited on Zenodo (ADD LINK)


## Preprocess fastq files into ASV and taxonomic tables

Preprocessing from raw fastq files, to ASV and taxonomic tables, is done with the [01_Dada2-Hugerth.Rmd](01_Dada2-Hugerth.Rmd) script. The output of that script is a phyloseq object, that contains an ASV table, taxonomic table, and metadata (can be found in the [phyloseq-without-phylotree](../../../data/phyloseq-objects/phyloseq-without-phylotree/) directory, within the "data" directory). Some general notes on the preprocessing are reported in the [01_Dada2-Hugerth_Notes.md](01_Dada2-Hugerth_Notes.md) file. An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.

To infer a phylogenetic tree, the [02_PhyloTree-Hugerth.R](02_PhyloTree-Hugerth.R) can be run on a server (or local computer, but will take a while). The input is the phyloseq object obtained from [01_Dada2-Hugerth.Rmd](01_Dada2-Hugerth.Rmd), and the output is a phyloseq object containing ASV+taxonomic+metadata tables and a phylogenetic tree (can be found in the [phyloseq-objects](../../../data/phyloseq-objects/) directory, within the "data" directory))


## Quick Exploratory Data Analysis (EDA)

For a quick EDA, you can take a look at [03_EDA-Hugerth.Rmd](03_EDA-Hugerth.Rmd). An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.