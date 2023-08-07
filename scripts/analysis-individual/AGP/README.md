# AGP dataset
American Gut Project (_mSystems_, 2018) - [American Gut: an Open Platform for Citizen Science Microbiome Research][1]

[1]: https://journals.asm.org/doi/full/10.1128/mSystems.00031-18#DC1


## Choose samples of interest to download

The AGP dataset is a repository with more than 30,000 microbiome samples. We are only interested in samples from healthy and IBS individuals, so first we identified the list of samples to download in the [00_Metadata-AGP.R](00_Metadata-AGP.R) script.
The outputs of this R script are:
- a [list of runs](../../../data/analysis-individual/AGP/raw_fastq/list_files_agp.txt) (e.g. ERR2313945, ERR2092276, ...) to download with `wget` from the ENA;
- a metadata dataframe exported as a .csv file (in the [00_Metadata-AGP](../../../data/analysis-individual/AGP/00_Metadata-AGP/) directory, as "Metadata-AGP.csv").


## Download samples

To download samples, you can either:
- from your terminal, go into the [data_empty/analysis-individual/AGP/raw_fastq/](data_empty/analysis-individual/AGP/raw_fastq/) directory and follow the instructions of the `README.md` file to download samples from the ENA;
- obtain directly the fastq files from the `data/` folder we have deposited in the on Zenodo (ADD LINK)


## Preprocess fastq files into ASV and taxonomic tables

Preprocessing from raw fastq files, to ASV and taxonomic tables, is done with the [01_Dada2-AGP.Rmd](01_Dada2-AGP.Rmd) script. The output of that script is a phyloseq object, that contains an ASV table, taxonomic table, and metadata (can be found in the [phyloseq-without-phylotree](../../../data/phyloseq-objects/phyloseq-without-phylotree/) directory, within the "data" directory). Some general notes on the preprocessing are reported in the [01_Dada2-AGP_Notes.md](01_Dada2-AGP_Notes.md) file. An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.

To infer a phylogenetic tree, the [02_PhyloTree-AGP.R](02_PhyloTree-AGP.R) can be run on a server (or local computer, but will take a while). The input is the phyloseq object obtained from [01_Dada2-AGP.Rmd](01_Dada2-AGP.Rmd), and the output is a phyloseq object containing ASV+taxonomic+metadata tables and a phylogenetic tree (can be found in the [phyloseq-objects](../../../data/phyloseq-objects/) directory, within the "data" directory))


## Quick Exploratory Data Analysis (EDA)

For a quick EDA, you can take a look at [03_EDA-AGP.Rmd](03_EDA-AGP.Rmd). An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.