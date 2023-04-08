# Lo Presti dataset
Lo Presti et al. (_Frontiers in Microbiology_, 2019) - [Fecal and mucosal microbiota profiling in irritable bowel syndrome and inflammatory bowel disease][1]

[1]: https://www.frontiersin.org/articles/10.3389/fmicb.2019.01655/full


## Download samples

The raw fastq files are accessible on the SRA or ENA with the PRJNA391149 accession number. To download samples, you can either:
- from your terminal, go into the [download-LoPresti-samples](download-LoPresti-samples/) directory and execute the [download_fastq_lopresti.sh](download-LoPresti-samples/download_fastq_lopresti.sh) file;
- download directly the fastq files we have deposited on Zenodo (ADD LINK)


## Preprocess fastq files into ASV and taxonomic tables

Preprocessing from raw fastq files, to ASV and taxonomic tables, is done with the [01_Dada2-LoPresti.Rmd](01_Dada2-LoPresti.Rmd) script. The output of that script is a phyloseq object, that contains an ASV table, taxonomic table, and metadata (can be found in the [phyloseq-without-phylotree](../../../data/phyloseq-objects/phyloseq-without-phylotree/) directory, within the "data" directory). Some general notes on the preprocessing are reported in the [01_Dada2-LoPresti_Notes.md](01_Dada2-LoPresti_Notes.md) file. An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.

To infer a phylogenetic tree, the [02_PhyloTree-LoPresti.R](02_PhyloTree-LoPresti.R) can be run on a server (or local computer, but will take a while). The input is the phyloseq object obtained from [01_Dada2-LoPresti.Rmd](01_Dada2-LoPresti.Rmd), and the output is a phyloseq object containing ASV+taxonomic+metadata tables and a phylogenetic tree (can be found in the [phyloseq-objects](../../../data/phyloseq-objects/) directory, within the "data" directory))


## Quick Exploratory Data Analysis (EDA)

For a quick EDA, you can take a look at [03_EDA-LoPresti.Rmd](03_EDA-LoPresti.Rmd). An HTML output of this script can be found in the [html-outputs](./html-outputs/) directory.