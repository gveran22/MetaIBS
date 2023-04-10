# 06_QCplot

This directory will contain the number of reads per sample before and after the DADA2 pipeline (dataframe saved as an `.rds` file for each dataset separately). These `.rds` files are generated when running the `01_Dada2-NameDataset.Rmd` files in [analysis-individual](../../../scripts/analysis-individual/).
- Before: number of reads in the raw .fastq files
- After: in the ASV table generated from the DADA2 scripts (found in each subdirectory of [analysis-individual](../../../scripts/analysis-individual/)). Eukaryotic or unknown ASVs have not been removed yet.

Afterwards, these dataframes will be used in [06_QCplot.R](../../../scripts/analysis-combined/06_QCplot.R) to make a plot of the number of reads per sample before/after DADA2 pipeline (figure S1).
