# data/analysis-combined/

This folder will auto-populate when you run the scripts in [scripts/analysis-combined/](../scripts/analysis-combined). These scripts were used to make the figures of the paper. When you run those scripts, you will be able to reproduce the paper figures, and you will be able to save them as `.jpeg` files in each sub-directory.

As a reminder, here are the scripts \& the corresponding figures they generate:

|                                **Script**                                |**Paper figure(s)**|               **Short description**         |
| ------------------------------------------------------------------------ | :-----------: | ---------------------------------------------- |
| [`01_TaxonomicTree.R`](../scripts/analysis-combined/01_TaxonomicTree.R)                             | Fig1A, FigS2  | Taxonomic tree of all ASVs inferred across datasets |
| [`02_LogRatio-FirmBact.R`](../scripts/analysis-combined/02_LogRatio-FirmBact.R)                     | Fig2, FigS5   | Log-ratio of Firmicutes:Bacteroidota abundance in healthy vs IBS samples |
| [`03_Heatmaps.R`](../scripts/analysis-combined/03_Heatmaps.R)                                       | Fig3A, FigS6A | Heatmap of microbial families relative abundances |
| [`04a_LogRatios-Taxa.R`](../scripts/analysis-combined/04a_LogRatios-Taxa.R)                         | Fig3B, FigS6B | Compute log-ratio between all combinations of microbial families, and save the sample x log-ratio dataframe (to be provided as input for UMAP) |
| [`04b_UMAP.R`](../scripts/analysis-combined/04b_UMAP.R)                                             | Fig3B, FigS6B | UMAP of log-ratios between microbial families across datasets |
| [`05_Common-ASVs.R`](../scripts/analysis-combined/05_Common-ASVs.R)                                   | Fig5A         | Find how many ASVs are identical across datasets (expectation is to find common ASVs between datasets that amplified the same variable regions) |
| [`06_QCplot.R`](../scripts/analysis-combined/06_QCplot.R)                                           | FigS1         | Plot number of reads per sample before/after quality filtering with DADA2 preprocessing |
| [`07_RelativAbund.R`](../scripts/analysis-combined/07_RelativAbund.R)                               | FigS3         | Plot relative abundance of 5 main phyla across datasets |
| [`08_AlphaDiversity.R`](../scripts/analysis-combined/08_AlphaDiversity.R)                           | FigS4         | Shannon and Simpson &alpha;-diversity indexes in healthy vs IBS samples |
| [`09_PCoA-BrayCurtis-BigDatasets.R`](../scripts/analysis-combined/09_PCoA-BrayCurtis-BigDatasets.R) | FigS7         | Compute Bray-Curtis dissimilarity in AGP, Pozuelo and Hugerth datasets (3 biggest datasets) and perform PCoA |
