# scripts/analysis-combined/

This repository contains `scripts` to perform statistical analyses on combined datasets. These scripts were used to make the figures in the paper.


# Structure of this directory

|                                    Script                                |  Paper figure |                 Short description               |
| ------------------------------------------------------------------------ | :-----------: | :---------------------------------------------- |
| [`01_TaxonomicTree.R`](./01_TaxonomicTree.R)                             | Fig1A, FigS2  | Taxonomic tree of all ASVs inferred across datasets |
| [`02_LogRatio-FirmBact.R`](./02_LogRatio-FirmBact.R)                     | Fig2, FigS5   | Log-ratio of Firmicutes:Bacteroidota abundance in healthy vs IBS samples |
| [`03_Heatmaps.R`](./03_Heatmaps.R)                                       | Fig3A, FigS6A | Heatmap of microbial families relative abundances |
| [`04a_LogRatios-Taxa.R`](./04a_LogRatios-Taxa.R)                         | Fig3B, FigS6B | Compute log-ratio between all combinations of microbial families, and save the sample x log-ratio dataframe (to be provided as input for UMAP) |
| [`04b_UMAP.R`](./04b_UMAP.R)                                             | Fig3B, FigS6B | UMAP of log-ratios between microbial families across datasets |
| [`05_CommonASVs.R`](./05_CommonASVs.R)                                   | Fig5A         | Find how many ASVs are identical across datasets (expectation is to find common ASVs between datasets that amplified the same variable regions) |
| [`06_QCplot.R`](./06_QCplot.R)                                           | FigS1         | Plot number of reads per sample before/after quality filtering with DADA2 preprocessing |
| [`07_RelativAbund.R`](./07_RelativAbund.R)                               | FigS3         | Plot relative abundance of 5 main phyla across datasets |
| [`08_AlphaDiversity.R`](./08_AlphaDiversity.R)                           | FigS4         | Shannon and Simpson &alpha;-diversity indexes in healthy vs IBS samples |
| [`09_PCoA-BrayCurtis-BigDatasets.R`](./09_PCoA-BrayCurtis-BigDatasets.R) | FigS7         | Compute Bray-Curtis dissimilarity in AGP, Pozuelo and Hugerth datasets (3 biggest datasets) and perform PCoA |