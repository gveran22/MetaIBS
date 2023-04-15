# 04a_LogRatios-Taxa

This directory is linked to the [04a_LogRatios-Taxa.R](../../../scripts/analysis-combined/04a_LogRatios-Taxa.R) script, which allows to compute sample x log-ratio between taxa, at each taxonomic level (Genus -> Phylum). For instance, a matrix at the phylum level will look as such (don't mind the values, they are made up):

| Samples    | Abditibacteriota/Acidobacteriota | ... | Bacteroidota/Firmicutes | Bacteroidota/Proteobacteria | ... |
| ---------- | :-------: | :-------: | :-------: | :-------: | :-------: |
| ERR1051013 |     0     |    ...    |  -0.737   |  2.815    |    ...    |
| ERR1051015 |     0     |    ...    |   1.097   |  1.329    |    ...    |
| ERR1051017 | -0.87     |    ...    |   0.184   |  2.936    |    ...    |
|   ...      |    ...    |    ...    |    ...    |    ...    |    ...    |

This matrix will then be used to compute a UMAP with the [04b_UMAP.R](../../../scripts/analysis-combined/04b_UMAP.R) script.

<br/>

We recommend computing the log-ratios between all pairwise taxa on an HPC cluster. The bashscript `logratios-taxa.sh` can be used to submit a job on your cluster (you should adapt it to your HPC structure).

We have saved our matrices in the `pseudocounts_aft-agg` subdirectory (we added pseudocounts after aggregation to a taxonomic level, as doing it the other way around would cause the pseudocounts to sum up considerably during aggregation). If you are trying to reproduce our data, you can create an `output` folder (name given here as example), and make sure to change the name of your output folder (`path.output`) in the [04a_LogRatios-Taxa.R](../../../scripts/analysis-combined/04a_LogRatios-Taxa.R) script before executing it!