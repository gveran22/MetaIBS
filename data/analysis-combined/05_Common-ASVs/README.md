# Common ASVs

When using the `merge_phyloseq()` function, ASVs that are **exactly** identical between datasets were automatically merged (see [05_Common-ASVs.R](../../../scripts/analysis-combined/05_Common-ASVs.R) script).
This plot shows the number of common ASVs that were identified between datasets:

<img src="./commonASV_merge-phyloseq-funct.jpg" alt="show common ASVs" width="400"/>

The `phyloseq` objects and `dataframes`  saved in this folder as `.rds` and `.csv` files (respectively) contain only the common ASVs identified between Nagel-Pozuelo, Liu-Zhuang, Hugerth-Zhu and LoPresti-Ringel.

In the `dataframes`, the "OTU" column contains the nucleotide sequence of these common ASVs. If you check the number of unique "OTU" sequence in each dataframe, you should find the same number shown on the plot (i.e. 806 common ASVs for Nagel-Pozuelo). The counts in the "Abundance" column correspond to the absolute counts (no normalization done) of these "common ASVs". 

Careful, it doesn't correspond to _all_ the counts from each dataset, but only the counts of the shared ASVs, so like this:
- Nagel-Pozuelo: shared ASVs represent ~80% of total counts
- Liu-Zhuang: shared ASVs represent ~81% of total counts
- Hugerth-Zhu: shared ASVs represent ~58% of total counts
- LoPresti-Ringel: shared ASVs represent ~0.6% of total counts