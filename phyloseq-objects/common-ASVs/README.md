# Common ASVs

When using the `merge_phyloseq()` function, ASVs that are **exactly** identical between datasets were automatically merged.
This plot shows the number of common ASVs that were identified between datasets:

<img src="./commonASV_merge-phyloseq-funct.jpg" alt="show common ASVs" width="400"/>

The phyloseq objects saved in this folder contain only the common ASVs identified between Nagel-Pozuelo, Liu-Zhuang, Hugerth-Zhu and LoPresti-Ringel.
Careful, it doesn't correspond to all the counts from each dataset, but only the counts of the shared ASVs, so like this:
- Nagel-Pozuelo: shared ASVs represent ~80% of total counts
- Liu-Zhuang: shared ASVs represent ~81% of total counts
- Hugerth-Zhu: shared ASVs represent ~58% of total counts
- LoPresti-Ringel: shared ASVs represent ~0.6% of total counts