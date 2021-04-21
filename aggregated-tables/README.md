# Taxa aggregation

This directory contains tables of aggregated taxa for each dataset. The code used to obtain these .csv files is in the [aggregation.R](aggregation.R) file.

For each dataset, taxa were aggregated at all taxonomic levels (from Genus to Phylum). The counts correspond to the **absolute count** (no normalization of total count per sample!!). Tables contain the available covariates, too.

For Johannes: try the **Zhu-2019** dataset first!


# Description of datasets

|   Dataset  | N samples |  Covariates | N phyla | N classes | N orders | N families | N genera |
| ---------- | :-------: | :---------: | :-----: | :-------: | :------: | :--------: | :------: |
| Zhu (2019) |    29     | age, gender |    11   |    17     |    40    |     57     |   141    |