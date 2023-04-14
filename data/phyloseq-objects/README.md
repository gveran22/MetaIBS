# Phyloseq objects

This directory contains all the phyloseq objects saved as `.rds` files, obtained from preprocessing the raw data for each dataset (refer to the [scripts/analysis-individual/](../../scripts/analysis-individual/) folder for more information on how the data was preprocessed).
To better understand how phyloseq objects are structured, you can refer to the [phyloseq tutorial](https://joey711.github.io/phyloseq/import-data.html).

Phyloseq objects in this directory contain:
- an ASV table (`otu_table()`);
- a taxonomic table (`tax_table()`);
- a metadata table (`sample_data()`);
- a phylogenetic tree (`phy_tree()`);
- a DNAStringSet (`refseq()`) that contains the nucleotide sequence for each ASV (e.g. ASV1 - ATTT...TTGA / ASV2 - TTAG...GGC / etc.).

Phyloseq objects in the [phyloseq-without-phylotree/](./phyloseq-without-phylotree/) subdirectory contain:
- an ASV table (`otu_table()`);
- a taxonomic table (`tax_table()`);
- a metadata table (`sample_data()`).
These phyloseq objects without a phylogenetic tree are useful for analyses on combined datasets, as you can't merge phylogenetic trees...


## Datasets description

|     Dataset    | N samples |                  Covariates               |                        TO BE CAREFUL                           |
| -------------- | :-------: | :---------------------------------------- | :------------------------------------------------------------- |
|  `AGP (2021)`  |  1,183    |age, BMI, bowel mvt, comorbidities         |                    -                                           |
| `Fukui (2020)` |    110    |                    -                      |                    -                                           |
|`Hugerth (2019)`|    525    |age, gender, BMI, psychology               | 2 `sample_type`: stool, sigmoid mucosa                         |
| `Labus (2017)` |     52    |age, gender, BMI, IBS subtype              |                    -                                           |
|  `Liu (2020)`  |    128    |age, gender, BMI, IBS subtype (all IBS-D)  | all IBS patients are IBS-diarrhea                              |
|`LoPresti (2019)`|    73    |age, gender, IBS subtype                   | 2 `sample_type`: stool, sigmoid mucosa                         |
|  `Mars (2020)` |     69    |age, gender, BMI, IBS subtype              |2 sample `Collection` per patient & only sigmoid mucosa samples |
| `Nagel (2016)` |     30    |age, gender, IBS subtype (all IBS-D)       | all IBS patients are IBS-diarrhea                              |
|`Pozuelo (2015)`|    273    |IBS subtype                                | 2 sample `Collection` per patient                              |
|`Ringel (2015)` |     76    | _No label IBS/Healthy_                    |                    -                                           |
| `Zeber (2016)` |     90    |gender, IBS subtype                        |                    -                                           |
|  `Zhu (2019)`  |     29    |age, gender                                |                    -                                           |
|`Zhuang (2018)` |     30    |IBS subtype (all IBS-D)                    | all IBS patients are IBS-diarrhea                              |
