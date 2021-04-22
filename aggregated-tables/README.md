# Taxa aggregation

This directory contains tables of aggregated taxa for each dataset. The code used to obtain these .csv files is in the [aggregation.R](aggregation.R) file.

For each dataset, taxa were aggregated at all taxonomic levels (from Genus to Phylum). The counts correspond to the **absolute count** (no normalization of total count per sample!!). Tables contain the available covariates, too.

For Johannes: try the **Zhu-2019** dataset first!


# Description of datasets

|     Dataset    | N samples |          Covariates          | N phylum | N class | N order | N family | N genus |
| -------------- | :-------: | :--------------------------: | :------: | :-----: | :-----: | :------: | :-----: |
|  Fukui (2020)  |    109    |             -                |    12    |   18    |   45    |    69    |   207   |
| Hugerth (2019) |    520    |age, gender, BMI              |    16    |   23    |   58    |    92    |   252   |
|  Labus (2017)  |     52    |age, gender, BMI, IBS subtype |     6    |    9    |   21    |    34    |    91   |
|LoPresti (2019) |     73    |age, gender, IBS subtype      |     8    |   12    |   29    |    42    |    97   |
|  Nagel (2016)  |     30    |age, gender                   |    12    |   18    |   44    |    64    |   161   |
| Pozuelo (2015) |    290    |             -                |    15    |   24    |   65    |   105    |   312   |
|  Zeber (2016)  |     90    |gender, IBS subtype           |    11    |   20    |   49    |    76    |   215   |
|   Zhu (2019)   |     29    |age, gender                   |    11    |   17    |   40    |    57    |   141   |
| Zhuang (2018)  |     30    |IBS subtype (all IBS-D)       |     8    |   13    |   27    |    38    |    98   |


!!! **TO NOTE** !!!

Datasets with fecal samples and sigmoid mucosa biopsy samples, specified in the `sample_data` column:
- Hugerth (2019)
- LoPresti (2019)

Dataset with 2 samples per patient (taken with 1-month interval), specified in the `Collection` column:
- Pozuelo (2015)

Dataset with **only IBS-diarrhea** patients (other datasets have some patients with IBS-constipation, IBS-diarrhea, or IBS-mixed):
- Zhuang(2018)