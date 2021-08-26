# Taxa aggregation

This directory contains tables of aggregated taxa for each dataset. The code used to obtain these .csv files is in the [aggregation.R](aggregation.R) file.

For each dataset, taxa were aggregated at all taxonomic levels (from Genus to Phylum). The counts correspond to the **absolute count** (no normalization of total count per sample!!). Tables contain the available covariates, too.


# Description of datasets

|     Dataset    | N samples |               Covariates          | N phylum | N class | N order | N family | N genus |
| -------------- | :-------: | :--------------------------------:| :------: | :-----: | :-----: | :------: | :-----: |
|   AGP (2021)   |  1,183    |age, BMI, bowel mvt, comorbidities |    31    |   69    |  141    |   225    |   619   |
|  Fukui (2020)  |    110    |                  -                |    12    |   18    |   45    |    70    |   210   |
| Hugerth (2019) |    525    |age, gender, BMI                   |    18    |   25    |   64    |   100    |   263   |
|  Labus (2017)  |     52    |age, gender, BMI, IBS subtype      |     6    |    9    |   21    |    34    |    91   |
|   Liu (2020)   |    128    |age, gender, BMI, Bristol & IBS-SSS|    44    |   96    |  217    |   315    |   697   |
|LoPresti (2019) |     57    |age, gender, IBS subtype           |     8    |   12    |   29    |    40    |    90   |
|   Mars (2020)  |     69    |age, gender, BMI, IBS subtype      |    18    |   29    |   66    |   108    |   269   |
|  Nagel (2016)  |     30    |age, gender,IBS subtype (all IBS-D)|    12    |   18    |   44    |    64    |   163   |
| Pozuelo (2015) |    273    |IBS subtype, 2 collection times    |    15    |   25    |   64    |   103    |   308   |
|  Zeber (2016)  |     90    |gender, IBS subtype                |    11    |   19    |   49    |    80    |   236   |
|   Zhu (2019)   |     29    |age, gender                        |    11    |   17    |   39    |    57    |   142   |
| Zhuang (2018)  |     30    |IBS subtype (all IBS-D)            |     9    |   14    |   32    |    44    |   116   |


!!! **TO NOTE** !!!

Datasets with fecal samples and sigmoid mucosa biopsy samples, specified in the `sample_data` column:
- Hugerth (2019)
- LoPresti (2019)
- Mars (2020): only sigmoid biopsy samples!

Dataset with 2 samples per patient, specified in the `Collection` column:
- Pozuelo (2015) - 1 month interval
- Mars (2020) - 6 months interval

Dataset with **only IBS-diarrhea** patients (other datasets have some patients with IBS-constipation, IBS-diarrhea, or IBS-mixed):
- Nagel (2016)
- Zhuang(2018)
- Liu (2020)