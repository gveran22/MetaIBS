# Phyloseq objects

## My intuition on which datasets will work best on trac

**Try the physeq_zhu.rds first!** I would get >90% accuracy with a logistic regression or random forest.

Then, other datasets that might work well:
- labus
- pozuelo
- zhuang
- nagel

Datasets for which I have little hope:
- hugerth
- lopresti
- zeber



## Datasets description

|     Dataset    | N samples |          Covariates          |             TO BE CAREFUL              |
| -------------- | :-------: | :--------------------------- | :------------------------------------- |
|  Fukui (2020)  |    109    |             -                |                    -                   |
| Hugerth (2019) |    520    |age, gender, BMI, psychology  | 2 `sample_type`: stool, sigmoid mucosa |
|  Labus (2017)  |     52    |age, gender, BMI, IBS subtype |                    -                   |
|LoPresti (2019) |     73    |age, gender, IBS subtype      | 2 `sample_type`: stool, sigmoid mucosa |
|  Nagel (2016)  |     30    |age, gender                   |                    -                   |
| Pozuelo (2015) |    290    |             -                | 2 sample `Collection` per patient      |
| Ringel (2015)  |     76    | _No label IBS/Healthy_       |                    -                   |
|  Zeber (2016)  |     90    |gender, IBS subtype           |                    -                   |
|   Zhu (2019)   |     29    |age, gender                   |                    -                   |
| Zhuang (2018)  |     30    |IBS subtype (all IBS-D)       | all IBS patients are IBS-diarrhea      |


