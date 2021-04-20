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



## Details to know

- fukui:
	- no covariates (except healthy/IBS)

- hugerth:
  - stool samples & sigmoid mucosa biopsy samples (specified in column "sample_type"), sometimes from the same patient (patient identified by "host_ID")
  - covariates on age, BMI, gender, psychological state

- labus:
  - covariates on age, BMI, gender, IBS subtype (healthy control; ibs-constipation; ibs-diarrhea; ibs-mixed)

- lopresti:
  - stool samples & sigmoid mucosa biopsy samples (specified in column "sample_type"), sometimes from the same patient (patient identified by "host_ID")
  - covariates on age, gender, IBS subtype

- nagel:
  - covariates on age, gender

- pozuelo:
  - some patients provided a second stool sample after 1 month (specified in column "Collection", as "1st" and "2nd")
  - no covariates (except healthy/IBS)

- ringel:
  - NO METADATA AVAILABLE (not even healthy/ibs... we kept it only in case we could build a good predictive model, and use it on this dataset later on).

- zeber:
  - covariates on gender, IBS subtype

- zhu:
  - covariates on age, gender

- zhuang:
  - no covariates (except healthy/IBS)


