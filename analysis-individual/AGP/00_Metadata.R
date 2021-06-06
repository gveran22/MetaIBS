##########################
# Purpose: select IBS and healthy samples from American Gut Project fecal data
# Date: April 2021
# Author: Salomé Carcy
##########################



####################
# IMPORT LIBRARIES #
####################

library(tidyverse)

path <- "~/Projects/IBS_Meta-analysis_16S"



###################
# IMPORT METADATA #
###################

# This csv file was downloaded from the SRA
# with the accession number PRJEB11419
sradf <- read.csv(file.path(path, "data/analysis-individual/AGP/SraRunTable.csv"))

# This tsv file was downloaded from the figshare
# link shared by the paper from McDonald et al., 2018
paperdf <- read_tsv(file.path(path, "data/analysis-individual/AGP/correctedt2.tsv"))

# Compare two tables
dim(sradf) # 34,552 samples
dim(paperdf) # 17,854 samples

# The paper was published in 2018, so it is possible that more samples have been added on the SRA since
table(str_sub(paperdf$collection_date, start=-4)) # latest sample from paper metadata table is 2017
table(str_sub(sradf$ENA.FIRST.PUBLIC..run., end=4)) # latest sample from SRA metadata is 2021

# => We will extract sample Runs from SRA metadata


####################
# CLEANUP METADATA #
####################

# For simplicity purposes, we will keep only samples and columns of interest

# Keep only gut metadata
sradf <- sradf %>%
  filter(target_gene..exp. == "16S rRNA") %>% # 16s rRNA seq
  filter(SAMPLE_TYPE %in% c("feces", "Stool", "stool") & Organism %in% c("feces metagenome", "human gut metagenome", "gut metagenome")) # gut metagenome

# Reformat AGE & BMI columns
sradf <- sradf %>%
  # AGE
  mutate(host_age = as.integer(as.numeric(age_corrected))) %>%
  # WEIGHT
  mutate(weight_kg = as.numeric(weight_kg)) %>%
  filter(is.na(weight_kg) | (weight_kg > 30 & weight_kg < 200)) %>%
  # HEIGHT
  mutate(height_m = as.numeric(height_cm)*0.01) %>%
  filter(is.na(height_m) | (height_m > 1 & height_m < 2.5)) %>%
  # BMI
  mutate(host_bmi = weight_kg / (height_m)**2)

# Check that the "subset_healthy" corresponds to:
# - ppl between 20-69 yo
# - BMI between 18.5-30
# - no antibiotic in past year
# - no IBD
# - no diabetes
sradf %>%
  # SUBSET_HEALTHY
  filter(subset_healthy == TRUE) %>%
  
  # SAME SUBSET
  # filter(subset_age == TRUE) %>%
  # filter(bmi_cat %in% c("Normal", "Overweight")) %>%
  # filter(diabetes == "I do not have this condition") %>%
  # filter(ibd == "I do not have this condition") %>%
  # filter(antibiotic_history == "I have not taken antibiotics in the past year.") %>%
  
  # MANUAL SUBSET
  # filter(Age_years > 18 & Age_years < 70) %>%
  # filter(bmi >= 16 & bmi <= 35) %>%
  # filter(diabetes == "I do not have this condition") %>%
  # filter(ibd == "I do not have this condition") %>%
  # filter(antibiotic_history %in% c("I have not taken antibiotics in the past year.", "6 months")) %>%
  dim
# 8,379 people labeled as 'healthy'


# Subset SRADF
# - ppl between 18-69 yo
# - BMI between 18.5-30
# - no antibiotic in past year
# - no IBD, Cdiff, celiac disease, fungal overgrowth
# - no diabetes

subsetdf <- sradf %>%
  # Demographics
  filter(host_age >= 18 & host_age < 70) %>%
  filter(host_bmi >= 16 & host_bmi <= 35) %>%
  # Drugs taken
  filter(antibiotic_history %in% c("I have not taken antibiotics in the past year.", "6 months")) %>%
  # Gut comorbidities/differential diagnoses
  filter(ibd == "I do not have this condition") %>%
  filter(cdiff == "I do not have this condition") %>%
  filter(!gluten %in% unique(grep("diagnosed", sradf$gluten, value=TRUE))) %>%
  filter(fungal_overgrowth == "I do not have this condition") %>%
  filter(diabetes == "I do not have this condition")

# Check age_category & bmi_category
unique(subsetdf$age_cat)
unique(subsetdf$bmi_cat)
# There are some missing data in age & bmi categories, so we will re-define them
subsetdf <- subsetdf %>%
  mutate(age_cat = ifelse(host_age < 20, "teen",
                          ifelse(host_age %in% 20:29, "20s",
                          ifelse(host_age %in% 30:39, "30s",
                          ifelse(host_age %in% 40:49, "40s",
                          ifelse(host_age %in% 50:59, "50s",
                          ifelse(host_age %in% 60:69, "60s","NA"))))))) %>%
  mutate(bmi_cat = ifelse(host_bmi < 18.5, "Underweight",
                          ifelse(host_bmi >= 18.5 & host_bmi < 25.0, "Normal",
                          ifelse(host_bmi >= 25.0 & host_bmi < 30.0, "Overweight",
                          ifelse(host_bmi >= 30.0, "Obese", "NA")))))

# Keep only relevant columns for downstream analyses
keep_columns <- c("Run",
                  "pcr_primers..exp.",
                  # "Assay.Type",
                  # "Instrument",
                  "host_age",
                  "age_cat",
                  "host_bmi",
                  "bmi_cat",
                  "country",
                  "RACE",
                  "exercise_frequency",
                  "alcohol_frequency",
                  "age_cat",
                  "probiotic_frequency",
                  "gluten",
                  "bowel_movement_frequency",
                  "bowel_movement_quality",
                  "collection_season",
                  "sibo",
                  "lactose",
                  "subset_healthy",
                  "ibs",
                  "lung_disease", "liver_disease", "kidney_disease", "clinical_condition")
subsetdf <- subsetdf %>%
  select(all_of(keep_columns))


# Split healthy & IBS samples
healthyDF <- subsetdf %>%
  filter(subset_healthy == TRUE & 
           ibs == "I do not have this condition")

ibsDF <- subsetdf %>%
  filter(ibs == "Diagnosed by a medical professional (doctor, physician assistant)")

# We have too many healthy samples, so we will randomly select ~645 healthy samples
healthyDF <- healthyDF %>%
  sample_n(size=nrow(ibsDF))

# Verify there is an equal distribution in age & bmi
t.test(healthyDF$host_age, ibsDF$host_age)
t.test(healthyDF$host_bmi, ibsDF$host_bmi)


######################################### TRYING OUT #########################################




# Cleanup columns: bowel_movement_frequency, bowel_movement_quality








