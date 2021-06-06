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
dim(sradf) # 33,610 samples
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

# Reformat age & bmi columns
sradf <- sradf %>%
  # AGE
  mutate(Age_years = as.numeric(Age_years)) %>%
  # WEIGHT
  mutate(weight_kg = as.numeric(weight_kg)) %>%
  filter(is.na(weight_kg) | (weight_kg > 30 & weight_kg < 200)) %>%
  # HEIGHT
  mutate(height_m = as.numeric(height_cm)*0.01) %>%
  filter(is.na(height_m) | (height_m > 1 & height_m < 2.5)) %>%
  # BMI
  mutate(bmi = weight_kg / (height_m)**2)

# Check that the "subset_healthy" corresponds to:
# - ppl between 20-69 yo
# - BMI between 18.5-30
# - no antibiotic in past year
# - no IBD
# - no diabetes
sradf %>%
  # SUBSET_HEALTHY
  # filter(subset_healthy == TRUE) %>%
  
  # SAME SUBSET
  # filter(subset_age == TRUE) %>%
  # filter(bmi_cat %in% c("Normal", "Overweight")) %>%
  # filter(diabetes == "I do not have this condition") %>%
  # filter(ibd == "I do not have this condition") %>%
  # filter(antibiotic_history == "I have not taken antibiotics in the past year.") %>%
  
  # MANUAL SUBSET
  filter(Age_years > 18 & Age_years < 70) %>%
  filter(bmi >= 16 & bmi <= 35) %>%
  filter(diabetes == "I do not have this condition") %>%
  filter(ibd == "I do not have this condition") %>%
  filter(antibiotic_history %in% c("I have not taken antibiotics in the past year.", "6 months")) %>%
  dim
# 8,353 people labeled as 'healthy'


# Subset SRADF
# - ppl between 18-69 yo
# - BMI between 18.5-30
# - no antibiotic in past year
# - no IBD, Cdiff, celiac disease, fungal overgrowth
# - no diabetes

subsetdf <- sradf %>%
  # Demographics
  filter(Age_years > 18 & Age_years < 70) %>%
  filter(bmi >= 16 & bmi <= 35) %>%
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
test <- subsetdf %>%
  mutate(age_cat = ifelse(Age_years < 20, "teen",
                          ifelse(Age_years %in% 20:29.9, "20s",
                          ifelse(Age_years %in% 30:39.9, "30s",
                          ifelse(Age_years %in% 40:49.9, "40s",
                          ifelse(Age_years %in% 50:59.9, "50s",
                          ifelse(Age_years %in% 60:70, "60s","NA")))))))

# Split healthy & IBS samples
healthyDF <- subsetdf %>%
  filter(subset_healthy == TRUE & 
           ibs == "I do not have this condition")

ibsDF <- subsetdf %>%
  filter(ibs == "Diagnosed by a medical professional (doctor, physician assistant)")

# We have too many healthy samples, so we will randomly select ~630 healthy samples




######################################### TRYING OUT #########################################


# Without any filter, there are 2153 IBS patients
# + no IBD + no C.diff infection + no celiac disease/gluten allergy + no ATB in the past 6 months : 1116

# Additional optional filters:
# + no fungal_overgrowth => 819 samples
# ++ age between 20 and 69 => 722 samples
# +++ BMI between 18-30 => 590 samples
# ++++ remove other comorbidities (diabetes, lung disease, etc.) => 331 samples



# Remove columns that will not be of use
skip_columns <- c("SAMPLE_TYPE", "Organism",
                  "BioProject",
                  "center_name..exp.",
                  "Center.Name",
                  "Consent",
                  "ENA.FIRST.PUBLIC..run.",
                  "LibrarySelection",
                  "LibrarySource",
                  "Platform",
                  "ReleaseDate",
                  "sequencing_meth..exp.",
                  "SRA.Study",
                  "well_description..exp.",
                  "run_prefix..exp.",
                  "run_date..exp.",
                  "Bytes",
                  "DATASTORE.provider",
                  "DATASTORE.region",
                  "Title",
                  "sample_plate..exp.",
                  "target_gene..exp.",
                  "target_subfragment..exp.",
                  "exercise_location",
                  "cosmetics_frequency",
                  "pool_frequency",
                  "last_travel",
                  "Elevation",
                  "contraceptive",
                  "LATITUDE",
                  "LONGITUDE",
                  "fermented_plant_frequency",
                  "dominant_hand",
                  "livingwith",
                  "flu_vaccine_date",
                  "types_of_plants",
                  "consume_animal_products_abx",
                  "last_move",
                  "age_corrected",
                  "primer_plate..exp.",
                  "well_id..exp.",
                  "roommates_in_study",
                  "collection_time")
sradf <- sradf %>%
  select(-all_of(skip_columns))



# Cleanup columns: bowel_movement_frequency, bowel_movement_quality




###############
# IBS SAMPLES #
###############

# Subset table to individuals with IBS medical diagnosis
table(sradf$ibs) # 2895 "Diagnosed by a medical professional (doctor, physician assistant)"
ibs.df <- sradf[grep("doctor", sradf$ibs),]
dim(ibs.df) # 2895 patients
table(ibs.df$bowel_movement_quality) 

# Check it is only fecal samples
table(ibs.df$Organism) # 2671 human gut metagenome
table(ibs.df$SAMPLE_TYPE) # 2383 stool


# Remove other conditions
table(ibs.df$ibd) # 2445 "I do not have this condition"
table(ibs.df$cdiff) # 2593 éI do not have this condition" remove C.difficile
table(ibs.df$fungal_overgrowth) # 1907 "I do not have this condition"
# remove diabetes??
# remove lung_disease, liver_disease, kidney_disease??





