##########################
# Purpose: select IBS and healthy samples from American Gut Project fecal data
# Date: April 2021
# Author: Salomé Carcy
##########################



####################
# IMPORT LIBRARIES #
####################

library(tidyverse)

path <- "~/Projects/MetaIBS"



###################
# IMPORT METADATA #
###################

# This csv file was downloaded from the SRA
# with the accession number PRJEB11419
sradf <- read.csv(file.path(path, "scripts/analysis-individual/AGP/00_Metadata-AGP/SraRunTable.csv"))
# head(sradf)

# This tsv file was downloaded from the figshare
# link shared by the paper from McDonald et al., 2018
paperdf <- read_tsv(file.path(path, "scripts/analysis-individual/AGP/00_Metadata-AGP/correctedt2.tsv"))
# head(paperdf)

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
  dim
# 8,379 people labeled as 'healthy'


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


# Subset SRADF
# - ppl between 18-69 yo
# - BMI between 16-35
# - no antibiotic in past 6 months
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



###########################
# SELECT IBS & HC SAMPLES #
###########################

# Split healthy & IBS samples
healthyDF <- subsetdf %>%
  filter(subset_healthy == TRUE & 
           ibs == "I do not have this condition") %>%
  mutate(host_disease = "Healthy")
# 4722 "healthy"

ibsDF <- subsetdf %>%
  filter(ibs == "Diagnosed by a medical professional (doctor, physician assistant)") %>%
  mutate(host_disease = "IBS")
# 645 IBS

# We have too many healthy samples, so we will randomly select ~645 healthy samples
healthyDF <- healthyDF %>%
  sample_n(size=nrow(ibsDF))

# Verify there is an equal distribution in age & bmi
t.test(healthyDF$host_age, ibsDF$host_age)
t.test(healthyDF$host_bmi, ibsDF$host_bmi)




# Merge back IBS & healthy samples within same df
# Do re-formatting of some columns for clarification/simplicity
metadata <- bind_rows(healthyDF, ibsDF) %>%
  # Probiotic frequency
  mutate(probiotic_frequency=replace(probiotic_frequency, probiotic_frequency %in% c("Not provided", "Unspecified"), "Unknown")) %>%
  # Gluten consumption (NB: celiac disease was excluded)
  mutate(gluten=replace(gluten, gluten == "No", "I eat gluten"),
         gluten=replace(gluten, gluten == "I do not eat gluten because it makes me feel bad", "I don't eat gluten"),
         gluten=replace(gluten, gluten %in% c("Not provided", "Unspecified"), "Unknown")) %>%
  # Bowel movement
  mutate(bowel_movement_frequency=replace(bowel_movement_frequency, bowel_movement_frequency %in% c("","Not provided","Unspecified"), "Unknown")) %>%
  mutate(bowel_movement_quality=replace(bowel_movement_quality,
                                        bowel_movement_quality %in% c("",grep("know", unique(metadata$bowel_movement_quality), value=TRUE), "Unspecified", "Not provided"),
                                        "Unknown"),
         bowel_movement_quality=replace(bowel_movement_quality,
                                        bowel_movement_quality %in% grep("constipated", unique(metadata$bowel_movement_quality), value=TRUE),
                                        "Constipated"),
         bowel_movement_quality=replace(bowel_movement_quality,
                                        bowel_movement_quality %in% grep("diarrhea", unique(metadata$bowel_movement_quality), value=TRUE),
                                        "Diarrhea"),
         bowel_movement_quality=replace(bowel_movement_quality,
                                        bowel_movement_quality %in% grep("normal", unique(metadata$bowel_movement_quality), value=TRUE),
                                        "Normal")) %>%
  # Small intestinal bacterial overgrowth (sibo)
  mutate(sibo=replace(sibo, sibo %in% c("", "Diagnosed by an alternative medicine practitioner", "Not provided", "Self-diagnosed", "Unspecified"), "Unknown"),
         sibo=replace(sibo, sibo=="Diagnosed by a medical professional (doctor, physician assistant)", "SIBO"),
         sibo=replace(sibo, sibo=="I do not have this condition", "Normal")) %>%
  # Lactose
  mutate(lactose=replace(lactose, lactose %in% c("", "Not provided", "Unspecified"), "Unknown"),
         lactose=replace(lactose, lactose %in% c("FALSE", "No"), "Normal"),
         lactose=replace(lactose, lactose %in% c("TRUE", "Yes"), "Lactose intolerant")) %>%
  # Lung, liver, kidney diseases
  mutate(lung_disease=replace(lung_disease, lung_disease=="Diagnosed by a medical professional (doctor, physician assistant)", "Lung disease"),
         lung_disease=replace(lung_disease, lung_disease=="I do not have this condition", "Normal"),
         lung_disease=replace(lung_disease, !lung_disease %in% c("Lung disease", "Normal"), "Unknown")) %>%
  mutate(liver_disease=replace(liver_disease, liver_disease=="Diagnosed by a medical professional (doctor, physician assistant)", "Liver disease"),
         liver_disease=replace(liver_disease, liver_disease=="I do not have this condition", "Normal"),
         liver_disease=replace(liver_disease, !liver_disease %in% c("Liver disease", "Normal"), "Unknown")) %>%
  mutate(kidney_disease=replace(kidney_disease, kidney_disease=="Diagnosed by a medical professional (doctor, physician assistant)", "Kidney disease"),
         kidney_disease=replace(kidney_disease, kidney_disease=="I do not have this condition", "Normal"),
         kidney_disease=replace(kidney_disease, !kidney_disease %in% c("Kidney disease", "Normal"), "Unknown")) %>%
  mutate(clinical_condition=replace(clinical_condition, clinical_condition=="Diagnosed by a medical professional (doctor, physician assistant)", "Clinical condition"),
         clinical_condition=replace(clinical_condition, clinical_condition=="I do not have this condition", "Normal"),
         clinical_condition=replace(clinical_condition, !clinical_condition %in% c("Clinical condition", "Normal"), "Unknown")) %>%
  
  # Remove useless columns
  select(-c(ibs, subset_healthy)) %>%
  relocate(host_disease, .before=host_age)
  


########
# SAVE #
########

# Export the list of Runs to download
write.table(metadata$Run, "./scripts/analysis-individual/AGP/download-samples/list_files.txt", sep="\t",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

# Export the metadata table
write.csv(metadata, "./scripts/analysis-individual/AGP/00_Metadata-AGP/Metadata-AGP.csv")
