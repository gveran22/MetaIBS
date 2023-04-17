##########################
# Purpose: Build metadata table for Liu dataset
# Date: March 2023
# Author: Salomé Carcy
##########################



#____________________________________________________________________
# IMPORT

# Libraries
library(data.table)
library(tidyverse)

# Data
path.root <- "~/Projects/MetaIBS" # CHANGE THIS PATH ON YOUR COMPUTER
path.scripts <- file.path(path.root, "scripts/analysis-individual/Liu-2020")
path.data    <- file.path(path.root, "data/analysis-individual/Liu-2020")

# Metadata table downloaded from the SRA
sra_metadata <- read.csv(file.path(path.data, "00_Metadata-Liu/LiuSraRunTable.txt"), header=TRUE, sep = ",")
head(sra_metadata)




#____________________________________________________________________
# DATAFRAME CLEANUP

clean.df <- sra_metadata %>%
  # Select relevant columns
  select(all_of(c("Run","Isolation_source", "sex", "AGE", "BMI", "Bristol", "IBS_SSS", "Collection_date"))) %>%
  # Reformat DF
  dplyr::rename(host_subtype = Isolation_source,
         host_sex = sex,
         host_age = AGE,
         host_bmi = BMI,
         collection_date = Collection_date) %>%
  mutate(host_disease = replace(host_subtype, host_subtype=="HC", "Healthy")) %>%
  mutate(host_disease = str_replace(host_disease, "IBS-D", "IBS")) %>%
  mutate(sample_type = "stool") %>%
  relocate(host_disease, .after=Run) %>%
  relocate(host_subtype, .after=IBS_SSS) %>%
  # Add some additional covariates useful for future combined analyses
  mutate(author="Liu", sequencing_tech="Illumina single-end", variable_region="V3-V4", Collection="1st")
rownames(clean.df) <- clean.df$Run
head(clean.df)




#____________________________________________________________________
# SAVE DATAFRAME 

# Export metadata table
write.csv(clean.df, file.path(path.data, "00_Metadata-Liu/Metadata-Liu.csv"))

# Export the list of Runs to download
write.table(clean.df$Run, file.path(path.scripts, "download-Liu-samples/list_files_liu.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)