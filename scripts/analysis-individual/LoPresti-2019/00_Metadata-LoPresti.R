##########################
# Purpose: Build metadata table for LoPresti dataset
# Date: March 2023
# Author: Salomé Carcy
##########################



#____________________________________________________________________
# IMPORT

# Libraries
library(tidyverse)

# Data
path.root <- "~/Projects/MetaIBS" # CHANGE THIS PATH ON YOUR COMPUTER
# path.scripts <- file.path(path.root, "scripts/analysis-individual/LoPresti-2019")
path.data    <- file.path(path.root, "data/analysis-individual/LoPresti-2019")

# Metadata table downloaded from the SRA
sra_metadata <- read.csv(file.path(path.data, "00_Metadata-LoPresti/LoPrestiSraRunTable.csv"))
head(sra_metadata)




#____________________________________________________________________
# DATAFRAME CLEANUP

clean.df <- sra_metadata %>%
  # Select relevant columns and rename some
  select(all_of(c("Run", "host_disease", "env_material", "Host_Age", "host_sex",
                  "gastrointest_disord", "Sample_Label", "Collection_Date"))) %>%
  dplyr::rename(sample_type=env_material,
                host_age=Host_Age,
                host_subtype=gastrointest_disord,
                host_ID=Sample_Label,
                collection_date=Collection_Date) %>%
  # Rename some variables
  mutate(host_disease = replace(host_disease, host_disease=="Healthy subject", "Healthy"),
         host_disease = replace(host_disease, host_disease=="Irritable Bowel Syndrome", "IBS"),
         sample_type  = replace(sample_type,  sample_type=="intestinal mucosa", "sigmoid"),
         host_subtype = replace(host_subtype, host_subtype=="Absence", "HC"),
         host_subtype = replace(host_subtype, host_subtype=="Alternating constipation and diarrhea", "IBS-M"),
         host_subtype = replace(host_subtype, host_subtype=="diarrhea", "IBS-D"),
         host_subtype = replace(host_subtype, host_subtype=="constipation", "IBS-C"),
         host_ID      = gsub("\\.A$", '', host_ID),
         host_ID      = gsub("\\A$", '', host_ID)) %>%
  # Add some variables
  mutate(Collection = "1st",
         author     = "LoPresti",
         sequencing_tech = "454 pyrosequencing",
         variable_region = "V1-V3")

# Sanity check that paired samples have same demographics
clean.df %>%
  group_by(host_ID) %>%
  filter(n_distinct(Run)>1) %>%
  filter(n_distinct(host_age)>1) %>%
  arrange(host_ID)
# a few samples (10 individuals => 20 samples) have different age, but only 1 year difference, so likely that
# the individuals just had their birthday in-between collection of stool & sigmoid

# Finish
rownames(clean.df) <- clean.df$Run
head(clean.df)




#____________________________________________________________________
# SAVE DATAFRAME

# Export metadata table
write.csv(clean.df, file.path(path.data, "00_Metadata-LoPresti/Metadata-LoPresti.csv"))

# Export the list of Runs to download
write.table(clean.df$Run, file.path(path.data, "raw_fastq/list_files_lopresti.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
