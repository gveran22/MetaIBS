##########################
# Purpose: Build metadata table for Zeber-Lubecka dataset
# Date: March 2023
# Author: Salomé Carcy
##########################


#____________________________________________________________________
# IMPORT

# Libraries
library(tidyverse)

# Data
path <- "~/Projects/MetaIBS/data/analysis-individual/Zeber-2016"
sraDF <- read.table(file.path(path, "00_Metadata-Zeber/ZeberSraRunTable.txt"), header = TRUE, sep = ",")
additional_metadata <- read.table(file.path(path, "00_Metadata-Zeber/metadata_sex.txt"), header=FALSE)
colnames(additional_metadata) <- c("title", "host_sex")


#____________________________________________________________________
# DATAFRAME CLEANUP

# Take a look at table
colnames(sraDF)
head(sraDF)

# Keep relevant covariates
table(sraDF$Alias %in% additional_metadata$title)
table(sraDF$title %in% additional_metadata$title) # this column better corresponds to the additional_metadata identifiers

# Clean up dataframe
sampledf <- sraDF %>%
  as_tibble() %>%
  select(Run, title, Description) %>%
  separate(col=Description, into=c("host_ID", "treatment"), sep=":") %>%
  mutate(treatment = replace_na(treatment, "0")) %>%
  # keep only patients who didn't get RIF treatment
  filter(treatment != "RIF") %>%
  select(!treatment) %>%
  # define host_disease and host_subtype based on host_ID column
  mutate(host_ID=str_remove(host_ID, " "),
         host_disease=ifelse(str_detect(host_ID, "PATIENT"), "IBS", "Healthy"),
         host_subtype=ifelse(str_detect(host_ID, "M")==T, "IBS-M",
                      ifelse(str_detect(host_ID, "D")==T, "IBS-D",
                      ifelse(str_detect(host_ID, "[0-9]C")==T, "IBS-C",
                      ifelse(str_detect(host_ID, "CONTROL")==T, "HC", "Problem"))))) %>%
  # add host_sex covariate
  left_join(additional_metadata, by="title") %>%
  mutate(host_sex=ifelse(str_detect(host_sex, "F")==T, "female",
                  ifelse(str_detect(host_sex, "M")==T, "male", "unknown"))) %>%
  # remove unwanted columns & add useful columns
  select(!title) %>%
  mutate(sample_type="stool",
         author="Zeber-Lubecka",
         sequencing_tech="Ion Torrent",
         variable_region="Ion 16S Metagenomics Kit",
         Collection="1st")

# Add rownames
sampledf <- as.data.frame(sampledf)
rownames(sampledf) <- sampledf$Run

# sanity check
table(sampledf$host_disease, useNA="ifany")


#____________________________________________________________________
# SAVE DATAFRAME 
write.csv(sampledf, file.path(path, "00_Metadata-Zeber/Metadata-Zeber.csv"))

# Export the list of Runs to download
# Open the filereport downloaded from the ENA
ena_table <- read.csv(file.path(path, "raw_fastq/filereport_read_run_PRJEB11252_tsv.txt"), header=T, sep="\t")
table(ena_table$run_accession %in% sampledf$Run, useNA="ifany") # 90 samples of interest and 85 that we don't want to download
ena_table <- ena_table[ena_table$run_accession %in% sampledf$Run,]
dim(ena_table) # 90 rows
# Export the list of links to download the Runs from ENA
write.table(ena_table$fastq_ftp, file.path(path, "raw_fastq/list_files_zeber.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
