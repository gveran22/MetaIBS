##########################
# Purpose: Build metadata table for Pozuelo dataset
# Date: August 2021
# Author: Salomé Carcy
##########################


#____________________________________________________________________
# IMPORT

# Libraries
library(tidyverse)

# Data
path <- "~/Projects/MetaIBS/data/analysis-individual/Pozuelo-2015"
sraDF <- read.table(file.path(path, "00_Metadata-Pozuelo/PozueloSraRunTable.txt"), header = TRUE, sep = ",")
subtypeDF <- read.csv(file.path(path, "00_Metadata-Pozuelo/Pozuelo_hostsubtype.csv"))


#____________________________________________________________________
# DATAFRAME CLEANUP

# Cleanup SRA metadata
sampledf <- sraDF %>%
  dplyr::select(c("Run", "Sample.Name", "Collection_Date")) %>%
  mutate(Sample.Name=as.character(Sample.Name),
         Collection=ifelse(gsub('.*(?=.{2}$)', '', Sample.Name, perl=T) == '.2', "2nd", "1st"),
         host_disease=ifelse(substr(Sample.Name, 1, 2) == "HC", "Healthy", "IBS"),
         sample_type="stool",
         author="Pozuelo",
         sequencing_tech="Illumina single-end",
         variable_region="V4") %>%
  # remodel the Sample.Name
  mutate(Sample.Name=gsub("MO", "IBS", Sample.Name),
         Sample.Name=gsub("\\.2$", "", Sample.Name)) %>%
  # rename columns
  dplyr::rename(host_ID=Sample.Name,
                collection_date=Collection_Date)

# Join with subtypeDF
sampledf <- subtypeDF %>%
  dplyr::select(Run, host_subtype) %>%
  left_join(sampledf, by="Run") %>%
  mutate(host_subtype = replace(host_subtype, !host_subtype %in% c("HC", "IBS-D", "IBS-C", "IBS-M"), "TO-DISCARD"),
         host_disease = replace(host_disease, !host_subtype %in% c("HC", "IBS-D", "IBS-C", "IBS-M"), "TO-DISCARD"))
table(sampledf$host_disease, useNA="ifany") # sanity check: we should have 88 HC and 185 IBS

# Keep only healthy/IBS samples
sampledf <- sampledf %>%
  filter(host_disease != "TO-DISCARD")
rownames(sampledf) <- sampledf$Run


#____________________________________________________________________
# SAVE DATAFRAME 

write.csv(sampledf, file.path(path, "00_Metadata-Pozuelo/Metadata-Pozuelo.csv"))

# Export the list of Runs to download
# Open the filereport downloaded from the ENA
ena_table <- read.csv(file.path(path, "raw_fastq/filereport_read_run_PRJNA268708_tsv.txt"), header=T, sep="\t")
table(ena_table$run_accession %in% sampledf$Run, useNA="ifany") # 273 samples of interest and 17 that we don't want to download
ena_table <- ena_table[ena_table$run_accession %in% sampledf$Run,]
dim(ena_table) # 273 rows
# Export the list of links to download the Runs from ENA
write.table(ena_table$fastq_ftp, file.path(path, "raw_fastq/list_files_pozuelo.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
