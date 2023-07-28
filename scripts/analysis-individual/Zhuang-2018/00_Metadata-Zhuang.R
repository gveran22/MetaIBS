##########################
# Purpose: Build metadata table for Zhuang dataset
# Date: March 2023
# Author: Salomé Carcy
##########################


#____________________________________________________________________
# IMPORT

# Libraries
library(tidyverse)

# Data
path <- "~/Projects/MetaIBS/data/analysis-individual/Zhuang-2018"
sraDF <- read.table(file.path(path, "00_Metadata-Zhuang/ZhuangSraRunTable.txt"), header = TRUE, sep = ",")


#____________________________________________________________________
# DATAFRAME CLEANUP

# Take a look at table
colnames(sraDF)

# Keep relevant data
sampledf <- sraDF %>%
  as_tibble() %>%
  dplyr::select(Run, Library.Name) %>%
  mutate(host_subtype=sub("_.*", "", Library.Name)) %>%
  dplyr::rename(host_ID=Library.Name) %>%
  # Keep only healthy/IBS samples
  filter(host_subtype %in% c("HC", "IBS")) %>%
  # Rename host_subtype values
  mutate(host_subtype=ifelse(host_subtype=="IBS", "IBS-D", "HC"), # only IBS-diarrhea patients in this cohort
         host_disease=ifelse(host_subtype=="IBS-D", "IBS", "Healthy"),
         Collection="1st",
         sample_type="stool",
         author="Zhuang",
         sequencing_tech="Illumina paired-end",
         variable_region="V3-V4") %>%
  as.data.frame()
rownames(sampledf) <- sampledf$Run


#____________________________________________________________________
# SAVE DATAFRAME 
write.csv(sampledf, file.path(path, "00_Metadata-Zhuang/Metadata-Zhuang.csv"))

# Open the filereport downloaded from the ENA
ena_table <- read.csv(file.path(path, "raw_fastq/filereport_read_run_PRJNA475187_tsv.txt"), header=T, sep="\t")
table(ena_table$run_accession %in% sampledf$Run, useNA="ifany") # 30 samples of interest
ena_table <- ena_table[ena_table$run_accession %in% sampledf$Run,]
dim(ena_table) # 30 rows, like sraDF
head(ena_table)
# the fastq_ftp column contains both forward & reverse reads separated by ";"
ena_tableF <- ena_table[,c("run_accession", "fastq_ftp")]
ena_tableF$fastq_ftp <- gsub(";.*", "", ena_tableF$fastq_ftp)
ena_tableR <- ena_table[,c("run_accession", "fastq_ftp")]
ena_tableR$fastq_ftp <- gsub(".*;", "", ena_tableR$fastq_ftp)
ena_table_final <- rbind(ena_tableF, ena_tableR)
head(ena_table_final[order(ena_table_final$run_accession),])
table(ena_table_final$run_accession, useNA="ifany") # should all appear twice

# Export the list of links to download the Runs from ENA
write.table(ena_table_final$fastq_ftp, file.path(path, "raw_fastq/list_files_zhuang.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
