##########################
# Purpose: Build metadata table for Zhu dataset
# Date: March 2023
# Author: Salomé Carcy
##########################


#____________________________________________________________________
# IMPORT

# Libraries
library(tidyverse)

# Data
path <- "~/Projects/MetaIBS/data/analysis-individual/Zhu-2019" # CHANGE THIS PATH ON YOUR COMPUTER
sraDF <- read.table(file.path(path, "00_Metadata-Zhu/ZhuSraRunTable.txt"), header = TRUE, sep = ",")


#____________________________________________________________________
# DATAFRAME CLEANUP

# Take a look at table
colnames(sraDF)
head(sraDF)

# Keep relevant data
sampledf <- sraDF[, c("Run", "Sample.Name", "Age", "sex")]
rownames(sampledf) <- sraDF$Run
colnames(sampledf) <- c("Run", "host_disease", "host_age", "host_sex")

# Simplify data
sampledf$host_disease <- as.character(sampledf$host_disease)
sampledf$host_disease[substr(sampledf$host_disease, start = 1, stop = 2) == "N_"] <- "Healthy" # Replace with better names
sampledf$host_disease[substr(sampledf$host_disease, start = 1, stop = 2) == "F_"] <- "IBS" # Replace with better names

sampledf[sampledf$host_disease == "Healthy", "host_subtype"] <- "HC"
sampledf[sampledf$host_disease == "IBS", "host_subtype"] <- "IBS-unspecified"

sampledf$Collection <- "1st" # only 1 collection point per individual
sampledf$sample_type <- "stool"
sampledf$author <- "Zhu"
sampledf$sequencing_tech <- "Illumina paired-end"
sampledf$variable_region <- "V4"


#____________________________________________________________________
# SAVE DATAFRAME 
write.csv(sampledf, file.path(path, "00_Metadata-Zhu/Metadata-Zhu.csv"))

# Open the filereport downloaded from the ENA
ena_table <- read.csv(file.path(path, "raw_fastq/filereport_read_run_PRJNA566284_tsv.txt"), header=T, sep="\t")
table(ena_table$run_accession %in% sampledf$Run, useNA="ifany") # all samples are present
dim(ena_table) # 29 rows
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
write.table(ena_table_final$fastq_ftp, file.path(path, "raw_fastq/list_files_zhu.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
