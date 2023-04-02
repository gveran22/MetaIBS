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
path <- "~/Projects/MetaIBS/scripts/analysis-individual/Zhu-2019"
sraDF <- read.table(file.path(path, "00_Metadata-Zhu/ZhuSraRunTable.txt"), header = TRUE, sep = ",")


#____________________________________________________________________
# DATAFRAME CLEANUP

# Take a look at table
colnames(sraDF)
head(sraDF)

# Keep relevant data
sampledf <- sraDF[, c("Run", "Sample.Name", "Age", "sex")]
rownames(sampledf) <- paste(sraDF$Run) # have the same sample.names
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

# Export list of files to download
write.table(sampledf$Run, file.path(path, "download-Zhu-samples/list_files_zhu.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
