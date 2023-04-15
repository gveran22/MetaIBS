##########################
# Purpose: Build metadata table for Fukui dataset
# Date: March 2023
# Author: Salomé Carcy
##########################



#____________________________________________________________________
# IMPORT

# Libraries
library(tidyverse)

# Data
path.root <- "~/Projects/MetaIBS" # CHANGE THIS PATH ON YOUR COMPUTER
path.scripts <- file.path(path.root, "scripts/analysis-individual/Fukui-2020")
path.data    <- file.path(path.root, "data/analysis-individual/Fukui-2020")
sraDF <- read.table(file.path(path.data, "00_Metadata-Fukui/FukuiSraRunTable.txt"), header = TRUE, sep = ",")



#____________________________________________________________________
# DATAFRAME CLEANUP

# Keep relevant columns
metadata_table <- sraDF[, c("Run", "Sample.Name", "Collection_date")]

metadata_table[grep("IBS*",metadata_table$Sample.Name), "host_disease"] <- "IBS"
metadata_table[grep("HC*",metadata_table$Sample.Name), "host_disease"] <- "Healthy"

metadata_table[metadata_table$host_disease == "Healthy", "host_subtype"] <- "HC"
metadata_table[metadata_table$host_disease == "IBS", "host_subtype"] <- "IBS-unspecified"

# Sanity check
table(metadata_table$host_disease) # 85 IBS & 26 Healthy
table(metadata_table$host_subtype) # 85 IBS-unspecified & 26 HC

metadata_table$Sample.Name <- NULL
colnames(metadata_table)[2] <- "collection_date"
rownames(metadata_table) <- metadata_table$Run
head(metadata_table)

metadata_table$sample_type <- "stool"
metadata_table$author <- "Fukui"
metadata_table$sequencing_tech <- "Illumina paired-end"
metadata_table$variable_region <- "V1-V2"
metadata_table$Collection <- "1st" # only 1 collection time point per individual
head(metadata_table)



#____________________________________________________________________
# SAVE DATAFRAME 
write.csv(metadata_table, file.path(path.data, "00_Metadata-Fukui/Metadata-Fukui.csv"))

# Export list of files to download
write.table(sampledf$Run, file.path(path.scripts, "download-Fukui-samples/list_files_fukui.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)