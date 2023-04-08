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
path <- "~/Projects/MetaIBS/scripts/analysis-individual/Pozuelo-2015"
sraDF <- read.table(file.path(path, "00_Metadata-Pozuelo/PozueloSraRunTable.txt"), header = TRUE, sep = ",")
subtypeDF <- read.csv(file.path(path, "00_Metadata-Pozuelo/host_subtype.csv"))


#____________________________________________________________________
# DATAFRAME CLEANUP

# Keep relevant covariates
sampledf <- sraDF[, c("Run", "Sample.Name", "Collection_Date")]
rownames(sampledf) <- sampledf$Run

sampledf$Sample.Name <- as.character(sampledf$Sample.Name)

# Collection date
sampledf$Collection <- "1st" # add column to specify if sample is from first or 2nd collection
sampledf[gsub('.*(?=.{2}$)', '', sampledf$Sample.Name, perl=T) == '.2', 'Collection'] <- "2nd"

# host_disease
sampledf$host_disease <- "Healthy" # add column for host_disease
sampledf[substr(sampledf$Sample.Name, 1, 2) == "MO", 'host_disease'] <- "IBS"

# Sample.Name (host_ID)
sampledf[, 'Sample.Name'] <- gsub("MO", "IBS", sampledf[, 'Sample.Name']) # replace 'MO' by 'IBS' in the sample.names
sampledf$Sample.Name <- gsub("\\.2$", "", sampledf[, 'Sample.Name']) # remove the '.2' at the end of some of the sample.names (indicating it is the 2nd collection)

colnames(sampledf) <- c("Run", "host_ID", "collection_date", "Collection", "host_disease")

sampledf$sample_type <- "stool"
sampledf$author <- "Pozuelo"
sampledf$sequencing_tech <- "Illumina single-end"
sampledf$variable_region <- "V4"

# Join with subtypeDF
sampledf <- subtypeDF %>%
  select(Run, host_subtype) %>%
  left_join(sampledf, by="Run") %>%
  mutate(host_subtype = replace(host_subtype, !host_subtype %in% c("HC", "IBS-D", "IBS-C", "IBS-M"), "TO-DISCARD"),
         host_disease = replace(host_disease, !host_subtype %in% c("HC", "IBS-D", "IBS-C", "IBS-M"), "TO-DISCARD"))
table(sampledf$host_disease) # sanity check

# Keep only healthy/IBS samples
sampledf <- sampledf %>%
  filter(host_disease != "TO-DISCARD")


#____________________________________________________________________
# SAVE DATAFRAME 

write.csv(sampledf, file.path(path, "00_Metadata-Pozuelo/Pozuelo-Metadata.csv"))

# Export list of files to download
write.table(sampledf$Run, file.path(path, "download-Pozuelo-samples/list_files_pozuelo.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
