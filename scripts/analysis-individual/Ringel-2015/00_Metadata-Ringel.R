##########################
# Purpose: Build metadata table for Ringel-Kulka dataset
# Date: March 2023
# Author: Salomé Carcy
##########################


#____________________________________________________________________
# IMPORT

# Libraries
library(tidyverse)

# Data
path <- "~/Projects/MetaIBS/scripts/analysis-individual/Ringel-2015"
sraDF <- read.table(file.path(path, "00_Metadata-Ringel/RingelSraRunTable.txt"), header = TRUE, sep = ",")


#____________________________________________________________________
# DATAFRAME CLEANUP

# Take a look at table
colnames(sraDF)
head(sraDF) # there isn't any covariates

# Keep relevant covariates
sampledf <- data.frame(sraDF[,"Run"])
rownames(sampledf) <- sraDF$Run
colnames(sampledf) <- "Run"

sampledf$host_disease <- "NA"
sampledf$host_subtype <- "NA"
sampledf$sample_type <- "stool"
sampledf$Collection <- "1st"
sampledf$author <- "Ringel"
sampledf$sequencing_tech <- "454 pyrosequencing"
sampledf$variable_region <- "V1-V2"


#____________________________________________________________________
# SAVE DATAFRAME 

write.csv(sampledf, file.path(path, "00_Metadata-Ringel/Ringel-Metadata.csv"))

# Export list of files to download
write.table(sampledf$Run, file.path(path, "download-Ringel-samples/list_files_ringel.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
