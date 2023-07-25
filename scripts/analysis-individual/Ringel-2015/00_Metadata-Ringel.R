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
path <- "~/Projects/MetaIBS/data/analysis-individual/Ringel-2015"
sraDF <- read.table(file.path(path, "00_Metadata-Ringel/RingelSraRunTable.txt"), header = TRUE, sep = ",")


#____________________________________________________________________
# DATAFRAME CLEANUP

# Take a look at table
colnames(sraDF)
head(sraDF) # there isn't any covariates
dim(sraDF) # there's double the number of samples we expect (2*76=152)
table(sraDF$Bases) # every nb of bases appears twice, further proof that samples are duplicated
table(sraDF[,c("LibrarySelection", "BioSample")]) # clearly showing samples are duplicated

# Keep only half of the samples
sraDF <- sraDF[sraDF$BioSample=="SAMN04272978",]

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

write.csv(sampledf, file.path(path, "00_Metadata-Ringel/Metadata-Ringel.csv"))

# Export the list of Runs to download
# Open the filereport downloaded from the ENA
ena_table <- read.csv(file.path(path, "raw_fastq/filereport_read_run_PRJNA302437_tsv.txt"), header=T, sep="\t")
table(ena_table$run_accession %in% sampledf$Run, useNA="ifany") # 76 samples of interest and 76 that we don't want to download
ena_table <- ena_table[ena_table$run_accession %in% sampledf$Run,]
dim(ena_table) # 76 rows
# Export the list of links to download the Runs from ENA
write.table(ena_table$fastq_ftp, file.path(path, "raw_fastq/list_files_ringel.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)





# Export list of files to download
write.table(sampledf$Run, file.path(path, "download-Ringel-samples/list_files_ringel.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
