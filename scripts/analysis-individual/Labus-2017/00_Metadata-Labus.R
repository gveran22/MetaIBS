##########################
# Purpose: Build metadata table for Labus dataset
# Date: March 2023
# Author: Salomé Carcy
##########################



#____________________________________________________________________
# IMPORT

# Libraries
library(tidyverse)
library(readxl)

# Data
path.root <- "~/Projects/MetaIBS" # CHANGE THIS PATH ON YOUR COMPUTER
path.scripts <- file.path(path.root, "scripts/analysis-individual/Labus-2017")
path.data    <- file.path(path.root, "data/analysis-individual/Labus-2017")
sraDF   <- read.table(file.path(path.data, "00_Metadata-Labus/LabusSraRunTable.txt"), header = TRUE, sep = ",")
extraDF <- read_xlsx(file.path(path.data, "00_Metadata-Labus/MetaFileLabusMicrobiome201716sMEtaBHrequest.xlsx"), sheet=1)[,1:7]



#____________________________________________________________________
# DATAFRAME CLEANUP

# Keep relevant columns in SRA df
sraDF <- sraDF[, c("Run", "Collection_Date", "host_disease", "Host_Age", "host_sex", "host_subject_ID", "host_body_mass_index")]
colnames(sraDF) <- c("Run", "collection_date", "host_disease", "host_age", "host_sex", "host_ID", "host_bmi")

# Keep relevant columns in extra df shared by the authors
colnames(extraDF) <- c("host_ID", "Cluster", "Group", "host_disease", "BH", "host_subtype", "host_sex")
extraDF <- extraDF[, c("host_ID", "host_subtype", "host_sex")]

# Merge DFs
metadata <- sraDF %>%
  left_join(extraDF, by="host_ID") %>%
  as_tibble() %>%
  dplyr::rename(host_sex=host_sex.x) %>%
  mutate(host_disease=ifelse(host_disease=="None", "Healthy", "IBS"),
         host_subtype=str_replace(host_subtype, "IBS-A", "IBS-M"))

# Sanity checks
# table(metadata[,c("host_disease", "host_subtype")])
# table(metadata[,c("host_sex", "host_sex.y")])
metadata <- metadata %>% select(!host_sex.y)


# Add some additional relevant covariates (for future merging with other datasets)
metadata <- metadata %>%
  mutate(author="Labus", sequencing_tech="454 pyrosequencing", variable_region="V3-V5", sample_type="stool", Collection="1st") %>%
  as.data.frame()
rownames(metadata) <- metadata$Run



#____________________________________________________________________
# SAVE DATAFRAME 
write.csv(metadata, file.path(path.data, "00_Metadata-Labus/Metadata-Labus.csv"))

# No need to export list of files to download, we'll download all .fastq files directly from the ENA
# write.table(metadata$Run, file.path(path.scripts, "download-Labus-samples/list_files_labus.txt"),
#             sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)