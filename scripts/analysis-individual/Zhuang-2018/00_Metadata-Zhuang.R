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
path <- "~/Projects/MetaIBS/scripts/analysis-individual/Zhuang-2018"
sraDF <- read.table(file.path(path, "00_Metadata-Zhuang/ZhuangSraRunTable.txt"), header = TRUE, sep = ",")


#____________________________________________________________________
# DATAFRAME CLEANUP

# Take a look at table
colnames(sraDF)

# Keep relevant data
sampledf <- sraDF %>%
  as_tibble() %>%
  select(Run, Library.Name) %>%
  mutate(host_subtype=sub("_.*", "", Library.Name)) %>%
  rename(host_ID=Library.Name) %>%
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

# Export list of files to download
write.table(sampledf$Run, file.path(path, "download-Zhuang-samples/list_files_zhuang.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
