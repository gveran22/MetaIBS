##########################
# Purpose: Build metadata table for Jacobs dataset
# Date: September 2023
# Author: Salomé Carcy
##########################



#____________________________________________________________________
# IMPORT

# Libraries
library(tidyverse)

# Data
path.root <- "~/Projects/MetaIBS" # CHANGE THIS PATH ON YOUR COMPUTER
path.scripts <- file.path(path.root, "scripts/analysis-individual/Jacobs-2023")
path.data    <- file.path(path.root, "data/analysis-individual/Jacobs-2023")
sraDF <- read.table(file.path(path.data, "00_Metadata-Jacobs/JacobsSraRunTable.txt"), header = TRUE, sep = ",")
bioDF <- readxl::read_xlsx(file.path(path.data, "00_Metadata-Jacobs/IBS_cross-sectional_cohort_16S_metadata.xlsx"), sheet = 1) # shared by author <3

# 318 IBS and 177 HC joined study
# 312 IBS and 174 HC included for 16S => total of 486 samples, sequenced in 2 batches (292+195)


#____________________________________________________________________
# DATAFRAME CLEANUP

useful_columns <- c("Run", "Assay.Type", "ETHNICITY", "gastrointest_disord", "Host_Age", "host_sex", "host_subject_id", "Library.Name",
                    "host_body_mass_index", "Host_Diet", "host_phenotype", "misc_param", "host_occupation")


sraDF.clean <- sraDF[,useful_columns] %>%
  as_tibble() %>%
  # keep only 16S samples
  # filter(Assay.Type=="AMPLICON") %>%
  # select(-Assay.Type) %>%
  # rename some columns
  dplyr::rename(ethnicity=ETHNICITY,
                host_disease=gastrointest_disord,
                host_age=Host_Age,
                host_ID=host_subject_id,
                sample_name=Library.Name,
                host_bmi=host_body_mass_index,
                host_diet=Host_Diet,
                host_subtype=host_phenotype,
                host_education=host_occupation) %>%
  # replace some variables
  mutate(host_disease=ifelse(host_disease=="Control", "Healthy", "IBS"),
         host_subtype=dplyr::case_when(host_subtype=="Constipation" ~ "IBS-C",
                                       host_subtype=="Diarrhea"     ~ "IBS-D",
                                       host_subtype=="Mixed"        ~ "IBS-M",
                                       host_subtype=="Unspecified"  ~ "IBS-unspecified",
                                       host_subtype=="Normal" & host_disease=="IBS" ~ "IBS-unspecified",
                                       host_disease=="Healthy" ~ "HC")) #%>%
  # keep only "baseline" samples (no "followup" post treatment)
  # filter(misc_param=="Baseline") %>%
  # select(-misc_param)


bioDF.clean <- bioDF[,c("*sample_name", "ethnicity", "gastrointest_disord", "host_age", "host_body_mass_index",
                        "host_diet", "host_occupation", "host_phenotype", "host_sex", "host_subject_id")] %>%
  as_tibble() %>%
  # rename some columns
  dplyr::rename("sample_name"="*sample_name",
                host_disease=gastrointest_disord,
                host_ID=host_subject_id,
                host_bmi=host_body_mass_index,
                host_subtype=host_phenotype,
                host_education=host_occupation) %>%
  # replace some variables
  mutate(host_disease=ifelse(host_disease=="Control", "Healthy", "IBS"),
         host_subtype=dplyr::case_when(host_subtype=="Constipation" ~ "IBS-C",
                                       host_subtype=="Diarrhea"     ~ "IBS-D",
                                       host_subtype=="Mixed"        ~ "IBS-M",
                                       host_subtype=="Unspecified"  ~ "IBS-unspecified",
                                       host_subtype=="Normal" & host_disease=="IBS" ~ "IBS-unspecified",
                                       host_disease=="Healthy" ~ "HC")) %>%
  arrange(sample_name)

# sanity check nb of HC/IBS sample
table(bioDF.clean$host_disease, useNA="ifany") # 175 healthy and 311 IBS
table(bioDF.clean$host_subtype, useNA="ifany")
length(unique(bioDF.clean$host_ID)) # every individual has 1 sample

# Check that all bioDF samples are present in sraDF
table(bioDF.clean$sample_name %in% sraDF.clean$sample_name, useNA="ifany")
bioDF.clean[!bioDF.clean$sample_name %in% sraDF.clean$sample_name,] # which sample is not present? A6414_Cross_sectional
sraDF.clean[sraDF.clean$host_ID=="A6414",] # seems like this sample is called "1UCLA_470" in the sra DF
# so we'll replace that name to A6414_Cross_sectional (as it seems all other metadata matches)

# Subset sraDF to samples of interest
sraDF.clean.sub <- sraDF.clean %>%
  mutate(sample_name=str_replace(sample_name, "1UCLA_470", "A6414_Cross_sectional")) %>%
  filter(sample_name %in% bioDF.clean$sample_name) %>%
  select(Run, sample_name, host_disease, host_age) %>%
  arrange(sample_name)

# sanity checks
# table(sraDF.clean.sub$sample_name==bioDF.clean$sample_name, useNA="ifany")
# table(sraDF.clean.sub$host_disease==bioDF.clean$host_disease, useNA="ifany")
# table(sraDF.clean.sub$host_age==bioDF.clean$host_age, useNA="ifany")
# all metadata matches, perfect!!

# Add Run to bioDF, and we'll use bioDF as the metadata DF
metadata.final <- bioDF.clean %>%
  left_join(sraDF.clean.sub, by=c("sample_name", "host_disease", "host_age")) %>%
  relocate(Run) %>%
  # finish adding important details to metadata dataframe
  select(-sample_name) %>%
  mutate(sample_type="stool",
         author="Jacobs",
         sequencing_tech="Illumina paired-end",
         variable_region="V4",
         Collection="1st") %>%
  as.data.frame()
rownames(metadata.final) <- metadata.final$Run



#____________________________________________________________________
# SAVE DATAFRAME
write.csv(metadata.final, file.path(path.data, "00_Metadata-Jacobs/Metadata-Jacobs.csv"))

# Export the list of Runs to download
# Open the filereport downloaded from the ENA
ena_table <- read.csv(file.path(path.data, "raw_fastq/filereport_read_run_PRJNA812699_tsv.txt"), header=T, sep="\t")
table(ena_table$run_accession %in% metadata.final$Run, useNA="ifany") # 486 samples of interest we want to download (698 we don't care about)
ena_table <- ena_table[ena_table$run_accession %in% metadata.final$Run,]
dim(ena_table) # 486 rows
# head(ena_table)
# the fastq_ftp column contains both forward & reverse reads separated by ";"
ena_tableF <- ena_table[,c("run_accession", "fastq_ftp")]
ena_tableF$fastq_ftp <- gsub(";.*", "", ena_tableF$fastq_ftp)
ena_tableR <- ena_table[,c("run_accession", "fastq_ftp")]
ena_tableR$fastq_ftp <- gsub(".*;", "", ena_tableR$fastq_ftp)
ena_table_final <- rbind(ena_tableF, ena_tableR)
head(ena_table_final[order(ena_table_final$run_accession),])
table(ena_table_final$run_accession, useNA="ifany") # should all appear twice

# Export the list of links to download the Runs from ENA
write.table(ena_table_final$fastq_ftp, file.path(path.data, "raw_fastq/list_files_jacobs.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
