##########################
# Purpose: select IBS and healthy samples from American Gut Project fecal data
# Date: April 2021
# Author: Salomé Carcy
##########################



####################
# IMPORT LIBRARIES #
####################

library(tidyverse)

path <- "~/Projects/IBS_Meta-analysis_16S"



###################
# IMPORT METADATA #
###################

# This csv file was downloaded from the SRA
# with the accession number PRJEB11419
sradf <- read.csv(file.path(path, "data/analysis-individual/AGP/SraRunTable.csv"))

# This tsv file was downloaded from the figshare
# link shared by the paper from McDonald et al., 2018
paperdf <- read_tsv(file.path(path, "data/analysis-individual/AGP/correctedt2.tsv"))

# Compare two tables
dim(sradf) # 33,610 samples
dim(paperdf) # 17,854 samples

# The paper was published in 2018, so it is possible that more samples have been added on the SRA since
table(str_sub(paperdf$collection_date, start=-4)) # latest sample from paper metadata table is 2017
table(str_sub(sradf$ENA.FIRST.PUBLIC..run., end=4)) # latest sample from SRA metadata is 2021

# => We will extract sample Runs from SRA metadata


###################
# CLEANUP METADATA #
###################

# For simplicity purposes, we will keep only samples and columns of interest

# Keep only gut metadata
sradf %>%
  filter(target_gene..exp. == "16S rRNA") %>% # 16s rRNA seq
  filter(SAMPLE_TYPE %in% c("feces", "Stool", "stool") & Organism %in% c("feces metagenome", "human gut metagenome", "gut metagenome")) # gut metagenome

# Check that the "subset_healthy" corresponds to:
# - ppl between 20-69 yo
# - BMI between 18.5-30
# - no antibiotic in past year
# - no IBD
# - no diabetes
sradf %>%
  filter(subset_healthy == TRUE) %>%
  # filter(subset_age == TRUE) %>%
  # filter(bmi_cat %in% c("Normal", "Overweight")) %>%
  # filter(diabetes == "I do not have this condition") %>%
  # filter(ibd == "I do not have this condition") %>%
  # filter(antibiotic_history == "I have not taken antibiotics in the past year.") %>%
  dim
# 11,261 people labeled as 'healthy'


# Remove columns that will not be of use
skip_columns <- c("SAMPLE_TYPE", "Organism",
                  "BioProject",
                  "center_name..exp.",
                  "Center.Name",
                  "Consent",
                  "ENA.FIRST.PUBLIC..run.",
                  "LibrarySelection",
                  "LibrarySource",
                  "Platform",
                  "ReleaseDate",
                  "sequencing_meth..exp.",
                  "SRA.Study",
                  "well_description..exp.",
                  "run_prefix..exp.",
                  "run_date..exp.",
                  "Bytes",
                  "DATASTORE.provider",
                  "DATASTORE.region",
                  "Title",
                  "sample_plate..exp.",
                  "target_gene..exp.",
                  "target_subfragment..exp.",
                  "exercise_location",
                  "cosmetics_frequency",
                  "pool_frequency",
                  "last_travel",
                  "Elevation",
                  "contraceptive",
                  "LATITUDE",
                  "LONGITUDE",
                  "fermented_plant_frequency",
                  "dominant_hand",
                  "livingwith",
                  "flu_vaccine_date",
                  "types_of_plants",
                  "consume_animal_products_abx",
                  "last_move",
                  "age_corrected",
                  "primer_plate..exp.",
                  "well_id..exp.",
                  "roommates_in_study",
                  "collection_time")
sradf <- sradf %>%
  select(-all_of(skip_columns))




# Cleanup columns: bowel_movement_frequency, bowel_movement_quality


table(sradf$acne_medication)
 

###############
# IBS SAMPLES #
###############

# Subset table to individuals with IBS medical diagnosis
table(sradf$ibs) # 2895 "Diagnosed by a medical professional (doctor, physician assistant)"
ibs.df <- sradf[grep("doctor", sradf$ibs),]
dim(ibs.df) # 2895 patients
table(ibs.df$bowel_movement_quality) 

# Check it is only fecal samples
table(ibs.df$Organism) # 2671 human gut metagenome
table(ibs.df$SAMPLE_TYPE) # 2383 stool


# Remove other conditions
table(ibs.df$ibd) # 2445 "I do not have this condition"
table(ibs.df$cdiff) # 2593 éI do not have this condition" remove C.difficile
table(ibs.df$fungal_overgrowth) # 1907 "I do not have this condition"
# remove diabetes??
# remove lung_disease, liver_disease, kidney_disease??





