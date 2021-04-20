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





