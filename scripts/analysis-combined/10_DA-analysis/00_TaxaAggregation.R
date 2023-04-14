# *********************************
# Purpose: Taxa aggregation in each dataset
# Date: August 2021
# Author: Salomé Carcy
# *********************************




# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(phyloseq)
library(tidyverse)

## 1.2. Data ####
path.root <- "~/Projects/MetaIBS" # CHANGE THIS ROOT DIRECTORY ON YOUR COMPUTER

path.phylobj    <- file.path(path.root, "data/phyloseq-objects")
datasets        <- list.files(path.phylobj, pattern=".rds", include.dirs=FALSE)
phyloseqobjects <- sapply(datasets, function(x) readRDS(file.path(path.phylobj, x)), USE.NAMES=T, simplify=F)
# names(phyloseqobjects) # sanity check

# Change name of phyloseq objects to make it easier later on
names(phyloseqobjects) # sanity check
names(phyloseqobjects) <- c("AGP-2021", "Fukui-2020", "Hugerth-2019", "Labus-2017", "Liu-2020", "LoPresti-2019", "Mars-2020", "Nagel-2016", "Pozuelo-2015", "Ringel-2015", "Zeber-2016", "Zhu-2019", "Zhuang-2018")
names(phyloseqobjects) # sanity check




# ************************
# 2. TAXA AGGREGATION ####
# ************************

#_________________________
# Check that level of aggregation changes the total count at higher taxonomic levels
test.genus <- phyloseqobjects[["Fukui-2020"]] %>%
  tax_glom(taxrank = "Genus") %>% # Genus aggregation
  psmelt() %>%
  group_by(Phylum) %>%
  summarise(genus_agg=sum(Abundance))

test.phylum <- phyloseqobjects[["Fukui-2020"]] %>%
  tax_glom(taxrank = "Phylum") %>% # Phylum aggregation
  psmelt() %>%
  group_by(Phylum) %>%
  summarise(phylum_agg=sum(Abundance))

# Compare the total count per phylum after genus/phylum aggregation
test.genus %>%
  right_join(test.phylum, by="Phylum")
# Different total count (some phyla have lower total count with the genus aggregation)



#_________________________
# Aggregation for all datasets

taxa.levels <- c("Phylum", "Class", "Order", "Family", "Genus")

# Init
i=0

# Iterate through each dataset
for(physeq.cohort in phyloseqobjects){
  
  i=i+1
  cat("++", gsub( "-.*$", "", names(phyloseqobjects[i])), "++\n")
  
  # Aggregate at each taxonomic level
  for (taxa in taxa.levels){
    # Aggregate
    df <- physeq.cohort %>%
      tax_glom(taxrank = taxa) %>%
      psmelt() %>%
      select(-c(Run, sequencing_tech, variable_region))
    cat(taxa, ":", length(unique(df[,taxa])), "\n") # print nb of each taxonomic level
    # Save as csv file
    cohort <- tolower(gsub( "-.*$", "", names(phyloseqobjects[i])))
    file.name <- file.path(path.root, "data/aggregated-tables/",
                           names(phyloseqobjects[i]), # sub-directory
                           paste(paste(cohort, paste(tolower(taxa), "agg", sep="-"), sep="_"), # i.e. fukui_phylum-agg
                                 ".csv", sep=""))
    write.csv(df, file.name)
  }

}
