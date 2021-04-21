##########################
# Purpose: Taxa aggregation in each dataset
# Date: April 2021
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(tidyverse)

# Data
path <- "~/Projects/IBS_Meta-analysis_16S"
physeq.fukui <- readRDS(file.path(path, "phyloseq-objects/physeq_fukui.rds"))
physeq.hugerth <- readRDS(file.path(path, "phyloseq-objects/physeq_hugerth.rds"))
physeq.labus <- readRDS(file.path(path, "phyloseq-objects/physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path, "phyloseq-objects/physeq_lopresti.rds"))
physeq.nagel <- readRDS(file.path(path, "phyloseq-objects/physeq_nagel.rds"))
physeq.pozuelo <- readRDS(file.path(path, "phyloseq-objects/physeq_pozuelo.rds"))
physeq.zeber <- readRDS(file.path(path, "phyloseq-objects/physeq_zeber.rds"))
physeq.zhu <- readRDS(file.path(path, "phyloseq-objects/physeq_zhu.rds"))
physeq.zhuang <- readRDS(file.path(path, "phyloseq-objects/physeq_zhuang.rds"))



####################
# TAXA AGGREGATION #
####################

#_________________________
# Check that level of aggregation changes the total count at higher taxonomic levels
test.genus <- physeq.fukui %>%
  tax_glom(taxrank = "Genus") %>% # Genus aggregation
  psmelt() %>%
  group_by(Phylum) %>%
  summarise(Abundance=sum(Abundance)) %>%
  rename(genus_agg = Abundance)

test.phylum <- physeq.fukui %>%
  tax_glom(taxrank = "Phylum") %>% # Phylum aggregation
  psmelt() %>%
  group_by(Phylum) %>%
  summarise(Abundance=sum(Abundance)) %>%
  rename(phylum_agg = Abundance)

# Compare the total count per phylum after genus/phylum aggregation
test.genus %>%
  right_join(test.phylum, by="Phylum")
# Different total count (some phyla have lower total count with the genus aggregation)













