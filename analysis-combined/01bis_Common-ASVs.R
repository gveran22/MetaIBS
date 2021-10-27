##########################
# Purpose: merge phyloseq object and keep only common ASVs
# Date: October 2021
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(tidyverse)
library(phyloseq)
library(reshape2)
library(ggupset)

# Data
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
physeq.ringel <- readRDS(file.path(path.phy, "physeq_ringel.rds"))
physeq.labus <- readRDS(file.path(path.phy, "physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path.phy, "physeq_lopresti.rds"))
physeq.pozuelo <- readRDS(file.path(path.phy, "physeq_pozuelo.rds"))
physeq.zhuang <- readRDS(file.path(path.phy, "physeq_zhuang.rds"))
physeq.zhu <- readRDS(file.path(path.phy, "physeq_zhu.rds"))
physeq.hugerth <- readRDS(file.path(path.phy, "physeq_hugerth.rds"))
physeq.fukui <- readRDS(file.path(path.phy, "physeq_fukui.rds"))
physeq.mars <- readRDS(file.path(path.phy, "physeq_mars.rds"))
physeq.liu <- readRDS(file.path(path.phy, "physeq_liu.rds"))
physeq.agp <- readRDS(file.path(path.phy, "physeq_agp.rds"))
physeq.nagel <- readRDS(file.path(path.phy, "physeq_nagel.rds"))
physeq.zeber <- readRDS(file.path(path.phy, "physeq_zeber.rds"))




##############################################
# MERGE PHYLOSEQ OBJECTS WITH MERGE_PHYLOSEQ #
##############################################

# Merge phyloseq objects
physeq.all <- merge_phyloseq(physeq.labus,
                             physeq.lopresti,
                             physeq.ringel,
                             physeq.agp,
                             physeq.liu,
                             physeq.pozuelo,
                             physeq.fukui,
                             physeq.mars,
                             physeq.hugerth,
                             physeq.zhu,
                             physeq.zhuang,
                             physeq.nagel,
                             physeq.zeber)

# Compare the number of ASVs before/after merging
sum(ntaxa(physeq.labus)+
      ntaxa(physeq.lopresti)+
      ntaxa(physeq.ringel)+
      ntaxa(physeq.agp)+
      ntaxa(physeq.liu)+
      ntaxa(physeq.pozuelo)+
      ntaxa(physeq.fukui)+
      ntaxa(physeq.mars)+
      ntaxa(physeq.hugerth)+
      ntaxa(physeq.zhu)+
      ntaxa(physeq.zhuang)+
      ntaxa(physeq.nagel)+
      ntaxa(physeq.zeber)) # 81,452 before

ntaxa(physeq.all) # 79,918 after


# Put datasets in a list
datasets <- list("labus"    = physeq.labus,
                 "lopresti" = physeq.lopresti,
                 "ringel"   = physeq.ringel,
                 "agp"      = physeq.agp,
                 "liu"      = physeq.liu,
                 "pozuelo"  = physeq.pozuelo,
                 "fukui"    = physeq.fukui,
                 "mars"     = physeq.mars,
                 "hugerth"  = physeq.hugerth,
                 "zhu"      = physeq.zhu,
                 "zhuang"   = physeq.zhuang,
                 "nagel"    = physeq.nagel,
                 "zeber"    = physeq.zeber)


# Let's get all ASVs in a dataframe and check if we can find common ones
asv.df <- melt(lapply(datasets, function(x) taxa_names(x)))
colnames(asv.df) <- c("asv", "author")
length(unique(asv.df$asv)) # we do find 79,918 unique sequences

# Let's see which datasets share the exact same ASV sequence
common.asv <- asv.df %>%
  group_by(asv) %>%
  summarize(Datasets = list(unique(author))) %>%
  filter(lengths(Datasets)>1)

# jpeg("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/01_Merge-Datasets/commonASV_merge-phyloseq-funct.jpg", width=2000, height=2000, res=400)
ggplot(common.asv, aes(x=Datasets))+
  geom_bar() +
  scale_x_upset()+
  labs(y="# of common ASVs")
# dev.off()




#################################################
# CREATE PHYLOSEQ OBJECTS WITH ONLY COMMON ASVs #
#################################################

# Get a phyloseq object with only the identified common ASVs
physeq.common <- prune_taxa(taxa_names(physeq.all) %in% common.asv$asv, physeq.all)


# Remove samples with 0 count
physeq.common <- prune_samples(sample_sums(physeq.common)>0, physeq.common)


# Create a phyloseq object for each datasets having common ASVs
physeq.common.NagPoz <- subset_samples(physeq.common, author %in% c("Nagel", "Pozuelo"))
physeq.common.NagPoz <- prune_taxa(taxa_sums(physeq.common.NagPoz)>0, physeq.common.NagPoz)

physeq.common.LiuZhug <- subset_samples(physeq.common, author %in% c("Liu", "Zhuang"))
physeq.common.LiuZhug <- prune_taxa(taxa_sums(physeq.common.LiuZhug)>0, physeq.common.LiuZhug)

physeq.common.HugZhu <- subset_samples(physeq.common, author %in% c("Hugerth", "Zhu"))
physeq.common.HugZhu <- prune_taxa(taxa_sums(physeq.common.HugZhu)>0, physeq.common.HugZhu)

physeq.common.LopRin <- subset_samples(physeq.common, author %in% c("LoPresti", "Ringel"))
physeq.common.LopRin <- prune_taxa(taxa_sums(physeq.common.LopRin)>0, physeq.common.LopRin)


# Sanity check
# ntaxa(physeq.common.NagPoz)
# ntaxa(physeq.common.LiuZhug)
# ntaxa(physeq.common.HugZhu)
# ntaxa(physeq.common.LopRin)

# psmelt(physeq.common.NagPoz) %>%
#   group_by(OTU) %>%
#   summarize(Datasets = list(unique(author))) %>%
#   filter(length(Datasets)>1) # should give 806 OTUs that are present in both datasets


# Save the phyloseq objects
path <- "~/Projects/IBS_Meta-analysis_16S/phyloseq-objects/common-ASVs"
saveRDS(physeq.common.NagPoz, file.path(path, "physeq_commonASV_Nagel-Pozuelo.rds"))
saveRDS(physeq.common.LiuZhug, file.path(path, "physeq_commonASV_Liu-Zhuang.rds"))
saveRDS(physeq.common.HugZhu, file.path(path, "physeq_commonASV_Hugerth-Zhu.rds"))
saveRDS(physeq.common.LopRin, file.path(path, "physeq_commonASV_Lopresti-Ringel.rds"))




#############################
# CREATE LONG SHAPED TABLES #
#############################

# If need to re-import data
# physeq.common.NagPoz <- readRDS(file.path(path, "physeq_commonASV_Nagel-Pozuelo.rds"))
# physeq.common.LiuZhug <- readRDS(file.path(path, "physeq_commonASV_Liu-Zhuang.rds"))
# physeq.common.HugZhu <- readRDS(file.path(path, "physeq_commonASV_Hugerth-Zhu.rds"))
# physeq.common.LopRin <- readRDS(file.path(path, "physeq_commonASV_Lopresti-Ringel.rds"))

# Columns we can remove (only AGP dataset has info on these columns, and "Run" is being duplicated in a "Sample" column)
col_to_remove <- c("Run", "age_cat", "bmi_cat", "country", "RACE", "exercise_frequency", "alcohol_frequency", "probiotic_frequency", "gluten",
                   "bowel_movement_frequency", "bowel_movement_quality", "collection_season", "sibo", "lactose", "lung_disease", "liver_disease",
                   "kidney_disease", "clinical_condition")


# Get long-shaped dataframes
df.common.NagPoz <- psmelt(physeq.common.NagPoz) %>% select(-all_of(col_to_remove))
df.common.LiuZhug <- psmelt(physeq.common.LiuZhug) %>% select(-all_of(col_to_remove))
df.common.HugZhu <- psmelt(physeq.common.HugZhu) %>% select(-all_of(col_to_remove))
df.common.LopRin <- psmelt(physeq.common.LopRin) %>% select(-all_of(col_to_remove))


# Sanity checks
# length(unique(df.common.NagPoz$OTU))
# length(unique(df.common.LiuZhug$OTU))
# length(unique(df.common.HugZhu$OTU))
# length(unique(df.common.LopRin$OTU))
# 
# sum(df.common.NagPoz$Abundance) == sum(otu_table(physeq.common.NagPoz))
# sum(df.common.LiuZhug$Abundance) == sum(otu_table(physeq.common.LiuZhug))
# sum(df.common.HugZhu$Abundance) == sum(otu_table(physeq.common.HugZhu))
# sum(df.common.LopRin$Abundance) == sum(otu_table(physeq.common.LopRin))
# 
# unique(df.common.NagPoz$author)
# unique(df.common.LiuZhug$author)
# unique(df.common.HugZhu$author)
# unique(df.common.LopRin$author)


# Save tables!
path <- "~/Projects/IBS_Meta-analysis_16S/aggregated-tables/common-ASVs"
write.csv(df.common.NagPoz, file.path(path, "commonASV_Nagel-Pozuelo.csv"))
write.csv(df.common.LiuZhug, file.path(path, "commonASV_Liu-Zhuang.csv"))
write.csv(df.common.HugZhu, file.path(path, "commonASV_Hugerth-Zhu.csv"))
write.csv(df.common.LopRin, file.path(path, "commonASV_Lopresti-Ringel.csv"))

