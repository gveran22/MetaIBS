# ********************************************************
# Purpose: merge phyloseq object and keep only common ASVs
# Date: October 2021
# Author: Salomé Carcy
# ********************************************************




# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(tidyverse)
library(phyloseq)
library(reshape2)
library(ggupset)
library(cowplot)

## 1.2. Data ####
path.root   <- "~/Projects/MetaIBS" # CHANGE THIS ROOT DIRECTORY ON YOUR CLUSTER
path.data  <- file.path(path.root, "data/analysis-combined/05_Common-ASVs") 

path.phylobj    <- file.path(path.root, "data/phyloseq-objects/phyloseq-without-phylotree")
datasets        <- list.files(path.phylobj, pattern=".rds")
phyloseqobjects <- sapply(datasets, function(x) readRDS(file.path(path.phylobj, x)), USE.NAMES=T, simplify=F)

# Change name of phyloseq objects to make it easier later on
names(phyloseqobjects) # sanity check
names(phyloseqobjects) <- gsub("physeq_", "", names(phyloseqobjects))
names(phyloseqobjects) <- gsub(".rds", "", names(phyloseqobjects))
names(phyloseqobjects)# sanity check


# If the phyloseq objects contain phylogenetic trees, they need to be removed
# for (i in 1:length(phyloseqobjects)){
#   print(names(phyloseqobjects)[i])
#   # Get the phyloseq object
#   physeq <- phyloseqobjects[[i]]
#   # Replace rownames in tax_table and colnames in otu_table by ASV sequences (and not "ASV1, ASV2, ...")
#   taxa_names(physeq) <- as.character(refseq(physeq))
#   # Remove phylogenetic tree
#   physeq <- phyloseq(otu_table(physeq),
#                      sample_data(physeq), 
#                      tax_table(physeq),
#                      refseq(physeq))
#   # Replace the phyloseq object in the list
#   phyloseqobjects[[i]] <- physeq
# }
# # Sanity checks
# print(phyloseqobjects[[1]]) # shouldn't have a phy_tree() slot
# print(taxa_names(phyloseqobjects[[1]])[1:3]) # should be DNA sequences (the ASV sequences) and not "ASV1, ASV2, ASV3"




# **************************************************
# 2. MERGE PHYLOSEQ OBJECTS WITH MERGE_PHYLOSEQ ####
# **************************************************

# Merge phyloseq objects
physeq.all <- merge_phyloseq(phyloseqobjects[[1]], phyloseqobjects[[2]]) # Merge first two phyloseq objects in the list
# if there are more than 2 phyloseq objects, merge the rest of them
if(length(phyloseqobjects)>2){
  for (i in 3:length(phyloseqobjects)){
    print(paste0("merging with phyloseq object #", i))
    physeq.all <- merge_phyloseq(physeq.all, phyloseqobjects[[i]])
  }
}

# Compare the number of ASVs before/after merging
sum_taxa <- 0
for(i in 1:length(phyloseqobjects)){
  sum_taxa <- sum_taxa+ntaxa(phyloseqobjects[[i]])
}
print(sum_taxa) # 81,474 before merging

ntaxa(physeq.all) # 79,943 after merging


# Let's get all ASVs in a dataframe and check if we can find common ones
asv.df <- melt(lapply(phyloseqobjects, function(x) taxa_names(x)))
colnames(asv.df) <- c("asv", "author")
length(unique(asv.df$asv)) # we do find 79,943 unique sequences

# Let's see which datasets share the exact same ASV sequence
common.asv <- asv.df %>%
  group_by(asv) %>%
  summarize(Datasets = list(unique(author))) %>%
  filter(lengths(Datasets)>1)

# jpeg(file.path(path.data, "commonASV_merge-phyloseq-funct.jpg"), width=2000, height=2200, res=400)
ggplot(common.asv, aes(x=Datasets))+
  geom_bar(fill="#737373") +
  scale_x_upset()+
  theme_cowplot()+
  labs(y="# common ASVs")
# dev.off()




# *****************************************************
# 3. CREATE PHYLOSEQ OBJECTS WITH ONLY COMMON ASVs ####
# *****************************************************

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
# ntaxa(physeq.common.NagPoz) # 806
# ntaxa(physeq.common.LiuZhug) # 475
# ntaxa(physeq.common.HugZhu) # 231
# ntaxa(physeq.common.LopRin) # 19

# psmelt(physeq.common.NagPoz) %>%
#   group_by(OTU) %>%
#   summarize(Datasets = list(unique(author))) %>%
#   filter(length(Datasets)>1) # should give 806 OTUs that are present in both datasets

# See how many counts the common ASVs represent
# sum(otu_table(physeq.common.NagPoz)) / ( sum(otu_table(phyloseqobjects$nagel)) + sum(otu_table(phyloseqobjects$pozuelo)) )
# sum(otu_table(subset_samples(physeq.common.NagPoz, author == "Nagel"))) / sum(otu_table(phyloseqobjects$nagel))
# sum(otu_table(subset_samples(physeq.common.NagPoz, author == "Pozuelo"))) / sum(otu_table(phyloseqobjects$pozuelo))

# sum(otu_table(physeq.common.LiuZhug)) / ( sum(otu_table(phyloseqobjects$liu)) + sum(otu_table(phyloseqobjects$zhuang)) )
# sum(otu_table(subset_samples(physeq.common.LiuZhug, author == "Liu"))) / sum(otu_table(phyloseqobjects$liu))
# sum(otu_table(subset_samples(physeq.common.LiuZhug, author == "Zhuang"))) / sum(otu_table(phyloseqobjects$zhuang))

# sum(otu_table(physeq.common.HugZhu)) / ( sum(otu_table(phyloseqobjects$hugerth)) + sum(otu_table(phyloseqobjects$zhu)) )
# sum(otu_table(physeq.common.LopRin)) / ( sum(otu_table(phyloseqobjects$lopresti)) + sum(otu_table(phyloseqobjects$ringel)) )


# Save the phyloseq objects
saveRDS(physeq.common.NagPoz,  file.path(path.data, "physeq_commonASV_Nagel-Pozuelo.rds"))
saveRDS(physeq.common.LiuZhug, file.path(path.data, "physeq_commonASV_Liu-Zhuang.rds"))
saveRDS(physeq.common.HugZhu,  file.path(path.data, "physeq_commonASV_Hugerth-Zhu.rds"))
saveRDS(physeq.common.LopRin,  file.path(path.data, "physeq_commonASV_Lopresti-Ringel.rds"))




# *********************************
# 4. CREATE LONG SHAPED TABLES ####
# *********************************

# If need to re-import data
# physeq.common.NagPoz  <- readRDS(file.path(path.data, "physeq_commonASV_Nagel-Pozuelo.rds"))
# physeq.common.LiuZhug <- readRDS(file.path(path.data, "physeq_commonASV_Liu-Zhuang.rds"))
# physeq.common.HugZhu  <- readRDS(file.path(path.data, "physeq_commonASV_Hugerth-Zhu.rds"))
# physeq.common.LopRin  <- readRDS(file.path(path.data, "physeq_commonASV_Lopresti-Ringel.rds"))

# Columns we can remove (only AGP dataset has info on these columns, and "Run" is being duplicated in a "Sample" column)
col_to_remove <- c("Run", "age_cat", "bmi_cat", "country", "RACE", "exercise_frequency", "alcohol_frequency", "probiotic_frequency", "gluten",
                   "bowel_movement_frequency", "bowel_movement_quality", "collection_season", "sibo", "lactose", "lung_disease", "liver_disease",
                   "kidney_disease", "clinical_condition")


# Get long-shaped dataframes
df.common.NagPoz  <- psmelt(physeq.common.NagPoz) %>% select(-all_of(col_to_remove))
df.common.LiuZhug <- psmelt(physeq.common.LiuZhug) %>% select(-all_of(col_to_remove))
df.common.HugZhu  <- psmelt(physeq.common.HugZhu) %>% select(-all_of(col_to_remove))
df.common.LopRin  <- psmelt(physeq.common.LopRin) %>% select(-all_of(col_to_remove))


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
write.csv(df.common.NagPoz,  file.path(path.data, "commonASV_Nagel-Pozuelo.csv"))
write.csv(df.common.LiuZhug, file.path(path.data, "commonASV_Liu-Zhuang.csv"))
write.csv(df.common.HugZhu,  file.path(path.data, "commonASV_Hugerth-Zhu.csv"))
write.csv(df.common.LopRin,  file.path(path.data, "commonASV_Lopresti-Ringel.csv"))

