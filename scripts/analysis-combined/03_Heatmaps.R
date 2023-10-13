# ************************************
# Purpose: Plotting annotated heatmaps
# Date: August 2021
# Author: Salomé Carcy
# ************************************




# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(pals)


## 1.2. Data ####
path.root <- "~/Projects/MetaIBS" # CHANGE THIS ROOT DIRECTORY ON YOUR COMPUTER
path.data <- file.path(path.root, "data/analysis-combined/03_Heatmaps")

path.phylobj    <- file.path(path.root, "data/phyloseq-objects/phyloseq-without-phylotree")
datasets        <- list.files(path.phylobj, pattern=".rds")
phyloseqobjects <- sapply(datasets, function(x) readRDS(file.path(path.phylobj, x)), USE.NAMES=T, simplify=F)
# names(phyloseqobjects) # sanity check




# ***********************
# 2. PREPROCESS DATA ####
# ***********************

## 2.1. Merge phyloseq objects ####
physeq.all <- merge_phyloseq(phyloseqobjects[[1]], phyloseqobjects[[2]]) # Merge first two phyloseq objects in the list
# if there are more than 2 phyloseq objects, merge the rest of them
if(length(phyloseqobjects)>2){
  for (i in 3:length(phyloseqobjects)){
    print(paste0("merging with phyloseq object #", i))
    physeq.all <- merge_phyloseq(physeq.all, phyloseqobjects[[i]])
  }
}


## 2.2. Separate fecal & sigmoid samples ####
physeq.fecal <- subset_samples(physeq.all, sample_type == 'stool') # 2,228 samples
physeq.fecal <- prune_taxa(taxa_sums(physeq.fecal)>0, physeq.fecal) # remove ASVs that are not present anymore
cat("Nb of fecal samples:", nsamples(physeq.fecal))

physeq.sigmoid <- subset_samples(physeq.all, sample_type == 'sigmoid') # 431 samples
physeq.sigmoid <- prune_taxa(taxa_sums(physeq.sigmoid)>0, physeq.sigmoid) # remove ASVs that are not present anymore
cat("Nb of sigmoid samples:", nsamples(physeq.sigmoid))


## 2.3. Covariates for heatmap labels ####
color.df <- data.frame(disease = sample_data(physeq.fecal)[,'host_disease'],
                       seq_tech = sample_data(physeq.fecal)[,'sequencing_tech'],
                       author = sample_data(physeq.fecal)[,'author'])
color.df[is.na(color.df$host_disease),"host_disease"] <- "NA"
author.order <- c('Labus', 'LoPresti', 'Ringel', # 454 pyrosequencing
                  'AGP', 'Liu', 'Pozuelo', # Illumina single end
                  'Fukui', 'Hugerth', 'Zhu', 'Zhuang', # Illumina paired end
                  'Nagel', 'Zeber-Lubecka') # Ion Torrent
color.df <- color.df %>%
  mutate(author = factor(color.df$author, levels = author.order)) %>%
  arrange(author, host_disease)

# Order of samples
sample.order <- rownames(color.df)
# table(color.df$author) # sanity check

# Colors for heatmap
annotationCol <- list(host_disease = c(Healthy = '#08519c', IBS = '#ef3b2c', "NA"='black'),
                      sequencing_tech = c('454 pyrosequencing' = '#6a51a3',
                                          'Illumina single-end' = '#a1d99b',
                                          'Illumina paired-end' = '#238b45',
                                          'Ion Torrent' = '#f16913'),
                      author = setNames(brewer.paired(n=12), levels(color.df$author)))




# *****************************
# 3. HEATMAP PHYLA AS ROWS ####
# *****************************

# # Agglomerate to Phylum level, keeping only ASVs present in at least 2 samples
# phylum.agg <- physeq.fecal %>%
#   tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
#   transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
#   psmelt()                                             # Melt to long format
# 
# 
# # Identify phyla present in at least 3 datasets
# list.phylum <- phylum.agg %>%
#   # is each phylum present in each dataset (T/F)?
#   group_by(Phylum, author) %>%
#   summarize(phy_present=sum(Abundance)>0) %>%
#   ungroup() %>%
#   # in how many datasets is each phylum present (n)?
#   group_by(Phylum) %>%
#   count(phy_present) %>%
#   filter(phy_present == TRUE) %>%
#   filter(n>2) %>%
#   ungroup()
# list.phylum <- list.phylum$Phylum
#   
# 
# # Agglomerate again at Phylum level, but keeping only phyla present in at least 3 datasets
# phylum.agg <- subset_taxa(physeq.fecal, Phylum %in% list.phylum) %>%
#   tax_glom(taxrank = "Phylum") %>%
#   transform_sample_counts(function(x) {x/sum(x)} ) %>%
#   psmelt()
# 
# # Get dataframe phylum x samples
# phylumTable <- acast(phylum.agg %>% filter(Phylum %in% list.phylum),
#                      Phylum ~ Sample, fun.aggregate=sum, value.var = 'Abundance')
# 
# # Sanity checks
# dim(phylumTable)
# table(is.na(phylumTable))
# table(colSums(phylumTable))
# table(rownames(phylumTable))
# table(rowSums(phylumTable) == 0)
# 
# # For coloring, add "pseudocounts"
# # min(phylumTable[phylumTable>0]) # min is 2e-5
# phylumTable[phylumTable == 0] <- 10e-6
# 
# 
# # Reorder samples
# phylumTable <- phylumTable[,sample.order] # reorder samples
# 
# jpeg(file.path(path.data, "phylum_heatmp.jpg"), height = 3000, width = 4000, res = 300)
# pheatmap(log10(phylumTable),
#          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
#          show_rownames = T,
#          show_colnames = F,
#          fontsize_row = 8,
#          cluster_rows = T,
#          cluster_cols = T,
#          #cutree_rows = 2,
#          cutree_cols = 2,
#          clustering_method = 'ward.D',
#          annotation_col = color.df,
#          annotation_colors = annotationCol,
#          main = "Hierarchical clustering with ward linkage (phyla as rows)")
# dev.off()




# ********************************
# 4. HEATMAP FAMILIES AS ROWS ####
# ********************************

##___________________________
## 4.1. Agglomerate to Family level ####
family.agg <- physeq.fecal %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()
# saveRDS(family.agg, file.path(path.data, "famGlom_fecal.rds"))  # recommend to save it as takes a while to compute
# family.agg <- readRDS(file.path(path.data, "famGlom_fecal.rds"))


# Identify families present in at least 3 datasets
list.family <- family.agg %>%
  filter(Abundance > 0) %>%
  group_by(Family) %>%
  summarise(nb_datasets = n_distinct(author)) %>%
  filter(nb_datasets>2) %>%
  ungroup()
list.family <- list.family$Family # 120 Families

# Get dataframe family x samples
familyTable <- acast(family.agg %>% filter(Family %in% list.family),
                     Family ~ Sample, fun.aggregate=sum, value.var = 'Abundance')

# Sanity checks
familyTable[1:5,1:5]
dim(familyTable)
table(is.na(familyTable))
table(colSums(familyTable)) # sum per sample
table(rownames(familyTable))
table(rowSums(familyTable) == 0)

# For coloring, add "pseudocounts"
min(familyTable[familyTable>0]) # min is 1e-5
familyTable[familyTable == 0] <- 1e-7


##___________________________
## 4.2. Plot heatmap ####
# Reorder samples
familyTable <- familyTable[,sample.order] # reorder samples

jpeg(file.path(path.data, "fecal_family_heatmp_ordered.jpg"), height = 4000, width = 4000, res = 300)
pheatmap(log10(familyTable),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50),
         show_rownames = T,
         show_colnames = F,
         fontsize_row = 6,
         cluster_rows = T,
         cluster_cols = F,
         # cutree_rows = 2,
         # cutree_cols = 2,
         clustering_method = 'ward.D2',
         annotation_col = color.df,
         annotation_colors = annotationCol,
         main = "Hierarchical clustering with ward linkage (families as rows)")
dev.off()




# **************************************************
# 5. HEATMAP FAMILIES AS ROWS (SIGMOID SAMPLES) ####
# **************************************************

##___________________________
## 5.1. Agglomerate to Family level ####
family.agg.sigmoid <- physeq.sigmoid %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()


# Get dataframe family x samples
familyTable.sigmoid <- acast(family.agg.sigmoid, Family ~ Sample, fun.aggregate=sum, value.var = 'Abundance')

# Sanity checks
familyTable.sigmoid[1:5,1:5]
dim(familyTable.sigmoid)
table(is.na(familyTable.sigmoid))
table(colSums(familyTable.sigmoid)) # sum per sample
table(rownames(familyTable.sigmoid))
table(rowSums(familyTable.sigmoid) == 0)

# For coloring, add "pseudocounts"
# min(familyTable.sigmoid[familyTable.sigmoid>0]) # min is 2e-5
familyTable.sigmoid[familyTable.sigmoid == 0] <- 1e-7


# Covariates for heatmap
color.sigmoid <- data.frame(disease = sample_data(physeq.sigmoid)[,'host_disease'],
                            seq_tech = sample_data(physeq.sigmoid)[,'sequencing_tech'],
                            author = sample_data(physeq.sigmoid)[,'author'])
color.sigmoid <- color.sigmoid %>%
  mutate(author = factor(color.sigmoid$author, levels = c('LoPresti','Hugerth','Mars'))) %>%
  arrange(author, host_disease)
sample.order.sigmoid <- rownames(color.sigmoid) # order of samples
# table(color.sigmoid$author) # sanity check

# Colors for heatmap
annotationCol.sigmoid <- list(host_disease = c(Healthy='#08519c', IBS='#ef3b2c'),
                              sequencing_tech = c('454 pyrosequencing'='#6a51a3', 'Illumina paired-end'='#238b45'),
                              author = setNames(c("#1F78B4", "#FF7F00", "#dfc27d"), levels(color.sigmoid$author)))


##___________________________
## 5.2. Plot heatmap ####

# Reorder samples
familyTable.sigmoid <- familyTable.sigmoid[,sample.order.sigmoid] # reorder samples

# Plot
jpeg(file.path(path.data, "sigmoid_family_heatmp_ordered.jpg"), height = 4000, width = 4000, res = 400)
pheatmap(log10(familyTable.sigmoid),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50),
         show_rownames = T,
         show_colnames = F,
         fontsize_row = 6,
         cluster_rows = T,
         cluster_cols = F,
         # cutree_rows = 2,
         # cutree_cols = 2,
         clustering_method = 'ward.D2',
         annotation_col = color.sigmoid,
         annotation_colors = annotationCol.sigmoid,
         main = "Hierarchical clustering with ward linkage (families as rows)")
dev.off()
