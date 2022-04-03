##########################
# Purpose: Plotting annotated heatmaps
# Date: August 2021
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(pals)

# Data
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
physeq.labus    <- readRDS(file.path(path.phy, "physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path.phy, "physeq_lopresti.rds"))
physeq.ringel   <- readRDS(file.path(path.phy, "physeq_ringel.rds"))
physeq.agp      <- readRDS(file.path(path.phy, "physeq_agp.rds"))
physeq.liu      <- readRDS(file.path(path.phy, "physeq_liu.rds"))
physeq.pozuelo  <- readRDS(file.path(path.phy, "physeq_pozuelo.rds"))
physeq.fukui    <- readRDS(file.path(path.phy, "physeq_fukui.rds"))
physeq.hugerth  <- readRDS(file.path(path.phy, "physeq_hugerth.rds"))
physeq.mars     <- readRDS(file.path(path.phy, "physeq_mars.rds"))
physeq.zhu      <- readRDS(file.path(path.phy, "physeq_zhu.rds"))
physeq.zhuang   <- readRDS(file.path(path.phy, "physeq_zhuang.rds"))
physeq.nagel    <- readRDS(file.path(path.phy, "physeq_nagel.rds"))
physeq.zeber    <- readRDS(file.path(path.phy, "physeq_zeber.rds"))




###################
# PREPROCESS DATA #
###################

# Merge phyloseq objects
physeq <- merge_phyloseq(physeq.labus,
                         physeq.lopresti,
                         physeq.ringel,
                         physeq.agp,
                         physeq.liu,
                         physeq.pozuelo,
                         physeq.fukui,
                         physeq.hugerth,
                         physeq.mars,
                         physeq.zhu,
                         physeq.zhuang,
                         physeq.nagel,
                         physeq.zeber)

# Separate fecal & sigmoid samples
physeq.fecal <- subset_samples(physeq, sample_type == 'stool') # 2,220 samples
physeq.fecal <- prune_taxa(taxa_sums(physeq.fecal)>0, physeq.fecal) # remove ASVs that are not present anymore
cat("Nb of fecal samples:", nsamples(physeq.fecal))

# physeq.sigmoid <- subset_samples(physeq, sample_type == 'sigmoid') # 431 samples
# physeq.sigmoid <- prune_taxa(taxa_sums(physeq.sigmoid)>0, physeq.sigmoid) # remove ASVs that are not present anymore
# cat("\nNb of sigmoid samples:", nsamples(physeq.sigmoid))


# Covariates for heatmap labels
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




#########################
# HEATMAP PHYLA AS ROWS #
#########################

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
# jpeg("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/09_Heatmaps/phylum_heatmp.jpg", height = 3000, width = 4000, res = 300)
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




############################
# HEATMAP FAMILIES AS ROWS #
############################

# Agglomerate to Family level
family.agg <- physeq.fecal %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()
# saveRDS(family.agg, "~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/09_Heatmaps/famGlom_fecal.rds")
# family.agg <- readRDS("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/09_Heatmaps/famGlom_fecal.rds")


# Identify families present in at least 3 datasets
list.family <- family.agg %>%
  filter(Abundance > 0) %>%
  group_by(Family) %>%
  summarise(nb_datasets = n_distinct(author)) %>%
  filter(nb_datasets>2) %>%
  ungroup()
list.family <- list.family$Family # 116 Families


# Agglomerate again at family level, but keeping only families present in at least 3 datasets
# family.agg <- subset_taxa(physeq.fecal, Family %in% list.family) %>%
#   tax_glom(taxrank = "Family") %>%
#   transform_sample_counts(function(x) {x/sum(x)} ) %>%
#   psmelt()

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
# min(familyTable[familyTable>0]) # min is 6e-6
familyTable[familyTable == 0] <- 1e-7

# Reorder samples
familyTable <- familyTable[,sample.order] # reorder samples

jpeg("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/09_Heatmaps/family_heatmp_ordered.jpg", height = 4000, width = 4000, res = 300)
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


# Try heatmap with ggplot?
# ggplot(family.agg %>% filter(Family %in% list.family),
#        aes(Sample, Family, fill=Abundance)) + 
#   geom_tile()



