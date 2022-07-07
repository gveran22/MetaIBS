##########################
# Purpose: Compute Bray-Curtis distance between all samples & plot by PCoA
# Date: July 2022
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(ggplot2)
library(cowplot)
library(tidyverse)


# Data
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
physeq.ringel   <- readRDS(file.path(path.phy, "physeq_ringel.rds"))
physeq.labus    <- readRDS(file.path(path.phy, "physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path.phy, "physeq_lopresti.rds"))
physeq.pozuelo  <- readRDS(file.path(path.phy, "physeq_pozuelo.rds"))
physeq.zhuang   <- readRDS(file.path(path.phy, "physeq_zhuang.rds"))
physeq.zhu      <- readRDS(file.path(path.phy, "physeq_zhu.rds"))
physeq.hugerth  <- readRDS(file.path(path.phy, "physeq_hugerth.rds"))
physeq.fukui    <- readRDS(file.path(path.phy, "physeq_fukui.rds"))
physeq.mars     <- readRDS(file.path(path.phy, "physeq_mars.rds"))
physeq.liu      <- readRDS(file.path(path.phy, "physeq_liu.rds"))
physeq.agp      <- readRDS(file.path(path.phy, "physeq_agp.rds"))
physeq.nagel    <- readRDS(file.path(path.phy, "physeq_nagel.rds"))
physeq.zeber    <- readRDS(file.path(path.phy, "physeq_zeber.rds"))




###################
# PREPROCESS DATA #
###################

# Merge phyloseq objects
physeq.all <- merge_phyloseq(physeq.ringel,
                             physeq.labus,
                             physeq.lopresti,
                             physeq.pozuelo,
                             physeq.zhuang,
                             physeq.zhu,
                             physeq.hugerth,
                             physeq.fukui,
                             physeq.mars,
                             physeq.liu,
                             physeq.agp,
                             physeq.nagel,
                             physeq.zeber)

# Sanity check
physeq.all # 2,651 samples and 79,943 taxa


# Separate fecal & sigmoid samples
physeq.fecal <- subset_samples(physeq.all, sample_type == 'stool') # 2,220 samples
physeq.fecal <- prune_taxa(taxa_sums(physeq.fecal)>0, physeq.fecal) # remove ASVs that are not present anymore
cat("Nb of fecal samples:", nsamples(physeq.fecal))

physeq.sigmoid <- subset_samples(physeq.all, sample_type == 'sigmoid') # 431 samples
physeq.sigmoid <- prune_taxa(taxa_sums(physeq.sigmoid)>0, physeq.sigmoid) # remove ASVs that are not present anymore
cat("\nNb of sigmoid samples:", nsamples(physeq.sigmoid))




#################
# STOOL SAMPLES #
#################

# *********************
# BRAY-CURTIS DISTANCE 
# *********************

# Agglomerate to Family level
physeq_fecal.family <- tax_glom(physeq.fecal, "Family") # 5 min
# View(tax_table(physeq_fecal.family)) #  sanity check

# Common-scale normalization
physeq_fecal.CSN <- transform_sample_counts(physeq_fecal.family, function(x) (x*min(sample_sums(physeq_fecal.family))) / sum(x) )
# table(rowSums(otu_table(physeq_fecal.CSN))) # check that all rows are summing to the same total
# saveRDS(physeq_fecal.CSN, "~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/02_DimensionReduction/physeq_fecal_CSN.rds")

# Compute Bray-Curtis
glom_fecal.bray <- phyloseq::distance(physeq_fecal.CSN, method = "bray") # real long to compute

# Compute principal coordinate analysis (PCoA)
set.seed(123)
iMDS_fecal.Bray <- ordinate(physeq=physeq_fecal.CSN, method="PCoA", distance=glom_fecal.bray)
plot_fecal.df <- plot_ordination(physeq_fecal.CSN, iMDS_fecal.Bray)$data



# **********
# PLOT PCoA
# **********

# Per author/dataset
ggplot(plot_fecal.df, aes(Axis.1, Axis.2, color=author))+
  geom_point(size=2, alpha=0.5) +
  scale_color_manual(values=pals::brewer.paired(12), name="")+ # author
  labs(title = "Dataset")+
  theme_cowplot()+
  theme(line = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Per host_disease
ggplot(plot_fecal.df, aes(Axis.1, Axis.2, color=host_disease))+
  geom_point(size=.5) +
  scale_color_manual(values = c('blue', 'red', 'black'), name="")+ # disease
  labs(title = "Disease phenotype")+
  theme_cowplot()+
  theme(line = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Per seq. technology
ggplot(plot_fecal.df, aes(Axis.1, Axis.2, color=sequencing_tech))+
  geom_point(size=.5) +
  scale_color_manual(values=c("#6a51a3", "#a1d99b", "#238b45", "#f16913"), name="")+ # seqtech
  labs(title = "Seq. technology")+
  theme_cowplot()+
  theme(line = element_blank(),
        plot.title = element_text(hjust = 0.5))




###################
# SIGMOID SAMPLES #
###################
