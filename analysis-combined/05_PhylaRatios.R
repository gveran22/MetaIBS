##########################
# Purpose: Plotting Firmicutes/Bacteroidota ratio
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

# Data
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
# physeq.ringel <- readRDS(file.path(path.phy, "physeq_ringel.rds"))
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



###################
# PREPROCESS DATA #
###################

# Merge phyloseq objects
physeq <- merge_phyloseq(physeq.labus,
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

# Separate fecal & sigmoid samples
physeq.fecal <- subset_samples(physeq, sample_type == 'stool') # 2,170 samples
physeq.sigmoid <- subset_samples(physeq, sample_type == 'sigmoid') # 431 samples
cat("Nb of fecal samples:", nsamples(physeq.fecal))
cat("\nNb of sigmoid samples:", nsamples(physeq.sigmoid))


##################################
# PLOT PHYLA RELATIVE ABUNDANCES #
##################################

# Obtain relative abundances
phylum.table <- physeq %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()

# Plot
ggplot(phylum.table, aes(x = reorder(Sample, Sample, function(x) mean(phylum.table[Sample == x & Phylum == 'Bacteroidota', 'Abundance'])),
                         y = Abundance, fill = Phylum))+
  facet_wrap(~ host_disease, scales = "free_x") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c(colorRampPalette(brewer.pal(8, "Set2"))(48)),
                    guide=guide_legend(nrow=24))+
  scale_y_continuous(expand = c(0, 0))+ # remove empty space between axis and plot
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15),
        #legend.position = "None",
        legend.text = element_text(size=7),
        legend.key.size = unit(0.2, 'cm'))+
  labs(x = "Samples", y = "Relative abundance")

# Save figure
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/05_Relative-Abund/phyla_relabund.jpg", width=12, height=5)




