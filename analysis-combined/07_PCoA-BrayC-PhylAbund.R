##########################
# Purpose: PCoA of Bray-Curtis dissimilarity with phylum rel. abundance
# Date: August 2021
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(ggplot2)
library(umap)
library(tidyverse)
library(reshape2)
library(gtools)

# Data
path.phy <- "~/IBS/PhyloTree/input"
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




###################
# PREPROCESS DATA #
###################

# Merge phyloseq objects
cat("\n++ MERGE PHYLOSEQ OBJECTS ++\n")
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

# Separate fecal & sigmoid samples
physeq.fecal <- subset_samples(physeq.all, sample_type == 'stool') # 2,228 samples
# physeq.sigmoid <- subset_samples(physeq, sample_type == 'sigmoid') # 431 samples
cat("Nb of fecal samples:", nsamples(physeq.fecal))
# cat("\nNb of sigmoid samples:", nsamples(physeq.sigmoid))




#################################
# GET PHYLA RELATIVE ABUNDANCES #
#################################

cat("\n++ GET RELATIVE ABUNDANCES ++\n")

# Agglomerate to phylum level and get relative abundance
phylumTable.rel <- physeq.fecal %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()

# Extract relative abundance of major phyla
bacter.rel <- phylumTable.rel %>%
  filter(Phylum == "Bacteroidota") %>%
  select(Run, Abundance, Phylum) %>%
  rename(Bacteroidota = Abundance) %>%
  select(-Phylum)

firmi.rel <- phylumTable.rel %>%
  filter(Phylum == "Firmicutes") %>%
  select(Run, Abundance, Phylum) %>%
  rename(Firmicutes = Abundance) %>%
  select(-Phylum)

actino.rel <- phylumTable.rel %>%
  filter(Phylum == "Actinobacteriota") %>%
  select(Run, Abundance, Phylum) %>%
  rename(Actinobacteriota = Abundance) %>%
  select(-Phylum)

proteo.rel <- phylumTable.rel %>%
  filter(Phylum == "Proteobacteria") %>%
  select(Run, Abundance, Phylum) %>%
  rename(Proteobacteria = Abundance) %>%
  select(-Phylum)

# Join tables
relAbundance <- left_join(x=bacter.rel, y=firmi.rel, by="Run") %>%
  left_join(actino.rel, by="Run") %>%
  left_join(proteo.rel, by="Run")

# Add the phyla rel. abundances to the sample_data()
sampledf <- as_tibble(sample_data(physeq.fecal)) %>%
  left_join(relAbundance, by="Run")
sampledf <- as.data.frame(sampledf)
rownames(sampledf) <- sampledf$Run
sample_data(physeq.fecal) <- sampledf




###################
# PLOT ORDINATION #
###################

cat("\n++ RUN PCoA ++\n")

# Perform common-scale normalization
physeq.CSN <- physeq.fecal
physeq.CSN <- transform_sample_counts(physeq.CSN, function(x) (x*min(sample_sums(physeq.CSN))) / sum(x) )

# Calculate distance
ord <- ordinate(physeq.CSN, "PCoA", "bray")
saveRDS(ord, "~/IBS/PCoA/data/ord_bray.rds")


cat("\n++ PLOT PCoA ++\n")
# Plot
plot_ordination(physeq.CSN, ord, type="samples", color="host_disease")
ggsave("~/IBS/PCoA/data/pcoa_bray_disease.jpg", width=6, height=4, type="cairo")

plot_ordination(physeq.CSN, ord, type="samples", color="author")
ggsave("~/IBS/PCoA/data/pcoa_bray_author.jpg", width=6, height=4, type="cairo")

plot_ordination(physeq.CSN, ord, type="samples", color="sequencing_tech")
ggsave("~/IBS/PCoA/data/pcoa_bray_seqtech.jpg", width=6, height=4, type="cairo")


cat("\n++ RUN PCoA WITH RELATIVE ABUNDANCES ++\n")
# Plot phyla rel. abundances
plot_ordination(physeq.CSN, ord, type="samples", color="Firmicutes")+
  geom_point(size=1) +
  scale_color_gradient2(midpoint=mean(sample_data(physeq.fecal)$Firmicutes), low="blue", mid="white", high="red")+
  theme(panel.grid = element_blank(),
        panel.background=element_blank())
ggsave("~/IBS/PCoA/data/pcoa_bray_firmicutes.jpg", width=6, height=4, type="cairo")

plot_ordination(physeq.CSN, ord, type="samples", color="Bacteroidota")+
  geom_point(size=1) +
  scale_color_gradient2(midpoint=mean(sample_data(physeq.fecal)$Bacteroidota), low="blue", mid="white", high="red")+
  theme(panel.grid = element_blank(),
        panel.background=element_blank())
ggsave("~/IBS/PCoA/data/pcoa_bray_bacteroidota.jpg", width=6, height=4, type="cairo")

plot_ordination(physeq.CSN, ord, type="samples", color="Actinobacteriota")+
  geom_point(size=1) +
  scale_color_gradient2(midpoint=mean(sample_data(physeq.fecal)$Bacteroidota), low="blue", mid="white", high="red")+
  theme(panel.grid = element_blank(),
        panel.background=element_blank())
ggsave("~/IBS/PCoA/data/pcoa_bray_actinobacteriota.jpg", width=6, height=4, type="cairo")

plot_ordination(physeq.CSN, ord, type="samples", color="Proteobacteria")+
  geom_point(size=1) +
  scale_color_gradient2(midpoint=mean(sample_data(physeq.fecal)$Bacteroidota), low="blue", mid="white", high="red")+
  theme(panel.grid = element_blank(),
        panel.background=element_blank())
ggsave("~/IBS/PCoA/data/pcoa_bray_proteobacteria.jpg", width=6, height=4, type="cairo")


