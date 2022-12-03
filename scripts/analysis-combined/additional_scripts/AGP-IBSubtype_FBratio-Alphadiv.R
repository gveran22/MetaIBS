##########################
# Purpose: Plot Firm/Bact ratio and Alpha-div in AGP dataset (with IBS subtypes)
# Date: April 2022
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggsignif)

# Data
physeq.agp <- readRDS("~/Projects/IBS_Meta-analysis_16S/phyloseq-objects/physeq_agp.rds")

# Agglomerate to phylum level and melt to long format
phylumTable <- physeq.agp %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt()




#################################################
# PLOT FIRMICUTES/BACTEROIDOTA LOG RATIO IN AGP #
#################################################

# Extract abundance of Bacteroidota, Firmicutes and Actinobacteriota
relevant.covariates <- c('Sample', 'Abundance', 'Phylum', 'host_disease', 'host_subtype',
                         'bowel_movement_quality','bowel_movement_frequency', 'host_age', 'host_bmi')

bacter <- phylumTable %>%
  filter(Phylum == "Bacteroidota") %>%
  select(all_of(relevant.covariates)) %>%
  rename(Bacteroidota = Abundance) %>%
  select(-Phylum)

firmi <- phylumTable %>%
  filter(Phylum == "Firmicutes") %>%
  select(all_of(relevant.covariates)) %>%
  rename(Firmicutes = Abundance) %>%
  select(-Phylum)


# Check if there are any count 0 (would prevent from calculating ratio)
table(bacter$Bacteroidota == 0) # 1
table(firmi$Firmicutes == 0) # 0
min(bacter[bacter$Bacteroidota > 0, "Bacteroidota"]) # 11
min(firmi[firmi$Firmicutes > 0, "Firmicutes"]) # 105

# Sanity check (should have 1,183 samples)
# nrow(bacter)
# nrow(firmi)


# COMPUTE LOG RATIOS
common.columns <- relevant.covariates[!relevant.covariates %in% c("Abundance", "Phylum")]

ratio.df <- left_join(x=bacter, y=firmi, by=common.columns) %>%
  relocate(Bacteroidota, .before=Firmicutes) %>%
  # Add 0.5 pseudocounts for the few 0 values
  mutate(Bacteroidota=replace(Bacteroidota, Bacteroidota==0, 0.5),
         Firmicutes=replace(Firmicutes, Firmicutes==0, 0.5)) %>%
  # Compute log ratios
  mutate(LogRatio_FirmBact = log2(Firmicutes/Bacteroidota)) %>%
  # Add IBS subtype based on bowel mvt
  mutate(host_subtype=replace(host_subtype,
                              host_disease=="IBS" & bowel_movement_quality=="Constipated",
                              "IBS-C")) %>%
  mutate(host_subtype=replace(host_subtype,
                              host_disease=="IBS" & bowel_movement_quality=="Diarrhea",
                              "IBS-D")) %>%
  mutate(host_subtype=replace(host_subtype,
                              host_disease=="Healthy" & bowel_movement_quality=="Normal",
                              "Healthy")) # healthy samples with non-normal bowel_mvt will stay "HC" so we will exclude them in the plot

# sanity check
# ratio.df %>%
#   group_by(host_disease, bowel_movement_quality, host_subtype) %>%
#   count()

ggplot(data=ratio.df %>% filter(host_subtype %in% c("Healthy", "IBS-C", "IBS-D")),
       aes(x = host_subtype, y = LogRatio_FirmBact))+
  geom_jitter(aes(color=host_subtype), width=0.05, size=0.5)+
  geom_boxplot(aes(fill=host_subtype), lwd=0.5, width=0.4, outlier.shape=NA, alpha=0.2)+
  geom_signif(comparisons = list(c("Healthy", "IBS-C"), c("Healthy", "IBS-D")), map_signif_level = TRUE, test="wilcox.test") +
  scale_color_manual(values=c("#3182bd", "#a50f15", "#fcae91", "#fb6a4a"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#a50f15", "#fcae91", "#fb6a4a"), guide="none")+
  # THEME
  theme_cowplot()+
  theme(strip.text.y = element_text(angle=0, hjust=0),
        strip.background = element_blank())+
  labs(x = '', y =  expression(Log[2](Firmicutes/Bacteroidota)))




#############################
# PLOT SHANNON INDEX IN AGP #
#############################

# Get Shannon values
plt.shannon <- plot_richness(physeq.agp, measures="Shannon")
shannon.df <- plt.shannon$data %>%
  select(c("samples", "value", "host_disease", "host_subtype", "bowel_movement_quality")) %>%
  dplyr::rename(shannon=value) %>%
  # Add IBS subtype based on bowel mvt
  mutate(host_subtype=replace(host_subtype,
                              host_disease=="IBS" & bowel_movement_quality=="Constipated",
                              "IBS-C")) %>%
  mutate(host_subtype=replace(host_subtype,
                              host_disease=="IBS" & bowel_movement_quality=="Diarrhea",
                              "IBS-D")) %>%
  mutate(host_subtype=replace(host_subtype,
                              host_disease=="Healthy" & bowel_movement_quality=="Normal",
                              "Healthy")) 
# sanity check
# shannon.df %>%
#   group_by(host_disease, bowel_movement_quality, host_subtype) %>%
#   count()

ggplot(data=shannon.df %>% filter(host_subtype %in% c("Healthy", "IBS-C", "IBS-D")),
       aes(x = host_subtype, y = shannon))+
  geom_jitter(aes(color=host_subtype), width=0.05, size=0.5)+
  geom_boxplot(aes(fill=host_subtype), lwd=0.5, width=0.4, outlier.shape=NA, alpha=0.2)+
  geom_signif(comparisons = list(c("Healthy", "IBS-C"), c("Healthy", "IBS-D")), map_signif_level = TRUE, test="wilcox.test") +
  scale_color_manual(values=c("#3182bd", "#a50f15", "#fcae91", "#fb6a4a"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#a50f15", "#fcae91", "#fb6a4a"), guide="none")+
  # THEME
  theme_cowplot()+
  theme(strip.text.y = element_text(angle=0, hjust=0),
        strip.background = element_blank())+
  labs(x = '', y = "Shannon index")


