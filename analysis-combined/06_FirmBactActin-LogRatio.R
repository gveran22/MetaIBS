##########################
# Purpose: Plotting Firmicutes/Bacteroidota log ratio
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
nsamples(physeq.fecal)
physeq.sigmoid <- subset_samples(physeq, sample_type == 'sigmoid') # 431 samples
nsamples(physeq.sigmoid)

# Agglomerate to phylum level and melt to long format
phylumTable <- physeq %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt()




##########################################
# PLOT FIRMICUTES/BACTEROIDOTA LOG RATIO #
##########################################

##________________________________________
## BUILD DATAFRAME
# Extract abundance of Bacteroidota and Firmicutes
bacter <- phylumTable %>%
  filter(Phylum == "Bacteroidota") %>%
  select(c('Sample', 'Abundance', 'Phylum', 'host_disease', 'host_subtype', 'sample_type', 'Collection',
           'author', 'sequencing_tech', 'bowel_movement_quality')) %>%
  rename(Bacteroidota = Abundance) %>%
  select(-Phylum)
firmi <- phylumTable %>%
  filter(Phylum == "Firmicutes") %>%
  select(c('Sample', 'Abundance', 'Phylum', 'host_disease', 'host_subtype', 'sample_type', 'Collection',
           'author', 'sequencing_tech', 'bowel_movement_quality')) %>%
  rename(Firmicutes = Abundance) %>%
  select(-Phylum)

# Compute log ratio
ratioFB <- left_join(x=bacter, y=firmi,
                  by=c("Sample", "host_disease", "host_subtype", "sample_type", "Collection", "author", "sequencing_tech", "bowel_movement_quality")) %>%
  mutate(author_disease=paste0(author, sep="_", host_disease),
         LogRatioFB = log2(Firmicutes/Bacteroidota)) %>%
  relocate(Bacteroidota, .after=Firmicutes)


##________________________________________
## PLOT
# Set order
author.order <- c('Labus', 'LoPresti',
                   'Pozuelo', 'Zhuang', 'Zhu', 'Hugerth', 'Fukui', 'Mars', 'Liu', 'AGP',
                   'Nagel', 'Zeber-Lubecka')
author_disease.order <- c("Labus_Healthy", "Labus_IBS",
                          "LoPresti_Healthy", "LoPresti_IBS",
                          "Pozuelo_Healthy", "Pozuelo_IBS",
                          "Zhuang_Healthy", "Zhuang_IBS",
                          "Zhu_Healthy", "Zhu_IBS",
                          "Hugerth_Healthy", "Hugerth_IBS",
                          "Fukui_Healthy", "Fukui_IBS",
                          "Mars_Healthy", "Mars_IBS",
                          "Liu_Healthy", "Liu_IBS",
                          "AGP_Healthy", "AGP_IBS",
                          "Nagel_Healthy", "Nagel_IBS",
                          "Zeber-Lubecka_Healthy", "Zeber-Lubecka_IBS")
ratioFB$author_disease <- factor(ratioFB$author_disease, levels=author_disease.order)

# Plot Firm/Bact ratio PER DATASET, HEALTHY/IBS, with FACET_WRAP
ggplot(ratioFB %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = host_disease, y = LogRatioFB,
           fill=factor(author, levels=author.order)))+
  facet_wrap(~factor(author, levels = author.order), scales="fixed")+
  geom_violin()+
  # geom_boxplot(outlier.size = 0, width = 0.25)+
  geom_jitter(width = 0.1, size = 0.3)+
  theme_classic()+
  scale_fill_brewer(palette="Set3")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))),
        legend.position="None")+
  labs(x = '', y = "Log2(Firmicutes/Bacteroidota)")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_01.jpg", width=8, height=6)

# Plot Firm/Bact ratio PER DATASET, HEALTHY/IBS, with BOXPLOTS
ggplot(ratioFB %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = author_disease, y = LogRatioFB,
       color=factor(author, levels = authors.order)))+
  # geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.4, lwd=1)+
  geom_jitter(width = 0.1, size = 0.3)+
  theme_classic()+
  scale_x_discrete(breaks=author_disease.order, labels= rep(c("Healthy", "IBS"), times = 12))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))+
  labs(x = '', y = "Log2(Firmicutes/Bacteroidota)", color="author")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_02.jpg", width=12, height=5)

# Plot Firm/Bact ratio PER SEQUENCING TECH, HEALTHY/IBS
ggplot(ratioFB %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = sequencing_tech, y = LogRatioFB,
           fill=host_disease))+
  geom_violin(alpha=0.5)+
  geom_boxplot(aes(group=interaction(host_disease,sequencing_tech)), position = position_dodge(0.9),
               lwd = 1, width = 0.2, fill = "white", outlier.shape=NA)+
  theme_classic()+
  scale_fill_manual(values = c("blue", "red"))+
  theme(axis.text.x = element_text(size=15, color="black", angle=35, hjust=1),
        axis.title = element_text(size=15),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))+
  labs(x = '', y = "Log2(Firmicutes/Bacteroidota)", fill="Disease")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_03.jpg", width=12, height=5)

# Plot Firm/Bact ratio pooled together (HEALTHY/IBS)
ggplot(ratioFB %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = host_disease, y = LogRatioFB, fill=host_disease))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(lwd = 1, width = 0.2, fill = "white")+
  # geom_jitter(width = 0.1, size = 0.3)+
  scale_fill_manual(values=c("blue", "red"))+
  labs(x = '', y = "Log2(Firmicutes/Bacteroidota)")+
  theme_classic()+
  theme(axis.text = element_text(size=15, color="black"),
        axis.title = element_text(size=15),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))),
        legend.position="None")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_04.jpg", width=5, height=4)

# Plot Firm/Bact ratio pooled together (IBS SUBTYPE)
ggplot(ratioFB %>% filter(Collection=="1st", host_subtype != "IBS-unspecified",
                          sample_type == "stool",
                          author %in% c("Labus", "LoPresti", "Zhuang", "Zeber-Lubecka", "Nagel", "Mars")),
       aes(x = host_subtype, y = LogRatioFB))+
  # geom_violin(alpha = 0.5)+
  geom_boxplot(lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color = factor(author, levels = authors_order)), width = 0.05, size=0.6)+
  scale_color_discrete(name = "author")+
  # scale_fill_manual(values=c("blue", "red"))+
  labs(title="Stool samples only",
       x = '', y = "Log2(Firmicutes/Bacteroidota)")+
  theme_classic()+
  theme(axis.text = element_text(size=10, color="black"),
        axis.title = element_text(size=10),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_06.jpg", width=5, height=4)

##________________________________________
## STATISTICS


##############################################
# PLOT FIRMICUTES/ACTINOBACTERIOTA LOG RATIO #
##############################################










################################################
# PLOT BACTEROIDOTA/ACTINOBACTERIOTA LOG RATIO #
################################################




