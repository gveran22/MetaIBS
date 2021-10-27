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
physeq.fecal <- subset_samples(physeq, sample_type == 'stool') # 2,145 samples
nsamples(physeq.fecal)
physeq.sigmoid <- subset_samples(physeq, sample_type == 'sigmoid') # 431 samples
nsamples(physeq.sigmoid)

# Agglomerate to phylum level and melt to long format
phylumTable <- physeq %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt()




#########################################
# COMPUTE FIRM, BACT, ACTINO LOG RATIOS #
#########################################

# Extract abundance of Bacteroidota, Firmicutes and Actinobacteriota
relevant.covariates <- c('Sample', 'Abundance', 'Phylum', 'host_disease', 'host_subtype', 'sample_type', 'Collection',
                         'author', 'sequencing_tech', 'bowel_movement_quality','bowel_movement_frequency', 'Bristol')

bacter <- phylumTable %>%
  filter(Phylum == "Bacteroidota") %>%
  select(relevant.covariates) %>%
  rename(Bacteroidota = Abundance) %>%
  select(-Phylum)

firmi <- phylumTable %>%
  filter(Phylum == "Firmicutes") %>%
  select(relevant.covariates) %>%
  rename(Firmicutes = Abundance) %>%
  select(-Phylum)

actino <- phylumTable %>%
  filter(Phylum == "Actinobacteriota") %>%
  select(relevant.covariates) %>%
  rename(Actinobacteriota = Abundance) %>%
  select(-Phylum)

proteo <- phylumTable %>%
  filter(Phylum == "Proteobacteria") %>%
  select(relevant.covariates) %>%
  rename(Proteobacteria = Abundance) %>%
  select(-Phylum)


# Check if there are any count 0 (would prevent from calculating ratio)
table(bacter$Bacteroidota == 0) # 18
table(firmi$Firmicutes == 0) # 3
table(actino$Actinobacteriota == 0) # 150
table(proteo$Proteobacteria == 0) # 116
min(bacter[bacter$Bacteroidota > 0, "Bacteroidota"]) # 3
min(firmi[firmi$Firmicutes > 0, "Firmicutes"]) # 40
min(actino[actino$Actinobacteriota > 0, "Actinobacteriota"]) # 2
min(proteo[proteo$Proteobacteria > 0, "Proteobacteria"]) # 2


# COMPUTE LOG RATIOS
common.columns <- c("Sample", "host_disease", "host_subtype", "sample_type", "Collection",
                    "author", "sequencing_tech", "bowel_movement_quality", "bowel_movement_frequency", "Bristol")

ratio.df <- left_join(x=bacter, y=firmi, by=common.columns) %>%
  left_join(actino, by=common.columns) %>%
  left_join(proteo, by=common.columns) %>%
  relocate(Bacteroidota, .before=Firmicutes) %>%
  # Add 0.5 pseudocounts for the few 0 values
  mutate(Bacteroidota=replace(Bacteroidota, Bacteroidota==0, 0.5),
         Firmicutes=replace(Firmicutes, Firmicutes==0, 0.5),
         Actinobacteriota=replace(Actinobacteriota, Actinobacteriota==0, 0.5),
         Proteobacteria=replace(Proteobacteria, Proteobacteria==0, 0.5)) %>%
  # Compute log ratios
  mutate(LogRatio_FirmBact = log2(Firmicutes/Bacteroidota),
         LogRatio_FirmAct = log2(Firmicutes/Actinobacteriota),
         LogRatio_BactAct = log2(Bacteroidota/Actinobacteriota),
         LogRatio_FirmProt = log2(Firmicutes/Proteobacteria),
         LogRatio_BactProt = log2(Bacteroidota/Proteobacteria),
         LogRatio_ActProt = log2(Actinobacteriota/Proteobacteria)) %>%
  # Add column author_disease
  mutate(author_disease=paste0(author, sep="_", host_disease)) %>%
  # Add IBS subtype to AGP
  mutate(host_subtype=replace(host_subtype,
                              author=="AGP" & host_disease=="IBS" & bowel_movement_quality=="Constipated",
                              "IBS-C")) %>%
  mutate(host_subtype=replace(host_subtype,
                              author=="AGP" & host_disease=="IBS" & bowel_movement_quality=="Diarrhea",
                              "IBS-D")) %>%
  mutate(host_subtype=replace(host_subtype,
                              author=="AGP" & host_disease=="Healthy" & bowel_movement_quality!="Normal",
                              "HC-unknown"))

# sanity check
# ratio.df %>%
#   filter(author=="AGP") %>%
#   group_by(host_subtype, bowel_movement_quality) %>%
#   count()


# Set dataset order for plots
author.order <- c('Labus', 'LoPresti', # 454 pyrosequencing
                  'AGP', 'Liu', 'Pozuelo', # Illumina single end
                  'Fukui', 'Hugerth', 'Mars', 'Zhu', 'Zhuang', # Illumina paired end
                  'Nagel', 'Zeber-Lubecka') # Ion Torrent
ratio.df$author <- factor(ratio.df$author, levels=author.order)

author_disease.order <- c("Labus_Healthy", "Labus_IBS",
                          "LoPresti_Healthy", "LoPresti_IBS",
                          # Illumina single end
                          "AGP_Healthy", "AGP_IBS",
                          "Liu_Healthy", "Liu_IBS",
                          "Pozuelo_Healthy", "Pozuelo_IBS",
                          # Illumina paired end
                          "Fukui_Healthy", "Fukui_IBS",
                          "Hugerth_Healthy", "Hugerth_IBS",
                          "Mars_Healthy", "Mars_IBS",
                          "Zhu_Healthy", "Zhu_IBS",
                          "Zhuang_Healthy", "Zhuang_IBS",
                          # Ion Torrent
                          "Nagel_Healthy", "Nagel_IBS",
                          "Zeber-Lubecka_Healthy", "Zeber-Lubecka_IBS")
ratio.df$author_disease <- factor(ratio.df$author_disease, levels=author_disease.order)

seqtech.order <- c("454 pyrosequencing", "Illumina single-end", "Illumina paired-end", "Ion Torrent")
ratio.df$sequencing_tech <- factor(ratio.df$sequencing_tech, levels=seqtech.order)


##########################################
# PLOT FIRMICUTES/BACTEROIDOTA LOG RATIO #
##########################################

# Plot Firm/Bact ratio PER DATASET, HEALTHY/IBS, with FACET_WRAP
ggplot(ratio.df %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = host_disease, y = LogRatio_FirmBact,
           fill=factor(author, levels=author.order)))+
  facet_wrap(~factor(author, levels = author.order), scales="fixed")+
  geom_violin()+
  geom_jitter(width = 0.1, size = 0.3)+
  theme_classic()+
  scale_fill_brewer(palette="Set3")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))),
        legend.position="None")+
  labs(x = '', y = "Log2(Firmicutes/Bacteroidota)")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_01.jpg", width=8, height=6)


# Plot Firm/Bact ratio PER DATASET, HEALTHY/IBS, with BOXPLOTS
ggplot(ratio.df %>% filter(Collection=="1st"),
       aes(x = author_disease, y = LogRatio_FirmBact,
       color=factor(author, levels = author.order)))+
  geom_boxplot(outlier.shape = NA, width = 0.4, lwd=1)+
  geom_jitter(width = 0.1, size = 0.3)+
  theme_classic()+
  scale_x_discrete(breaks=author_disease.order, labels= rep(c("Healthy", "IBS"), times = 12))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))+
  labs(x = '', y = "Log2(Firmicutes/Bacteroidota)", color="author", title="All samples (sigmoid & stool)")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_02.jpg", width=12, height=5)


# Plot Firm/Bact ratio PER SEQUENCING TECH, HEALTHY/IBS
ggplot(ratio.df %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = sequencing_tech, y = LogRatio_FirmBact,
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


# Plot Firm/Bact ratio POOLED together (HEALTHY/IBS)
ggplot(ratio.df %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = host_disease, y = LogRatio_FirmBact, fill=host_disease))+
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


# Plot Firm/Bact ratio POOLED together by IBS SUBTYPE
ggplot(ratio.df %>% filter(sample_type == "stool", Collection=="1st",
                           host_subtype != "IBS-unspecified", host_subtype != "HC-unknown",
                           author %in% c("Labus", "LoPresti", "Liu", "Pozuelo", "Mars", "Zhuang", "Nagel", "Zeber-Lubecka")),
       aes(x = host_subtype, y = LogRatio_FirmBact))+
  geom_boxplot(lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color = factor(author, levels = author.order)), width = 0.05, size=0.6)+
  labs(title="Stool samples only",
       x = '', y = "Log2(Firmicutes/Bacteroidota)",
       color="author")+
  theme_classic()+
  theme(axis.text = element_text(size=10, color="black"),
        axis.title = element_text(size=10))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_05.jpg", width=5, height=4)


# Plot Firm/Bact ratio POOLED together by IBS SUBTYPE, with FACET_WRAP
ggplot(ratio.df %>% filter(Collection=="1st",
                           host_subtype != "IBS-unspecified", host_subtype != "HC-unknown",
                           author %in% c("Labus", "LoPresti", "Liu", "Pozuelo", "Mars", "Zhuang", "Nagel", "Zeber-Lubecka")),
       aes(x = host_subtype, y = LogRatio_FirmBact))+
  geom_violin()+
  facet_wrap(~author, scales="fixed", ncol=2)+
  ylim(c(-3,10))+
  #geom_boxplot(lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color = factor(author, levels = author.order)), width = 0.05, size=2, alpha=0.8)+
  labs(title="All sample types (stool & sigmoid)",
       x = '', y = "Log2(Firmicutes/Bacteroidota)",
       color="author")+
  theme_classic()+
  theme(axis.text = element_text(size=10, color="black"),
        axis.title = element_text(size=10))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_06.jpg", width=6, height=8)


# LIU BRISTOL STOOL SCALE (bacteroidota/firmicutes)
ggplot(ratio.df %>% filter(author == "Liu") %>% mutate(LogRatioBF = log2(Bacteroidota/Firmicutes)),
       aes(x = Bristol, y = LogRatioBF))+
  # geom_violin(alpha = 0.5)+
  geom_boxplot(aes(group=Bristol), lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color=host_disease), alpha=0.5, width = 0.05, size=3)+
  scale_color_manual(name = "disease", values=c("blue", "red"))+
  labs(title="Liu dataset (2020)",
       x = 'Bristol stool scale', y = "Log2(Bacteroidota/Firmicutes)")+
  theme_classic()+
  theme(axis.text = element_text(size=15, color="black"),
        axis.title = element_text(size=15),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_07.jpg", width=6, height=6)


# AGP bowel_movement_quality
ggplot(ratio.df %>% filter(author == "AGP" & host_subtype!="HC-unknown" & host_subtype!="IBS-unspecified") %>%
                    mutate(bowel_movement_quality=factor(bowel_movement_quality, levels=c("Normal", "Constipated", "Diarrhea"))),
       aes(x = bowel_movement_quality, y = LogRatio_FirmBact))+
  #geom_violin(alpha = 0.5)+
  geom_boxplot(aes(group=bowel_movement_quality), lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color=host_disease), alpha=0.5, width = 0.05, size=2)+
  scale_color_manual(name = "disease", values=c("blue", "red"))+
  labs(title="AGP dataset (2021)",
       x = 'Bowel movement quality', y = "Log2(Firmicutes/Bacteroidota)")+
  theme_classic()+
  theme(axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=15))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_08.jpg", width=6, height=6)


# AGP bowel_movement_frequency
freq.order <- c("Less than one", "One", "Two", "Three", "Four", "Five or more")
ggplot(ratio.df %>% filter(author == "AGP" & bowel_movement_frequency %in% freq.order) %>%
                    mutate(bowel_movement_frequency=factor(bowel_movement_frequency, levels=freq.order)),
       aes(x = bowel_movement_frequency, y = LogRatio_FirmBact))+
  #geom_violin(alpha = 0.5)+
  geom_boxplot(aes(group=bowel_movement_frequency), lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color=factor(bowel_movement_quality, levels=c("Normal", "Constipated", "Diarrhea"))), alpha=0.5, width = 0.05, size=2)+
  scale_color_manual(values=c("#6699FF", "#CC0000","#CC33FF"))+
  labs(title="AGP dataset (2021)", color="Bowel mvt quality",
       x = 'Bowel movement frequency (per day)', y = "Log2(Firmicutes/Bacteroidota)")+
  theme_classic()+
  theme(axis.text = element_text(size=12, color="black", angle=40, hjust=1),
        axis.title = element_text(size=15))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_09.jpg", width=6, height=6)




##############################################
# PLOT FIRMICUTES/ACTINOBACTERIOTA LOG RATIO #
##############################################

# Plot Firm/Actino ratio PER DATASET, HEALTHY/IBS, with FACET_WRAP
ggplot(ratio.df %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = host_disease, y = LogRatio_FirmAct,
           fill=factor(author, levels=author.order)))+
  facet_wrap(~factor(author, levels = author.order), scales="fixed")+
  geom_violin()+
  geom_jitter(width = 0.1, size = 0.3)+
  theme_classic()+
  scale_fill_brewer(palette="Set3")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))),
        legend.position="None")+
  labs(x = '', y = "Log2(Firmicutes/Actinobacteriota)")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_actino_01.jpg", width=8, height=6)


# Plot Firm/Actino ratio PER DATASET, HEALTHY/IBS, with BOXPLOTS
ggplot(ratio.df %>% filter(Collection=="1st"),
       aes(x = author_disease, y = LogRatio_FirmAct,
           color=factor(author, levels = author.order)))+
  geom_boxplot(outlier.shape = NA, width = 0.4, lwd=1)+
  geom_jitter(width = 0.1, size = 0.3)+
  theme_classic()+
  scale_x_discrete(breaks=author_disease.order, labels= rep(c("Healthy", "IBS"), times = 12))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))+
  labs(x = '', y = "Log2(Firmicutes/Actinobacteriota)", color="author", title="All samples (sigmoid & stool)")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_actino_02.jpg", width=12, height=5)


# Plot Firm/Actino ratio PER SEQUENCING TECH, HEALTHY/IBS
ggplot(ratio.df %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = sequencing_tech, y = LogRatio_FirmAct,
           fill=host_disease))+
  geom_violin(alpha=0.5)+
  geom_boxplot(aes(group=interaction(host_disease,sequencing_tech)), position = position_dodge(0.9),
               lwd = 1, width = 0.2, fill = "white", outlier.shape=NA)+
  theme_classic()+
  scale_fill_manual(values = c("blue", "red"))+
  theme(axis.text.x = element_text(size=15, color="black", angle=35, hjust=1),
        axis.title = element_text(size=15),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))+
  labs(x = '', y = "Log2(Firmicutes/Actinobacteriota)", fill="Disease")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_actino_03.jpg", width=12, height=5)


# Plot Firm/Actino ratio POOLED together (HEALTHY/IBS)
ggplot(ratio.df %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = host_disease, y = LogRatio_FirmAct, fill=host_disease))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(lwd = 1, width = 0.2, fill = "white")+
  # geom_jitter(width = 0.1, size = 0.3)+
  scale_fill_manual(values=c("blue", "red"))+
  labs(x = '', y = "Log2(Firmicutes/Actinobacteriota)")+
  theme_classic()+
  theme(axis.text = element_text(size=15, color="black"),
        axis.title = element_text(size=15),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))),
        legend.position="None")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_actino_04.jpg", width=5, height=4)


# Plot Firm/Actino ratio POOLED together by IBS SUBTYPE
ggplot(ratio.df %>% filter(sample_type == "stool", Collection=="1st",
                           host_subtype != "IBS-unspecified", host_subtype != "HC-unknown",
                           author %in% c("Labus", "LoPresti", "Liu", "Pozuelo", "Mars", "Zhuang", "Nagel", "Zeber-Lubecka")),
       aes(x = host_subtype, y = LogRatio_FirmAct))+
  geom_boxplot(lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color = factor(author, levels = author.order)), width = 0.05, size=0.6)+
  labs(title="Stool samples only",
       x = '', y = "Log2(Firmicutes/Actinobacteriota)",
       color="author")+
  theme_classic()+
  theme(axis.text = element_text(size=10, color="black"),
        axis.title = element_text(size=10))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_actino_05.jpg", width=5, height=4)


# Plot Firm/Actino ratio POOLED together by IBS SUBTYPE, with FACET_WRAP
ggplot(ratio.df %>% filter(Collection=="1st",
                           host_subtype != "IBS-unspecified", host_subtype != "HC-unknown",
                           author %in% c("Labus", "LoPresti", "Liu", "Pozuelo", "Mars", "Zhuang", "Nagel", "Zeber-Lubecka")),
       aes(x = host_subtype, y = LogRatio_FirmAct))+
  geom_violin()+
  facet_wrap(~author, scales="fixed", ncol=2)+
  #ylim(c(0,15))+
  #geom_boxplot(lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color = factor(author, levels = author.order)), width = 0.05, size=2, alpha=0.8)+
  labs(title="All sample types (stool & sigmoid)",
       x = '', y = "Log2(Firmicutes/Actinobacteriota)",
       color="author")+
  theme_classic()+
  theme(axis.text = element_text(size=10, color="black"),
        axis.title = element_text(size=10))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_actino_06.jpg", width=6, height=8)


# LIU BRISTOL STOOL SCALE (Firmicutes/Actinobacteriota)
ggplot(ratio.df %>% filter(author == "Liu"),
       aes(x = Bristol, y = LogRatio_FirmAct))+
  # geom_violin(alpha = 0.5)+
  geom_boxplot(aes(group=Bristol), lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color=host_disease), alpha=0.5, width = 0.05, size=3)+
  scale_color_manual(name = "disease", values=c("blue", "red"))+
  labs(title="Liu dataset (2020)",
       x = 'Bristol stool scale', y = "Log2(Firmicutes/Actinobacteriota)")+
  theme_classic()+
  theme(axis.text = element_text(size=15, color="black"),
        axis.title = element_text(size=15),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_actino_07.jpg", width=6, height=6)


# AGP bowel_movement_quality
ggplot(ratio.df %>% filter(author == "AGP" & host_subtype!="HC-unknown" & host_subtype!="IBS-unspecified") %>%
         mutate(bowel_movement_quality=factor(bowel_movement_quality, levels=c("Normal", "Constipated", "Diarrhea"))),
       aes(x = bowel_movement_quality, y = LogRatio_FirmAct))+
  #geom_violin(alpha = 0.5)+
  geom_boxplot(aes(group=bowel_movement_quality), lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color=host_disease), alpha=0.5, width = 0.05, size=2)+
  scale_color_manual(name = "disease", values=c("blue", "red"))+
  labs(title="AGP dataset (2021)",
       x = 'Bowel movement quality', y = "Log2(Firmicutes/Actinobacteriota)")+
  theme_classic()+
  theme(axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=15))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_actino_08.jpg", width=6, height=6)


# AGP bowel_movement_frequency
ggplot(ratio.df %>% filter(author == "AGP" & bowel_movement_frequency %in% freq.order) %>%
         mutate(bowel_movement_frequency=factor(bowel_movement_frequency, levels=freq.order)),
       aes(x = bowel_movement_frequency, y = LogRatio_FirmAct))+
  #geom_violin(alpha = 0.5)+
  geom_boxplot(aes(group=bowel_movement_frequency), lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color=factor(bowel_movement_quality, levels=c("Normal", "Constipated", "Diarrhea"))), alpha=0.5, width = 0.05, size=2)+
  scale_color_manual(values=c("#6699FF", "#CC0000","#CC33FF"))+
  labs(title="AGP dataset (2021)", color="Bowel mvt quality",
       x = 'Bowel movement frequency (per day)', y = "Log2(Firmicutes/Actinobacteriota)")+
  theme_classic()+
  theme(axis.text = element_text(size=12, color="black", angle=40, hjust=1),
        axis.title = element_text(size=15))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_actino_09.jpg", width=6, height=6)




################################################
# PLOT BACTEROIDOTA/ACTINOBACTERIOTA LOG RATIO #
################################################

# Plot Bacter/Actino ratio PER DATASET, HEALTHY/IBS, with FACET_WRAP
ggplot(ratio.df %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = host_disease, y = LogRatio_BactAct,
           fill=factor(author, levels=author.order)))+
  facet_wrap(~factor(author, levels = author.order), scales="fixed")+
  geom_violin()+
  geom_jitter(width = 0.1, size = 0.3)+
  theme_classic()+
  scale_fill_brewer(palette="Set3")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))),
        legend.position="None")+
  labs(x = '', y = "Log2(Bacteroidota/Actinobacteriota)")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/bacter_actino_01.jpg", width=8, height=6)


# Plot Bacter/Actino ratio PER DATASET, HEALTHY/IBS, with BOXPLOTS
ggplot(ratio.df %>% filter(Collection=="1st"),
       aes(x = author_disease, y = LogRatio_BactAct,
           color=factor(author, levels = author.order)))+
  geom_boxplot(outlier.shape = NA, width = 0.4, lwd=1)+
  geom_jitter(width = 0.1, size = 0.3)+
  theme_classic()+
  scale_x_discrete(breaks=author_disease.order, labels= rep(c("Healthy", "IBS"), times = 12))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))+
  labs(x = '', y = "Log2(Bacteroidota/Actinobacteriota)", color="author", title="All samples (sigmoid & stool)")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/bacter_actino_02.jpg", width=12, height=5)


# Plot Bacter/Actino ratio PER SEQUENCING TECH, HEALTHY/IBS
ggplot(ratio.df %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = sequencing_tech, y = LogRatio_BactAct,
           fill=host_disease))+
  geom_violin(alpha=0.5)+
  geom_boxplot(aes(group=interaction(host_disease,sequencing_tech)), position = position_dodge(0.9),
               lwd = 1, width = 0.2, fill = "white", outlier.shape=NA)+
  theme_classic()+
  scale_fill_manual(values = c("blue", "red"))+
  theme(axis.text.x = element_text(size=15, color="black", angle=35, hjust=1),
        axis.title = element_text(size=15),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))+
  labs(x = '', y = "Log2(Bacteroidota/Actinobacteriota)", fill="Disease")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/bacter_actino_03.jpg", width=12, height=5)


# Plot Bacter/Actino ratio POOLED together (HEALTHY/IBS)
ggplot(ratio.df %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = host_disease, y = LogRatio_BactAct, fill=host_disease))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(lwd = 1, width = 0.2, fill = "white")+
  # geom_jitter(width = 0.1, size = 0.3)+
  scale_fill_manual(values=c("blue", "red"))+
  labs(x = '', y = "Log2(Bacteroidota/Actinobacteriota)")+
  theme_classic()+
  theme(axis.text = element_text(size=15, color="black"),
        axis.title = element_text(size=15),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))),
        legend.position="None")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/bacter_actino_04.jpg", width=5, height=4)


# Plot Bacter/Actino ratio POOLED together by IBS SUBTYPE
ggplot(ratio.df %>% filter(sample_type == "stool", Collection=="1st",
                           host_subtype != "IBS-unspecified", host_subtype != "HC-unknown",
                           author %in% c("Labus", "LoPresti", "Liu", "Pozuelo", "Mars", "Zhuang", "Nagel", "Zeber-Lubecka")),
       aes(x = host_subtype, y = LogRatio_BactAct))+
  geom_boxplot(lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color = factor(author, levels = author.order)), width = 0.05, size=0.6)+
  labs(title="Stool samples only",
       x = '', y = "Log2(Bacteroidota/Actinobacteriota)",
       color="author")+
  theme_classic()+
  theme(axis.text = element_text(size=10, color="black"),
        axis.title = element_text(size=10))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/bacter_actino_05.jpg", width=5, height=4)


# Plot Bacter/Actino ratio POOLED together by IBS SUBTYPE, with FACET_WRAP
ggplot(ratio.df %>% filter(Collection=="1st",
                           host_subtype != "IBS-unspecified", host_subtype != "HC-unknown",
                           author %in% c("Labus", "LoPresti", "Liu", "Pozuelo", "Mars", "Zhuang", "Nagel", "Zeber-Lubecka")),
       aes(x = host_subtype, y = LogRatio_BactAct))+
  geom_violin()+
  facet_wrap(~author, scales="fixed", ncol=2)+
  #ylim(c(0,15))+
  #geom_boxplot(lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color = factor(author, levels = author.order)), width = 0.05, size=2, alpha=0.8)+
  labs(title="All sample types (stool & sigmoid)",
       x = '', y = "Log2(Bacteroidota/Actinobacteriota)",
       color="author")+
  theme_classic()+
  theme(axis.text = element_text(size=10, color="black"),
        axis.title = element_text(size=10))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/bacter_actino_06.jpg", width=6, height=8)


# LIU BRISTOL STOOL SCALE (Bacteroidota/Actinobacteriota)
ggplot(ratio.df %>% filter(author == "Liu"),
       aes(x = Bristol, y = LogRatio_BactAct))+
  # geom_violin(alpha = 0.5)+
  geom_boxplot(aes(group=Bristol), lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color=host_disease), alpha=0.5, width = 0.05, size=3)+
  scale_color_manual(name = "disease", values=c("blue", "red"))+
  labs(title="Liu dataset (2020)",
       x = 'Bristol stool scale', y = "Log2(Bacteroidota/Actinobacteriota)")+
  theme_classic()+
  theme(axis.text = element_text(size=15, color="black"),
        axis.title = element_text(size=15),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/bacter_actino_07.jpg", width=6, height=6)


# AGP bowel_movement_quality
ggplot(ratio.df %>% filter(author == "AGP" & host_subtype!="HC-unknown" & host_subtype!="IBS-unspecified") %>%
         mutate(bowel_movement_quality=factor(bowel_movement_quality, levels=c("Normal", "Constipated", "Diarrhea"))),
       aes(x = bowel_movement_quality, y = LogRatio_BactAct))+
  #geom_violin(alpha = 0.5)+
  geom_boxplot(aes(group=bowel_movement_quality), lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color=host_disease), alpha=0.5, width = 0.05, size=2)+
  scale_color_manual(name = "disease", values=c("blue", "red"))+
  labs(title="AGP dataset (2021)",
       x = 'Bowel movement quality', y = "Log2(Bacteroidota/Actinobacteriota)")+
  theme_classic()+
  theme(axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=15))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/bacter_actino_08.jpg", width=6, height=6)


# AGP bowel_movement_frequency
ggplot(ratio.df %>% filter(author == "AGP" & bowel_movement_frequency %in% freq.order) %>%
         mutate(bowel_movement_frequency=factor(bowel_movement_frequency, levels=freq.order)),
       aes(x = bowel_movement_frequency, y = LogRatio_BactAct))+
  #geom_violin(alpha = 0.5)+
  geom_boxplot(aes(group=bowel_movement_frequency), lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color=factor(bowel_movement_quality, levels=c("Normal", "Constipated", "Diarrhea"))), alpha=0.5, width = 0.05, size=2)+
  scale_color_manual(values=c("#6699FF", "#CC0000","#CC33FF"))+
  labs(title="AGP dataset (2021)", color="Bowel mvt quality",
       x = 'Bowel movement frequency (per day)', y = "Log2(Bacteroidota/Actinobacteriota)")+
  theme_classic()+
  theme(axis.text = element_text(size=12, color="black", angle=40, hjust=1),
        axis.title = element_text(size=15))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/bacter_actino_09.jpg", width=6, height=6)



