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
           'author', 'sequencing_tech', 'bowel_movement_quality','bowel_movement_frequency', 'Bristol')) %>%
  rename(Bacteroidota = Abundance) %>%
  select(-Phylum)
firmi <- phylumTable %>%
  filter(Phylum == "Firmicutes") %>%
  select(c('Sample', 'Abundance', 'Phylum', 'host_disease', 'host_subtype', 'sample_type', 'Collection',
           'author', 'sequencing_tech', 'bowel_movement_quality','bowel_movement_frequency', 'Bristol')) %>%
  rename(Firmicutes = Abundance) %>%
  select(-Phylum)

# Check if there are any count 0 (would prevent from calculating ratio)
table(bacter$Bacteroidota == 0)
table(firmi$Firmicutes == 0)
min(bacter[bacter$Bacteroidota > 0, "Bacteroidota"]) # 3
min(firmi[bacter$Firmicutes > 0, "Firmicutes"]) # 40

# Compute log ratio (add 0.5 pseudocounts for the few 0 counts in Firmicutes and Bacteroidota)
ratioFB <- left_join(x=bacter, y=firmi,
                  by=c("Sample", "host_disease", "host_subtype", "sample_type", "Collection", "author", "sequencing_tech",
                       "bowel_movement_quality", "bowel_movement_frequency","Bristol")) %>%
  # mutate(Firmicutes=replace(Firmicutes, Firmicutes==0, 0.5),
  #        Bacteroidota=replace(Bacteroidota, Bacteroidota==0, 0.5)) %>%
  mutate(author_disease=paste0(author, sep="_", host_disease),
         LogRatioFB = log2(Firmicutes/Bacteroidota)) %>%
  relocate(Bacteroidota, .after=Firmicutes)

# Add IBS subtype to AGP
# table(ratioFB$host_subtype)
ratioFB <- ratioFB %>%
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
# ratioFB %>%
#   filter(author=="AGP") %>%
#   group_by(host_subtype, bowel_movement_quality) %>%
#   count()


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
ratioFBA$author_disease <- factor(ratioFBA$author_disease, levels=author_disease.order)

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
ggplot(ratioFBA %>% filter(Collection=="1st"),
       aes(x = author_disease, y = LogRatioFA,
       color=factor(author, levels = author.order)))+
  # geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.4, lwd=1)+
  geom_jitter(width = 0.1, size = 0.3)+
  theme_classic()+
  scale_x_discrete(breaks=author_disease.order, labels= rep(c("Healthy", "IBS"), times = 12))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))+
  labs(x = '', y = "Log2(Firmicutes/Actino)", color="author", title="All samples (sigmoid & stool")
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
ggplot(ratioFBA %>% filter(sample_type == "stool", Collection=="1st"),
       aes(x = host_disease, y = LogRatioBA, fill=host_disease))+
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
ggplot(ratioFB %>% filter(#Collection=="1st",
                          host_subtype != "IBS-unspecified", host_subtype != "HC-unknown",
                          #sample_type == "stool",
                          author %in% c("Labus", "LoPresti", "Zhuang", "Zeber-Lubecka", "Nagel", "Mars", "Liu", "AGP")),
       aes(x = host_subtype, y = LogRatioFB))+
  # geom_violin(alpha = 0.5)+
  facet_wrap(~author, scales="fixed")+
  ylim(c(-3,7))+
  geom_boxplot(lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color = factor(author, levels = author.order)), width = 0.05, size=0.6)+
  labs(title="",
       x = '', y = "Log2(Firmicutes/Bacteroidota)",
       color="author")+
  theme_classic()+
  theme(axis.text = element_text(size=10, color="black"),
        axis.title = element_text(size=10)
        #axis.line = element_line(arrow = arrow(length = unit(0.1, "inches")))
        )
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/06_FirmBactActin-LogRatio/firm_bact_06.jpg", width=5, height=4)

#XXXXXXXXX TO REMOVE XXXXXXXXX
# Liu Bristol stool scale
ggplot(ratioFB %>% filter(author == "Liu") %>% mutate(LogRatioBF = log2(Bacteroidota/Firmicutes)),
       aes(x = Bristol, y = LogRatioBF))+
  # geom_violin(alpha = 0.5)+
  geom_boxplot(aes(group=Bristol), lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color=host_disease), alpha=0.5, width = 0.05, size=3)+
  scale_color_manual(name = "disease", values=c("blue", "red"))+
  # scale_fill_manual(values=c("blue", "red"))+
  labs(title="Liu dataset (2020)",
       x = 'Bristol stool scale', y = "Log2(Bacteroidota/Firmicutes)")+
  theme_classic()+
  theme(axis.text = element_text(size=20, color="black"),
        axis.title = element_text(size=20),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))

# AGP
ggplot(ratioFB %>% filter(author == "AGP" & bowel_movement_quality != "Unknown") %>%
         filter((host_disease=="Healthy" & bowel_movement_quality=="Normal") | 
                  (host_disease=="IBS" & bowel_movement_quality=="Constipated") |
                  (host_disease=="IBS" & bowel_movement_quality=="Diarrhea")) %>%
         mutate(LogRatioBF = log2(Bacteroidota/Firmicutes),
                #bowel_movement_frequency=factor(bowel_movement_frequency, levels=c("Less than one", "One", "Two", "Three", "Four", "Five or more"))),
                bowel_movement_quality=factor(bowel_movement_quality, levels=c("Normal", "Constipated", "Diarrhea"))
                ),
       aes(x = bowel_movement_quality, y = LogRatioFB))+
  #geom_violin(alpha = 0.5)+
  geom_boxplot(aes(group=bowel_movement_quality), lwd = 1, width=0.3, outlier.shape=NA, fill = "white")+
  geom_jitter(aes(color=host_disease), alpha=0.5, width = 0.05, size=2)+
  scale_color_manual(name = "disease", values=c("blue", "red"))+
  # scale_fill_manual(values=c("blue", "red"))+
  labs(title="AGP dataset (2021)",
       x = 'bowel_movement_quality', y = "Log2(Firmicutes/Bacteroidota)")+
  theme_classic()+
  theme(axis.text = element_text(size=20, color="black"),
        axis.title = element_text(size=20),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))))
#XXXXXXXXX TO REMOVE XXXXXXXXX

##________________________________________
## STATISTICS


##############################################
# PLOT FIRMICUTES/ACTINOBACTERIOTA LOG RATIO #
##############################################

##________________________________________
## BUILD DATAFRAME
# Extract abundance of Bacteroidota and Firmicutes
actino <- phylumTable %>%
  filter(Phylum == "Actinobacteriota") %>%
  select(c('Sample', 'Abundance', 'Phylum', 'host_disease', 'host_subtype', 'sample_type', 'Collection',
           'author', 'sequencing_tech', 'bowel_movement_quality', 'Bristol')) %>%
  rename(Actinobacteriota = Abundance) %>%
  select(-Phylum)

# Compute log ratio
ratioFBA <- left_join(x=ratioFB, y=actino,
                     by=c("Sample", "host_disease", "host_subtype", "sample_type", "Collection", "author", "sequencing_tech", "bowel_movement_quality", "Bristol")) %>%
  mutate(LogRatioFA = log2(Firmicutes/Actinobacteriota),
         LogRatioBA = log2(Bacteroidota/Actinobacteriota)) %>%
  relocate(Actinobacteriota, .after=Bacteroidota)








################################################
# PLOT BACTEROIDOTA/ACTINOBACTERIOTA LOG RATIO #
################################################




