##########################
# Purpose: Plotting Firmicutes/Bacteroidota log ratio
# Date: December 2021
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
                         'author', 'sequencing_tech', 'bowel_movement_quality','bowel_movement_frequency', 'Bristol',
                         'host_age', 'host_bmi')

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
table(bacter$Bacteroidota == 0) # 17
table(firmi$Firmicutes == 0) # 3
min(bacter[bacter$Bacteroidota > 0, "Bacteroidota"]) # 4
min(firmi[firmi$Firmicutes > 0, "Firmicutes"]) # 105


# Sanity check (should all have 2,576 samples)
# nrow(bacter)
# nrow(firmi)
# nrow(actino)
# nrow(proteo)


# COMPUTE LOG RATIOS
common.columns <- c("Sample", "host_disease", "host_subtype", "sample_type", "Collection",
                    "author", "sequencing_tech", "bowel_movement_quality", "bowel_movement_frequency", "Bristol",
                    "host_age", "host_bmi")

ratio.df <- left_join(x=bacter, y=firmi, by=common.columns) %>%
  relocate(Bacteroidota, .before=Firmicutes) %>%
  # Add 0.5 pseudocounts for the few 0 values
  mutate(Bacteroidota=replace(Bacteroidota, Bacteroidota==0, 0.5),
         Firmicutes=replace(Firmicutes, Firmicutes==0, 0.5)) %>%
  # Compute log ratios
  mutate(LogRatio_FirmBact = log2(Firmicutes/Bacteroidota)) %>%
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

seqtech.order <- c("454 pyrosequencing", "Illumina single-end", "Illumina paired-end", "Ion Torrent")
ratio.df$sequencing_tech <- factor(ratio.df$sequencing_tech, levels=seqtech.order)




####################################################
# PLOT FIRMICUTES/BACTEROIDOTA LOG RATIO FOR PAPER #
####################################################

# Plot Firm/Bact ratio for fecal samples
a <- ggplot(ratio.df %>% filter(Collection=="1st" & sample_type=="stool"),
       aes(x = author, y = LogRatio_FirmBact, fill=host_disease))+
  # facet_wrap(~author, nrow=1, strip.position="bottom")+
  geom_boxplot(position=position_dodge(width=0.75), outlier.shape = NA, width = 0.4, lwd=0.5, alpha=0.2)+
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), aes(color=host_disease), size=0.2)+
  scale_color_manual(values=c("#3182bd", "#de2d26"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#de2d26"))+
  theme_cowplot()+
  # scale_x_discrete(breaks=author_disease.order, labels= rep(c("Healthy", "IBS"), times = 12))+
  theme(axis.text.x = element_text(angle = 45, color="black", hjust=1))+
  labs(x = '', y = expression(Log[2](Firmicutes/Bacteroidota)), fill="", title="Stool samples")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/plots_paper/firm_bact_01.jpg", width=10, height=5)


# Plot Firm/Bact ratio for sigmoid samples
b <- ggplot(ratio.df %>% filter(Collection=="1st" & sample_type=="sigmoid"),
       aes(x = author, y = LogRatio_FirmBact, fill=host_disease))+
  # facet_wrap(~author, nrow=1, strip.position="bottom")+
  geom_boxplot(position=position_dodge(width=0.75), outlier.shape = NA, width = 0.4, lwd=0.5, alpha=0.2)+
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), aes(color=host_disease), size=0.2)+
  scale_color_manual(values=c("#3182bd", "#de2d26"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#de2d26"))+
  theme_cowplot()+
  # scale_x_discrete(breaks=author_disease.order, labels= rep(c("Healthy", "IBS"), times = 12))+
  theme(axis.text.x = element_text(angle = 45, color="black", hjust=1))+
  labs(x = '', y =  expression(Log[2](Firmicutes/Bacteroidota)), fill="", title="Sigmoid samples")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/plots_paper/firm_bact_02.jpg", width=5, height=5)


# Plot Firm/Bact ratio per collection time point
c <- ggplot(data = ratio.df %>%
                group_by(author) %>%
                filter(n_distinct(Collection)>1) %>%
                ungroup() %>%
                mutate(Collection=replace(Collection, Collection=="1st", "t0"),
                       Collection=replace(Collection, Collection=="2nd" & author=="Mars", "t0+6mo"),
                       Collection=replace(Collection, Collection=="2nd" & author=="Pozuelo", "t0+1mo")) %>%
                mutate(author=recode(author, "Mars"="Mars (sigmoid)"),
                       author=recode(author, "Pozuelo"="Pozuelo (stool)")),
       aes(x = Collection, y = LogRatio_FirmBact, fill=host_disease))+
  facet_wrap(~author, nrow=1, strip.position="top", scales="free_x")+
  geom_boxplot(position=position_dodge(width=0.75), outlier.shape = NA, width = 0.4, lwd=0.5, alpha=0.2)+
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), aes(color=host_disease), size=0.2)+
  scale_color_manual(values=c("#3182bd", "#de2d26"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#de2d26"))+
  theme_cowplot()+
  # scale_x_discrete(breaks=author_disease.order, labels= rep(c("Healthy", "IBS"), times = 12))+
  theme(axis.text.x = element_text(angle = 45, color="black", hjust=1))+
  labs(x = '', y = expression(Log[2](Firmicutes/Bacteroidota)), fill="", title="Collection time point")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/plots_paper/firm_bact_03.jpg", width=6, height=5)


# Plot Firm/Bact ratio per IBS subtype
d <- ggplot(
    #+++++++++
    # Keep only data of interest
    data = ratio.df %>% filter(Collection=="1st",
                               host_subtype!="IBS-unspecified", host_subtype!="HC-unknown",
                               author %in% c("Labus", "LoPresti", "Liu", "Pozuelo", "Mars", "Zhuang", "Nagel", "Zeber-Lubecka"),
                               !(author=="LoPresti" & sample_type=="sigmoid") # remove the (very) few sigmoid samples of LoPresti
                               ) %>%
                        mutate(author=recode(author, "Mars"="Mars (sigmoid)")) %>%
                        mutate(host_subtype=replace(host_subtype, host_subtype=="HC", "Healthy")),
    aes(x = host_subtype, y = LogRatio_FirmBact, fill=host_subtype))+
    #+++++++++
  facet_wrap(~author, scales="fixed", ncol=1, strip.position="right")+
  # facet_grid(author~sample_type)+
  #ylim(c(-3,10))+
  # geom_violin(lwd=0.5, alpha=0.2)+
  geom_boxplot(lwd=0.5, width=0.4, outlier.shape=NA, alpha=0.2)+
  geom_jitter(aes(color=host_subtype), width=0.05, size=0.5)+
  scale_color_manual(values=c("#3182bd", "#a50f15", "#fcae91", "#fb6a4a"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#a50f15", "#fcae91", "#fb6a4a"), guide="none")+
  # THEME
  theme_cowplot()+
  theme(strip.text.y = element_text(angle=0, hjust=0),
        strip.background = element_blank())+
  labs(title="IBS subtypes", x = '', y =  expression(Log[2](Firmicutes/Bacteroidota)))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/plots_paper/firm_bact_04.jpg", width=5, height=9)


# Combine plots
ggdraw() +
  draw_plot(a, x = 0, y = .5, width = .65, height = .5) +
  draw_plot(b, x = 0, y = 0, width = .30, height = .5) +
  draw_plot(c, x = 0.30, y = 0, width = .35, height = .5) +
  draw_plot(d, x = 0.65, y = 0, width = 0.35, height = 1) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0, 0.30, 0.65), y = c(1, 0.5, 0.5, 1))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/plots_paper/firm_bact_05.jpg", width=15, height=10)



