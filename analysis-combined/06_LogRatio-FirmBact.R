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
library(ggsignif)

# Data
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
physeq.labus    <- readRDS(file.path(path.phy, "physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path.phy, "physeq_lopresti.rds"))
# physeq.ringel <- readRDS(file.path(path.phy, "physeq_ringel.rds"))
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
                         physeq.pozuelo,
                         physeq.zhuang,
                         physeq.zhu,
                         physeq.hugerth,
                         physeq.fukui,
                         physeq.mars,
                         physeq.liu,
                         physeq.agp,
                         physeq.nagel,
                         physeq.zeber) # 2,576 samples

# Agglomerate to phylum level and melt to long format
phylumTable <- physeq %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt()




###############################
# COMPUTE FIRM/BACT LOG RATIO #
###############################

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


# Sanity check (should have 2,576 samples)
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
  mutate(LogRatio_FirmBact = log2(Firmicutes/Bacteroidota))


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
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), aes(color=host_disease), size=0.2)+
  geom_boxplot(position=position_dodge(width=0.75), outlier.shape = NA, width = 0.4, lwd=0.5, alpha=0.2)+
  # add p-values (computed in scripts 03_EDA_[dataset-name])
  geom_signif(y_position = c(8,8,8,8),
              xmin = c(0.8,1.8,4.8,7.8),
              xmax = c(1.2,2.2,5.2,8.2),
              annotation = c("*","*","**","**"), tip_length = 0)+
  # theme and color
  scale_color_manual(values=c("#3182bd", "#de2d26"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#de2d26"))+
  ylim(c(-5,17))+ # note: we're cutting one point from Fukui (at y=-12)
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, color="black", hjust=1))+
  labs(x = '', y = expression(Log[2](Firmicutes/Bacteroidota)), fill="", title="Stool samples")


# Plot Firm/Bact ratio for sigmoid samples
b <- ggplot(ratio.df %>% filter(Collection=="1st" & sample_type=="sigmoid"),
       aes(x = author, y = LogRatio_FirmBact, fill=host_disease))+
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), aes(color=host_disease), size=0.2)+
  geom_boxplot(position=position_dodge(width=0.75), outlier.shape = NA, width = 0.4, lwd=0.5, alpha=0.2)+
  scale_color_manual(values=c("#3182bd", "#de2d26"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#de2d26"))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, color="black", hjust=1))+
  labs(x = '', y =  expression(Log[2](Firmicutes/Bacteroidota)), fill="", title="Sigmoid samples")


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
            aes(x = Collection, y = LogRatio_FirmBact))+
  facet_wrap(~author, nrow=1, strip.position="top", scales="free_x")+
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), aes(color=host_disease), size=0.2)+
  geom_boxplot(aes(fill=host_disease), position=position_dodge(width=0.75), outlier.shape = NA, width = 0.4, lwd=0.5, alpha=0.2)+
  geom_signif(y_position=c(8,8), xmin=c(0.8,1.8), xmax=c(1.2,2.2), annotations = c("**","**"), tip_length = 0)+
  # CAREFUL: will need to remove manually (on adobe illustrator) the ** on the Mars facet
  scale_color_manual(values=c("#3182bd", "#de2d26"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#de2d26"))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, color="black", hjust=1))+
  labs(x = '', y = expression(Log[2](Firmicutes/Bacteroidota)), fill="", title="Collection time point")


# Plot Firm/Bact ratio per IBS subtype
annotation.df <- data.frame(author=c("Labus", "LoPresti", "Pozuelo", "Mars (sigmoid)", "Zeber-Lubecka"),
                             start=c("Healthy", "Healthy", "Healthy", "Healthy", "Healthy"),
                             end=c("IBS-D", "IBS-C", "IBS-D", "IBS-C","IBS-C"),
                             label=c("*", "**", "*", "*", "**"),
                             y=c(8,8,3.6,3.5,4))
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
    aes(x = host_subtype, y = LogRatio_FirmBact))+
    #+++++++++
  geom_signif(data=annotation.df, aes(y_position=y, xmin=start, xmax=end, annotations=label), manual=T, tip_length=0, vjust=0.4)+
  facet_wrap(~factor(author, levels=c("Labus", "LoPresti", "Liu", "Pozuelo", "Mars (sigmoid)", "Zhuang", "Nagel", "Zeber-Lubecka")),
             scales="free_y", ncol=1, strip.position="right")+
  # ylim(c(-3,10))+
  geom_jitter(aes(color=host_subtype), width=0.05, size=0.5)+
  geom_boxplot(aes(fill=host_subtype), lwd=0.5, width=0.4, outlier.shape=NA, alpha=0.2)+
  scale_color_manual(values=c("#3182bd", "#a50f15", "#fcae91", "#fb6a4a"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#a50f15", "#fcae91", "#fb6a4a"), guide="none")+
  # THEME
  theme_cowplot()+
  theme(strip.text.y = element_text(angle=0, hjust=0),
        strip.background = element_blank())+
  labs(title="IBS subtypes", x = '', y =  expression(Log[2](Firmicutes/Bacteroidota)))


# Combine plots
ggdraw() +
  draw_plot(a, x = 0, y = .5, width = .65, height = .5) +
  draw_plot(b, x = 0, y = 0, width = .30, height = .5) +
  draw_plot(c, x = 0.30, y = 0, width = .35, height = .5) +
  draw_plot(d, x = 0.65, y = 0, width = 0.35, height = 1) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0, 0.30, 0.65), y = c(1, 0.5, 0.5, 1))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/plots_paper/firm_bact_05.jpg", width=15, height=10)




#################################################
# PLOT FIRMICUTES/BACTEROIDOTA LOG RATIO IN AGP #
#################################################

ratio.df.AGP <- ratio.df %>%
  filter(author=="AGP") %>%
  # Update the host_subtype column to bowel movement
  mutate(host_subtype=replace(host_subtype,
                              host_disease=="IBS" & bowel_movement_quality=="Constipated",
                              "Constipated")) %>%
  mutate(host_subtype=replace(host_subtype,
                              host_disease=="IBS" & bowel_movement_quality=="Diarrhea",
                              "Diarrhea")) %>%
  mutate(host_subtype=replace(host_subtype,
                              host_disease=="IBS" & bowel_movement_quality=="Unknown",
                              "Unknown")) %>%
  mutate(host_subtype=replace(host_subtype,
                              host_disease=="Healthy" & bowel_movement_quality=="Normal",
                              "Healthy"))

# Plot Firm/Bact ratio per stool morphology
bwlmvt.order <- c("Healthy", "Constipated", "Diarrhea", "Unknown")
ratio.df.AGP$host_subtype <- factor(ratio.df.AGP$host_subtype, levels=bwlmvt.order)

a1 <- ggplot(ratio.df.AGP %>% filter(!is.na(host_subtype)),
       aes(x = host_subtype, y = LogRatio_FirmBact, fill=host_disease))+
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), aes(color=host_disease), size=0.2)+
  geom_boxplot(position=position_dodge(width=0.75), outlier.shape = NA, width = 0.4, lwd=0.5, alpha=0.2)+
  # add p-values
  # geom_signif(comparisons = list(c("Healthy", "Constipated"),
  #                                c("Healthy", "Diarrhea"),
  #                                c("Healthy", "Unknown")),
  #             map_signif_level = TRUE, test="wilcox.test", step_increase=.13, tip_length=0) +
  # theme and color
  scale_color_manual(values=c("#3182bd", "#de2d26"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#de2d26"))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, color="black", hjust=1))+
  labs(x = '', y = expression(Log[2](Firmicutes/Bacteroidota)), fill="", title="Bowel movement quality")


# Plot Firm/Bact ratio per stool frequency
ratio.df.AGP$bowel_movement_frequency <- factor(ratio.df.AGP$bowel_movement_frequency,
                                                levels=c("Less than one", "One", "Two", "Three", "Four", "Five or more", "Unknown"))
b1 <- ggplot(ratio.df.AGP,
       aes(x = bowel_movement_frequency, y = LogRatio_FirmBact, fill=host_disease))+
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), aes(color=host_disease), size=0.2)+
  geom_boxplot(position=position_dodge(width=0.75), outlier.shape = NA, width = 0.4, lwd=0.5, alpha=0.2)+
  # theme and color
  scale_color_manual(values=c("#3182bd", "#de2d26"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#de2d26"))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, color="black", hjust=1))+
  labs(x = '', y = expression(Log[2](Firmicutes/Bacteroidota)), fill="", title="Bowel movement frequency")


# Combine plots
ggdraw() +
  draw_plot(a1, x = 0, y = .5, width = 1, height = .5) +
  draw_plot(b1, x = 0, y = 0, width = 1, height = .5) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0), y = c(1, 0.5))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/plots_paper/SuppFig5.jpg", width=6, height=9)

