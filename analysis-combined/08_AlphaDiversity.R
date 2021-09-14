##########################
# Purpose: Alpha-diversity
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
library(ggpubr)

# Data
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input/"
physeq.labus <- readRDS(file.path(path.phy, "physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path.phy, "physeq_lopresti.rds"))
physeq.agp <- readRDS(file.path(path.phy, "physeq_agp.rds"))
physeq.liu <- readRDS(file.path(path.phy, "physeq_liu.rds"))
physeq.pozuelo <- readRDS(file.path(path.phy, "physeq_pozuelo.rds"))
physeq.fukui <- readRDS(file.path(path.phy, "physeq_fukui.rds"))
physeq.hugerth <- readRDS(file.path(path.phy, "physeq_hugerth.rds"))
physeq.mars <- readRDS(file.path(path.phy, "physeq_mars.rds"))
physeq.zhu <- readRDS(file.path(path.phy, "physeq_zhu.rds"))
physeq.zhuang <- readRDS(file.path(path.phy, "physeq_zhuang.rds"))
physeq.nagel <- readRDS(file.path(path.phy, "physeq_nagel.rds"))
physeq.zeber <- readRDS(file.path(path.phy, "physeq_zeber.rds"))

# Merge
physeq.all <- merge_phyloseq(physeq.labus,
                             physeq.lopresti,
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
author.order <- c('Labus', 'LoPresti', # 454 pyrosequencing
                  'AGP', 'Liu', 'Pozuelo', # Illumina single end
                  'Fukui', 'Hugerth', 'Mars', 'Zhu', 'Zhuang', # Illumina paired end
                  'Nagel', 'Zeber-Lubecka') # Ion Torrent
sample_data(physeq.all)$author <- factor(sample_data(physeq.all)$author, levels=author.order)

# physeq_list <- list("Labus"    = physeq.labus,
#                     "LoPresti" = physeq.lopresti,
#                     "AGP"      = physeq.agp,
#                     "Liu"      = physeq.liu,
#                     "Pozuelo"  = physeq.pozuelo,
#                     "Fukui"    = physeq.fukui,
#                     "Hugerth"  = physeq.hugerth,
#                     "Zhu"      = physeq.zhu,
#                     "Zhuang"   = physeq.zhuang,
#                     "Nagel"    = physeq.nagel,
#                     "Zeber"    = physeq.zeber)

# Infer IBS subtype in AGP data
sample_data(physeq.agp)[sample_data(physeq.agp)$host_disease=="IBS" & sample_data(physeq.agp)$bowel_movement_quality=="Constipated", "host_subtype"] <- "IBS-C"
sample_data(physeq.agp)[sample_data(physeq.agp)$host_disease=="IBS" & sample_data(physeq.agp)$bowel_movement_quality=="Diarrhea", "host_subtype"] <- "IBS-D"
sample_data(physeq.agp)[sample_data(physeq.agp)$host_disease=="Healthy" & sample_data(physeq.agp)$bowel_movement_quality!="Normal", "host_subtype"] <- "HC-unknown"




###########
# SHANNON #
###########

# plotDiversity <- function(physeq.list, xAxis="host_disease", alphaIndex="Shannon"){
#   
#   # List of plots
#   plist <- NULL
#   plist <- vector("list", length(physeq.list))
#   names(plist) <- names(physeq.list)
#   
#   # Iterate dataset per dataset
#   for(author in names(physeq.list)){
#     print(author)
#     physeq <- physeq.list[[author]]
#     physeq <- subset_samples(physeq, sample_type=="stool" & Collection=="1st")
#     plist[[author]] <- plot_richness(physeq, x=xAxis, measures=alphaIndex) +
#       geom_boxplot(fill=NA) +
#       theme_bw() +
#       labs(title=author, x="")
#   }
#   
#   return(plist)
# }
# 
# # Plot each dataset separately
# chao1.list <- plotDiversity(physeq_list)
# ggarrange(plotlist=chao1.list)


# Plot merged phyloseq object (host_disease)
plot_richness(subset_samples(physeq.all, sample_type=="stool" & Collection == "1st"),
              x="host_disease", measures="Shannon") +
  geom_boxplot(fill=NA, width=0.3) +
  facet_wrap(~author, scales="fixed") +
  theme_bw() +
  labs(x="", y="Shannon")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/08_AlphaDiversity/shannon_disease.jpg", width=6, height=8)

# Plot merged phyloseq object (host_subtype)
plot_richness(subset_samples(physeq.all, sample_type=="stool" & Collection == "1st"),
              x="host_subtype", measures="Shannon") +
  geom_boxplot(fill=NA, width=0.5) +
  facet_wrap(~author, scales="fixed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  labs(x="", y="Shannon")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/08_AlphaDiversity/shannon_subtype.jpg", width=6, height=8)

# Plot sigmoid samples
plot_richness(subset_samples(physeq.all, sample_type=="sigmoid" & Collection == "1st" & author != "LoPresti"),
              x="host_disease", measures="Shannon") +
  geom_boxplot(fill=NA, width=0.3) +
  facet_wrap(~author, scales="fixed") +
  theme_bw() +
  labs(x="", y="Shannon", title="Sigmoid samples")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/08_AlphaDiversity/shannon_sigmoid.jpg", width=4, height=4)




###########
# SIMPSON #
###########

# Plot merged phyloseq object (host_disease)
plot_richness(subset_samples(physeq.all, sample_type=="stool" & Collection == "1st"),
              x="host_disease", measures="Simpson") +
  geom_boxplot(fill=NA, width=0.3) +
  facet_wrap(~author, scales="fixed") +
  ylim(c(0.5,1)) +
  theme_bw() +
  labs(x="", y="Simpson")

ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/08_AlphaDiversity/simpson_disease.jpg", width=6, height=8)

# Plot merged phyloseq object (host_subtype)
plot_richness(subset_samples(physeq.all, sample_type=="stool" & Collection == "1st"),
              x="host_subtype", measures="Simpson") +
  geom_boxplot(fill=NA, width=0.5) +
  facet_wrap(~author, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  labs(x="", y="Simpson")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/08_AlphaDiversity/simpson_subtype.jpg", width=6, height=8)

# Plot sigmoid samples
plot_richness(subset_samples(physeq.all, sample_type=="sigmoid" & Collection == "1st" & author != "LoPresti"),
              x="host_disease", measures="Simpson") +
  geom_boxplot(fill=NA, width=0.3) +
  facet_wrap(~author, scales="fixed") +
  theme_bw() +
  labs(x="", y="Simpson", title="Sigmoid samples")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/08_AlphaDiversity/simpson_sigmoid.jpg", width=4, height=4)





###############
# AGP SUBTYPE #
###############


plot_richness(subset_samples(physeq.agp, host_subtype != "HC-unknown" & host_subtype != "IBS-unspecified"),
              x="host_subtype", measures=c("Shannon", "Simpson")) +
  geom_boxplot(fill=NA, width=0.3) +
  theme_bw() +
  labs(x="", y="", title="AGP")
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/08_AlphaDiversity/agp_subtype.jpg", width=4, height=4)


