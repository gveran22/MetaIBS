##########################
# Purpose: Correlation seqdepth & alpha-diversity
# Date: March 2022
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


# Data (phyloseq object for shannon index)
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
physeq.labus    <- readRDS(file.path(path.phy, "physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path.phy, "physeq_lopresti.rds"))
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


# Data (seqdepth)
path <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/DADA2-FILT"
seqdepth.df <- rbind(readRDS(file.path(path, "seqdepth_labus.rds")),
                     readRDS(file.path(path, "seqdepth_lopresti.rds")),
                     readRDS(file.path(path, "seqdepth_agp.rds")),
                     readRDS(file.path(path, "seqdepth_liu.rds")),
                     readRDS(file.path(path, "seqdepth_pozuelo.rds")),
                     readRDS(file.path(path, "seqdepth_fukui.rds")),
                     readRDS(file.path(path, "seqdepth_mars.rds")),
                     readRDS(file.path(path, "seqdepth_hugerth.rds")),
                     readRDS(file.path(path, "seqdepth_zhu.rds")),      
                     readRDS(file.path(path, "seqdepth_zhuang.rds")),
                     readRDS(file.path(path, "seqdepth_nagel.rds")),
                     readRDS(file.path(path, "seqdepth_zeber.rds")))
rownames(seqdepth.df) <- gsub(".fastq.*", "", rownames(seqdepth.df))
rownames(seqdepth.df) <- gsub("_F.*", "", rownames(seqdepth.df))




###############
# CORRELATION #
###############

# Get Shannon values
plt.shannon <- plot_richness(physeq.all, measures="Shannon")
shannon.df <- plt.shannon$data %>%
  select(c("samples", "value", "host_disease", "host_subtype", "sample_type", "Collection", "author")) %>%
  dplyr::rename(shannon=value)


# as some samples got deleted post-QC (samples < 500 count), we will remove them from the seqdepth.df
seqdepth.df <- seqdepth.df[rownames(seqdepth.df) %in% shannon.df$samples,-1]
seqdepth.df$samples <- rownames(seqdepth.df)

# Get dataframe with sample + seqdepth + shannon
corr.df <- merge(seqdepth.df, shannon.df, by="samples")
head(corr.df)
corr.df <- corr.df[,-c(8,9)]
colnames(corr.df)[2] <- "seqdepth"

# Plot
ggplot(corr.df, aes(x=seqdepth, y=shannon, color=host_disease))+
  geom_point(size=.2)+
  scale_x_continuous(trans='log10')+
  scale_color_manual(values=c("blue", "red"))+
  theme_cowplot()
