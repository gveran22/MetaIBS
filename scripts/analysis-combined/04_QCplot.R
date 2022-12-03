##########################
# Purpose: Plotting seq depth pre/post DADA2
# Date: August 2021
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(ggplot2)
library(cowplot)
library(scales)
library(RColorBrewer)

# Data
path <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/DADA2-FILT"
df <- reshape2::melt(rbind(readRDS(file.path(path, "seqdepth_labus.rds")),
                              readRDS(file.path(path, "seqdepth_lopresti.rds")),
                              readRDS(file.path(path, "seqdepth_ringel.rds")),
                              readRDS(file.path(path, "seqdepth_agp.rds")),
                              readRDS(file.path(path, "seqdepth_liu.rds")),
                              readRDS(file.path(path, "seqdepth_pozuelo.rds")),
                              readRDS(file.path(path, "seqdepth_fukui.rds")),
                              readRDS(file.path(path, "seqdepth_mars.rds")),
                              readRDS(file.path(path, "seqdepth_hugerth.rds")),
                              readRDS(file.path(path, "seqdepth_zhu.rds")),      
                              readRDS(file.path(path, "seqdepth_zhuang.rds")),
                              readRDS(file.path(path, "seqdepth_nagel.rds")),
                              readRDS(file.path(path, "seqdepth_zeber.rds"))),
                       id.vars="dataset", variable.name="QC", value.name="seqdepth")

# Sanity checks
dim(df[df$QC=="before",]) # 2948 samples
# Note: in the flowchart, we say that we processed 2951 raw FASTQ files
# but in the Hugerth dataset, 3 FASTQ files had a weird header and couldn't be processed => 2948 samples


########
# PLOT #
########

# Set order
authors_order <- c('Labus', 'LoPresti', 'Ringel',
                   'AGP', 'Liu', 'Pozuelo',
                   'Fukui', 'Hugerth', 'Mars','Zhu', 'Zhuang',
                   'Nagel', 'Zeber')
df$dataset <- factor(df$dataset, levels=authors_order)

# Plot
ggplot(df, aes(x = dataset, y = seqdepth, fill=QC))+
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), aes(color=QC), size=0.2)+
  geom_boxplot(position=position_dodge(width=0.75), outlier.shape = NA, width = 0.4, lwd=0.5, alpha=0.3)+
  theme_cowplot()+
  scale_y_continuous(labels=comma, limits=c(0,2e5))+
  scale_color_manual(values=c("#bf812d", "#35978f"), guide="none")+
  scale_fill_manual(values=c("#dfc27d", "#35978f"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x = '', y = "# reads per sample")

# Save figure
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/04_QCplot/seqdepth.jpg", width=10, height=5)

