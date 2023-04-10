# ******************************************
# Purpose: Plotting seq depth pre/post DADA2
# Date: August 2021
# Author: Salomé Carcy
# ******************************************




# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(dplyr)
library(ggplot2)
library(cowplot)
library(scales)
library(RColorBrewer)


## 1.2. Data ####
path.root <- "~/Projects/MetaIBS" # CHANGE THIS ROOT DIRECTORY ON YOUR CLUSTER
# path.data <- file.path(path.root, "data_local/analysis-individual/DADA2-FILT") # for authors of MetaIBS paper
path.data <- file.path("data/analysis-combined/06_QCplot")

datasets     <- list.files(path.data)
print(datasets) # should all be named "nbreads_NameDataset.rds"
nbreads_list <- sapply(datasets, function(x) readRDS(file.path(path.data, x)), USE.NAMES=T, simplify=F)

# rbind all the dataframes
df <- bind_rows(nbreads_list)
df <- reshape2::melt(df, id.vars="dataset", variable.name="QC", value.name="seqdepth")

# Sanity checks
dim(df[df$QC=="before",]) # 2948 samples
# Note: in the flowchart, we say that we processed 2951 raw FASTQ files
# but in the Hugerth dataset, 3 FASTQ files had a weird header and couldn't be processed => 2948 samples




# ************
# 2. PLOT ####
# ************

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
ggsave(file.path(path.data, "seqdepth.jpg"), width=10, height=5)

