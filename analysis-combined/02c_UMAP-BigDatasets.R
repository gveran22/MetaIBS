#
# Purpose: PCoA on Bray-Curtis in big datasets
# Date: July 2022
# Author: Salomé Carcy
#




# 1. IMPORT ----------------------------------------------------------------------------------------

# Libraries
library(phyloseq)
library(ggplot2)
library(cowplot)
library(umap)
library(tidyverse)
library(reshape2)
library(gtools)


# Phyloseq objects
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
physeq.agp      <- readRDS(file.path(path.phy, "physeq_agp.rds"))
physeq.pozuelo  <- readRDS(file.path(path.phy, "physeq_pozuelo.rds"))
physeq.hugerth  <- readRDS(file.path(path.phy, "physeq_hugerth.rds"))




# 2. Bray-Curtis on PCoA -------------------------------------------------------------------------------

# 2.1. Functions ===================================================

getBray <- function(physeq, tax_rank){
  
  # If want to agglomerate to a certain taxonomic level, do it before normalization
  if(tax_rank != "none"){
    cat("\n++ Agglomerate to", tax_rank, "level ++\n")
    physeq <- tax_glom(physeq, taxrank=tax_rank)
  }
  
  # Normalize
  cat("++ Normalize ++\n")
  physeq.CSN <- transform_sample_counts(physeq, function(x) (x*min(sample_sums(physeq))) / sum(x) )
  print(table(rowSums(otu_table(physeq.CSN)))) # check that all rows are summing to the same total
  
  # Compute BrayCurtis
  cat("\n++ Compute Bray-Curtis distance ++\n")
  bray.dist <- phyloseq::distance(physeq.CSN, method = "bray")
  
  # Perform PCoA
  cat("++ Perform PCoA ++\n")
  set.seed(123)
  iMDS.Bray <- ordinate(physeq=physeq.CSN, method="PCoA", distance=bray.dist)
  plt.data <- plot_ordination(physeq.CSN, iMDS.Bray)$data
  
  return(plt.data)
}



# 2.2. Compute Bray-Curtis ==========================================

bray.agp <- getBray(physeq=physeq.agp, tax_rank="none")
bray.pozuelo <- getBray(physeq=physeq.pozuelo, tax_rank="none")
# bray.pozuelo <- getBray(physeq=subset_samples(physeq.pozuelo, Collection=="1st"), tax_rank="none")
bray.hugerth <- getBray(physeq=physeq.hugerth, tax_rank="none")
# bray.hugerth <- getBray(physeq=subset_samples(physeq.hugerth, sample_type=="stool"), tax_rank="none")



# 2.3. Plot PCoA ==================================================

# 2.3.1. AGP ####
d <- ggplot(bray.agp, aes(x = Axis.1, y = Axis.2, color = host_disease))+
  geom_point(size = 1.5)+
  scale_color_manual(values = c('blue', 'red', 'black'), name="")+ # disease
  labs(title = "AGP")+
  theme_cowplot()+
  theme(line = element_blank(),
        axis.text = element_text(size=8),
        legend.position = c(0.8,1),
        plot.title = element_text(hjust = 0.5))


# 2.3.2. Pozuelo ####
e1 <- ggplot(bray.pozuelo, aes(x = Axis.1, y = Axis.2, color = host_disease))+
  geom_line(aes(group=host_ID), color="black", lwd=0.1) +
  geom_point(size = 1)+
  scale_color_manual(values = c('blue', 'red', 'black'), name="")+ # disease
  labs(title = "Pozuelo")+
  theme_cowplot()+
  theme(line = element_blank(),
        axis.text = element_text(size=8),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.1,0.95),
        plot.title = element_text(hjust = 0.5))

e2 <- ggplot(bray.pozuelo, aes(x = Axis.1, y = Axis.2, color = Collection))+
  geom_line(aes(group=host_ID), color="black", lwd=0.1) +
  geom_point(size = 1)+
  scale_color_manual(values = c('#f1a340', '#998ec3'), name="Time point")+
  theme_cowplot()+
  theme(line = element_blank(),
        axis.text = element_text(size=8),
        legend.position = c(0.1,0.95),
        plot.title = element_text(hjust = 0.5))

# 2.3.3. Hugerth ####
f1 <- ggplot(bray.hugerth, aes(x = Axis.1, y = Axis.2, color = host_disease))+
  geom_line(aes(group=host_ID), color="black", lwd=0.1) +
  geom_point(size = 1)+
  scale_color_manual(values = c('blue', 'red', 'black'), name="")+ # disease
  labs(title = "Hugerth")+
  theme_cowplot()+
  theme(line = element_blank(),
        axis.text = element_text(size=8),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.1,0.95),
        plot.title = element_text(hjust = 0.5))

f2 <- ggplot(bray.hugerth, aes(x = Axis.1, y = Axis.2, color = sample_type))+
  geom_line(aes(group=host_ID), color="black", lwd=0.1) +
  geom_point(size = 1)+
  scale_color_manual(values = c('#7fbf7b', '#af8dc3'), name="Sample type")+
  theme_cowplot()+
  theme(line = element_blank(),
        axis.text = element_text(size=8),
        legend.position = c(0.1,0.95),
        plot.title = element_text(hjust = 0.5))


# 2.3.4. Combine plots ####
ggdraw() +
  draw_plot(d,  x = 0,    y = .1, width = .39, height = .8) +
  draw_plot(e1, x = .4,   y = .5, width = .29, height = .5) +
  draw_plot(e2, x = .4,   y = 0,  width = .29, height = .5) +
  draw_plot(f1, x = .7,   y = .5, width = .29, height = .5) +
  draw_plot(f2, x = .7,   y = 0,  width = .29, height = .5)
ggsave("~/Projects/IBS_Meta-analysis_16S/data/plots_paper/AGP-Poz-Hug_bray-PCoA.jpeg", width=15, height=5)





