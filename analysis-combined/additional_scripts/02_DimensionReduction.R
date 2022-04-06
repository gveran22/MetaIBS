##########################
# Purpose: Dimension reduction of all samples
# Date: January 2022
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(ggplot2)
library(cowplot)
#library(umap)
library(tidyverse)
#library(reshape2)
#library(gtools)

# Data
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
physeq.ringel   <- readRDS(file.path(path.phy, "physeq_ringel.rds"))
physeq.labus    <- readRDS(file.path(path.phy, "physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path.phy, "physeq_lopresti.rds"))
physeq.pozuelo  <- readRDS(file.path(path.phy, "physeq_pozuelo.rds"))
physeq.zhuang   <- readRDS(file.path(path.phy, "physeq_zhuang.rds"))
physeq.zhu      <- readRDS(file.path(path.phy, "physeq_zhu.rds"))
physeq.hugerth  <- readRDS(file.path(path.phy, "physeq_hugerth.rds"))
physeq.fukui    <- readRDS(file.path(path.phy, "physeq_fukui.rds"))
physeq.mars     <- readRDS(file.path(path.phy, "physeq_mars.rds"))
physeq.liu      <- readRDS(file.path(path.phy, "physeq_liu.rds"))
physeq.agp      <- readRDS(file.path(path.phy, "physeq_agp.rds"))
physeq.nagel    <- readRDS(file.path(path.phy, "physeq_nagel.rds"))
physeq.zeber    <- readRDS(file.path(path.phy, "physeq_zeber.rds"))




###################
# PREPROCESS DATA #
###################

# Merge phyloseq objects
physeq.all <- merge_phyloseq(physeq.ringel,
                             physeq.labus,
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

# Sanity check
physeq.all # 2,651 samples and 79,917 taxa

# Separate fecal & sigmoid samples
# physeq.fecal <- subset_samples(physeq.all, sample_type == 'stool') # 2,228 samples
# physeq.sigmoid <- subset_samples(physeq, sample_type == 'sigmoid') # 431 samples
# cat("Nb of fecal samples:", nsamples(physeq.fecal))
# cat("\nNb of sigmoid samples:", nsamples(physeq.sigmoid))




######################
# DATA NORMALIZATION #
######################

# table(sample_sums(physeq.all)<500) # sanity check
path <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/02_DimensionReduction"
# physeq.NZcomp <- readRDS(file.path(path, "physeq_all_NZcomp.rds"))
# physeq.CSN <- readRDS(file.path(path, "physeq_all_CSN.rds"))
# physeq.clr <- readRDS(file.path(path, "physeq_all_clr.rds"))


# NON-ZERO COMPOSITIONS
physeq.NZcomp <- physeq.all
otu_table(physeq.NZcomp)[otu_table(physeq.NZcomp) == 0] <- 0.5 # pseudocounts
# otu_table(physeq.all)[1:5,1:5]
# otu_table(physeq.NZcomp)[1:5,1:5]
physeq.NZcomp <- transform_sample_counts(physeq.NZcomp, function(x) x / sum(x) )
# table(rowSums(otu_table(physeq.NZcomp))) # check if there is any row not summing to 1
# saveRDS(physeq.NZcomp, file.path(path, "physeq_all_NZcomp.rds"))


# COMMON-SCALE NORMALIZATION
physeq.CSN <- physeq.all
physeq.CSN <- transform_sample_counts(physeq.CSN, function(x) (x*min(sample_sums(physeq.CSN))) / sum(x) )
# table(rowSums(otu_table(physeq.CSN))) # check that all rows are summing to the same total
# saveRDS(physeq.CSN, file.path(path, "data/analysis-combined/02_DimensionReduction/physeq_all_CSN.rds"))


# CENTERED LOG RATIO (CLR) COUNT
physeq.clr <- microbiome::transform(physeq.all, "clr") # the function adds pseudocounts itself
otu_table(physeq.all)[1:5, 1:5] # should contain absolute counts
otu_table(physeq.clr)[1:5, 1:5] # should all be relative
# saveRDS(physeq.clr, file.path(path, "data/analysis-combined/02_DimensionReduction/physeq_all_clr.rds"))




#######################
# DIMENSION REDUCTION #
#######################

# Compute distances
getDistances <- function(){
  print("Aitchison")
  glom.ait <- phyloseq::distance(physeq.clr, method = 'euclidean') # aitchison
  print("Bray")
  glom.bray <- phyloseq::distance(physeq.CSN, method = "bray") # bray-curtis
  print("Canberra")
  glom.can <- phyloseq::distance(physeq.NZcomp, method = "canberra") # canberra
  dist.list <- list("Ait" = glom.ait, "Bray" = glom.bray, "Canb" = glom.can)
  return(dist.list)
}

# Function to plot distances
plotDistances2D <- function(dlist, ordination="MDS", coloring="host_disease"){
  plist <- NULL
  plist <- vector("list", 3)
  names(plist) <- c("Aitchison", "Bray-Curtis", "Canberra")
  
  print("Aitchison")
  # Aitchison
  set.seed(123)
  iMDS.Ait <- ordinate(physeq=physeq.clr, method=ordination, distance=dlist$Ait)
  plist[[1]] <- plot_ordination(physeq.clr, iMDS.Ait)
  
  print("Bray")
  # Bray-Curtis
  set.seed(123)
  iMDS.Bray <- ordinate(physeq=physeq.CSN, method=ordination, distance=dlist$Bray)
  plist[[2]] <- plot_ordination(physeq.CSN, iMDS.Bray)
  
  print("Canberra")
  # Canberra
  set.seed(123)
  iMDS.Can <- ordinate(physeq=physeq.NZcomp, method=ordination, distance=dlist$Can)
  plist[[3]] <- plot_ordination(physeq.NZcomp, iMDS.Can)
  
  # Creating a dataframe to plot everything
  print("Store data")
  plot.df <- plyr::ldply(plist, function(x) x$data)
  names(plot.df)[1] <- "distance"

  return(plot.df)
}




########
# PLOT #
########

# dist.all <- readRDS(file.path(path, "output_ait-bray-can-distances.rds"))

dist.all <- getDistances()
plot.df <- plotDistances2D(dlist=dist.all, ordination="NMDS")
plot.df <- readRDS("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/02_DimensionReduction/output_plot-df-pcoa.rds")

# Plot
ggplot(plot.df, aes(Axis.1, Axis.2, color=author))+
  geom_point(size=2, alpha=0.5)  +# scale_color_manual(values = c('blue', 'red', 'black'))+
  facet_wrap(distance~., scales='free', nrow=1)+
  theme_cowplot()+
  theme(strip.text.x = element_text(size=20))+
  labs(color="Disease")



# XXXXXXXXXXX TEST XXXXXXXXXXX

# Aitchison distance
glom.ait <- phyloseq::distance(physeq.clr, method = 'euclidean') # aitchison (takes ~1h30 to compute)
# saveRDS(glom.ait, "~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/02_DimensionReduction/dist_aitchison_allSamples.rds")
iMDS.Ait <- ordinate(physeq=physeq.clr, method="NMDS", distance=glom.ait, k=3)
plot.df  <- plot_ordination(physeq.clr, iMDS.Ait)$data


# Plot
ggplot(plot.df, aes(Axis.1, Axis.2, color=author))+
  geom_point(size=2, alpha=0.5)  + #scale_color_manual(values = c('blue', 'red', 'black'))+
  theme_cowplot()+
  theme(strip.text.x = element_text(size=20))+
  labs(color="Disease")

ggplot(plot.df, aes(NMDS1, NMDS2, color=host_disease))+
  geom_point(size=2, alpha=0.5) +
  scale_color_manual(values = c('blue', 'red', 'black'))+
  theme_cowplot()+
  theme(strip.text.x = element_text(size=20))+
  labs(color="Disease")

ggplot(plot.df, aes(NMDS1, NMDS2, color=sequencing_tech))+
  geom_point(size=2, alpha=0.5)  + #scale_color_manual(values = c('blue', 'red', 'black'))+
  theme_bw()+
  theme(strip.text.x = element_text(size=20))+
  labs(color="Disease")
