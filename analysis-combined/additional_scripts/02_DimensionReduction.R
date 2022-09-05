##########################
# Purpose: Compute Aitchison, Bray-Curtis, Canberra distances on all samples
# Date: July 2022
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
physeq.all # 2,651 samples and 79,943 taxa


# Separate fecal & sigmoid samples
# physeq.fecal <- subset_samples(physeq.all, sample_type == 'stool') # 2,220 samples
# physeq.fecal <- prune_taxa(taxa_sums(physeq.fecal)>0, physeq.fecal) # remove ASVs that are not present anymore
# cat("Nb of fecal samples:", nsamples(physeq.fecal))
# 
# physeq.sigmoid <- subset_samples(physeq.all, sample_type == 'sigmoid') # 431 samples
# physeq.sigmoid <- prune_taxa(taxa_sums(physeq.sigmoid)>0, physeq.sigmoid) # remove ASVs that are not present anymore
# cat("\nNb of sigmoid samples:", nsamples(physeq.sigmoid))




######################
# DATA NORMALIZATION #
######################

# If want to agglomerate to a certain taxonomic level, do it before normalization
# physeq.all <- tax_glom(physeq.all, "Family")

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
# saveRDS(physeq.CSN, file.path(path, "physeq_all_CSN.rds"))


# CENTERED LOG RATIO (CLR) COUNT
physeq.clr <- microbiome::transform(physeq.all, "clr") # the function adds pseudocounts itself
otu_table(physeq.all)[1:5, 1:5] # should contain absolute counts
otu_table(physeq.clr)[1:5, 1:5] # should all be relative
# saveRDS(physeq.clr, file.path(path, "physeq_all_clr.rds"))




#############
# FUNCTIONS #
#############

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
plotDistances2D <- function(dlist, ordination="MDS"){
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




##################################
# COMPUTE DISTANCES & ORDINATION #
##################################

# Compute distances
dist.all <- getDistances()
# saveRDS(dist.all, file.path(path, "output_ait-bray-can-distances.rds"))
# dist.all <- readRDS(file.path(path, "output_ait-bray-can-distances.rds"))


# Compute ordination (MDS)
plot.df.mds <- plotDistances2D(dlist=dist.all, ordination="MDS")
# saveRDS(plot.df.mds, file.path(path, "output_plot-df-MDS.rds"))
# plot.df.mds <- readRDS(file.path(path, "output_plot-df-MDS.rds"))


# Compute ordination (NMDS)
plot.df.nmds <- plotDistances2D(dlist=dist.all, ordination="NMDS")
# saveRDS(plot.df.nmds, file.path(path, "output_plot-df-NMDS.rds"))
# plot.df.nmds <- readRDS(file.path(path, "output_plot-df-NMDS.rds"))
# CAREFUL !!!! NMDS DID NOT CONVERGE !!!!!


# Set order for authors & seqtech
ordered_author <- c('Labus', 'LoPresti', 'Ringel', # 454 pyrosequencing
                    'AGP', 'Liu', 'Pozuelo', # Illumina single end
                    'Fukui', 'Hugerth', 'Mars', 'Zhu', 'Zhuang', # Illumina paired end
                    'Nagel', 'Zeber-Lubecka')
ordered_seqtech <- c("454 pyrosequencing", "Illumina single-end",
                     "Illumina paired-end", "Ion Torrent")




########
# PLOT #
########


plot <- function(plot.df, dist){
  
  # Set order for authors & seqtech
  plot.df <- plot.df %>%
    mutate(author = factor(plot.df$author, levels = ordered_author)) %>%
    mutate(sequencing_tech = factor(plot.df$sequencing_tech, levels= ordered_seqtech))
  
  # Keep only the values from the distance (Aitchison, Bray or Canberra)
  plot.df <- plot.df %>% filter(distance==dist)
  
  # Per host_disease
  a <- ggplot(plot.df, aes(Axis.1, Axis.2, color=host_disease))+
    geom_point(size=.5) +
    scale_color_manual(values = c('blue', 'red', 'black'), name="")+ # disease
    labs(title = "Disease phenotype")+
    theme_cowplot()+
    theme(line = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  # Per author/dataset
  b <- ggplot(plot.df, aes(Axis.1, Axis.2, color=author))+
    geom_point(size=.5) +
    scale_color_manual(values=pals::brewer.paired(13), name="")+ # author
    labs(title = "Dataset")+
    theme_cowplot()+
    theme(line = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  # Per seq. technology
  c <- ggplot(plot.df, aes(Axis.1, Axis.2, color=sequencing_tech))+
    geom_point(size=.5) +
    scale_color_manual(values=c("#6a51a3", "#a1d99b", "#238b45", "#f16913"), name="")+ # seqtech
    labs(title = "Seq. technology")+
    theme_cowplot()+
    theme(line = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  ggdraw() +
    draw_plot(a, x = 0,     y = 0, width = .32, height = 1) +
    draw_plot(b, x = .33,   y = 0, width = .32, height = 1) +
    draw_plot(c, x = .66,   y = 0, width = .32, height = 1)
  
}


# *******
#   MDS
# *******

plot(plot.df=plot.df.mds, dist="Aitchison")
plot(plot.df=plot.df.mds2, dist="Bray-Curtis")
plot(plot.df=plot.df.mds, dist="Canberra")


# *******
#  NMDS
# *******

# CAREFUL !!!! NMDS DID NOT CONVERGE !!!!!
# plot.df.nmds <- plot.df.nmds %>%
#   rename("Axis.1" = "NMDS1",
#          "Axis.2" = "NMDS2")
# 
# plot(plot.df=plot.df.nmds, dist="Aitchison")
# plot(plot.df=plot.df.nmds, dist="Bray-Curtis")
# plot(plot.df=plot.df.nmds, dist="Canberra")
