##########################
# Purpose: UMAP plotting at all taxonomic levels
# Date: August 2021
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(ggplot2)
library(umap)
library(tidyverse)
library(reshape2)
library(gtools)

# Data
path.phy <- "~/IBS/PhyloTree/input"
physeq.ringel <- readRDS(file.path(path.phy, "physeq_ringel.rds"))
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
cat("\n++ MERGE PHYLOSEQ OBJECTS ++\n")
physeq <- merge_phyloseq(physeq.ringel,
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

# Separate fecal & sigmoid samples
physeq.fecal <- subset_samples(physeq, sample_type == 'stool') # 2,245 samples
physeq.sigmoid <- subset_samples(physeq, sample_type == 'sigmoid') # 431 samples
cat("Nb of fecal samples:", nsamples(physeq.fecal))
cat("\nNb of sigmoid samples:", nsamples(physeq.sigmoid))

# Have pseudocounts
physeq_fecal.pseudocts <- physeq.fecal
otu_table(physeq_fecal.pseudocts)[otu_table(physeq_fecal.pseudocts) == 0] <- 0.5

physeq_sigmoid.pseudocts <- physeq.sigmoid
otu_table(physeq_sigmoid.pseudocts)[otu_table(physeq_sigmoid.pseudocts) == 0] <- 0.5



###############
# AGGLOMERATE #
###############

# Function to obtain table with log-ratios between taxa
LogRatios <- function(abundanceTable){
  
  # Get combinations between all taxa, 1st column numerator, 2nd column denominator, 3rd column taxa1/taxa2
  comb <- as.data.frame(combinations(nrow(abundanceTable), 2, rownames(abundanceTable), repeats.allowed = FALSE))
  comb[,3] <- paste0(comb[,1], "/", comb[,2])
  
  cat("Number of combinations:", nrow(comb), "\n")
  
  # Compute the ratios
  ratios <- as.data.frame(abundanceTable) # taxa as rows, samples as columns
  test <- ratios[1:3,1:3] # for later sanity check
  
  ratios[ comb[,3] ,] <- mapply(function(x, y) log2(x/y), # compute the log ratio
                                abundanceTable[ as.character(comb[,1]) ,], # between first column
                                abundanceTable[ as.character(comb[,2]) ,]) # and second column of all combinations
  
  ratios <- ratios[-c(1:nrow(abundanceTable)),] # remove first x rows that are not ratios (it is simply the taxa)
  
  # sanity check
  cat("Sanity check: log-ratio is", log2(test[1,1]/test[2,1]),
      "and the corresponding value in the log-ratio table is", ratios[1,1], "\n")
  cat("Sanity check2: log-ratio is", log2(test[1,2]/test[3,2]),
      "and the corresponding value in the log-ratio table is", ratios[2,2], "\n")
  
  cat("We have", dim(ratios)[2], "samples and", dim(ratios)[1], "predictors (log-ratios) \n")
  
  # Return table with samples as rows and logratios as columns
  return(t(ratios))
}


# Get the familyTable
fam.agg <- physeq.pseudocts %>%
  tax_glom(taxrank = "Family") %>%
  psmelt()
familyTable <- acast(fam.agg, Family ~ Sample, value.var = 'Abundance')
# saveRDS(fam.agg, "~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/02_UMAP/family_agg.rds")

# sanity checks
# dim(familyTable)
# table(familyTable == 0)

# Get log-ratios for family table
ratiosFam <- LogRatios(familyTable)
# saveRDS(ratiosFam, "~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/02_UMAP/ratiosFam.rds")

# Mean centering
ratiosFam.scaled <- scale(ratiosFam, center = TRUE, scale = FALSE)
table(is.na(ratiosFam.scaled)) # sanity check










