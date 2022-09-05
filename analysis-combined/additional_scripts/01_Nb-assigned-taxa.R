##########################
# Purpose: Nb ASVs vs nb taxa with assigned taxonomy
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


# Phyloseq objects
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
physeq.labus    <- readRDS(file.path(path.phy, "physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path.phy, "physeq_lopresti.rds"))
physeq.ringel   <- readRDS(file.path(path.phy, "physeq_ringel.rds"))
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
                             physeq.ringel,
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
cat("Nb of samples:", nsamples(physeq.all), "\n") # should be 2651
cat("Nb of taxa:", ntaxa(physeq.all), "\n") # should be 79,943




###########################################################
# GET PROPORTION OF ASVs ASSIGNED TO EACH TAXONOMIC LEVEL #
###########################################################

taxtable <- as.data.frame(tax_table(physeq.all))
nrow(taxtable) # 79,943 taxa


df <- apply(taxtable, 2, function(x) sum(!is.na(x))*100/nrow(taxtable))
df


