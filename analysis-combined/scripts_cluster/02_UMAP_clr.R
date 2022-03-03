##########################
# Purpose: UMAP plotting at Genus level
# Date: August 2021
# Author: Salomé Carcy
##########################




##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(microbiome)
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

# Separate fecal & sigmoid samples
physeq.fecal <- subset_samples(physeq.all, sample_type == 'stool') # 2,228 samples
# physeq.sigmoid <- subset_samples(physeq, sample_type == 'sigmoid') # 431 samples
cat("Nb of fecal samples:", nsamples(physeq.fecal))
# cat("\nNb of sigmoid samples:", nsamples(physeq.sigmoid))

# CLR-TRANSFORM DATA
physeq_fecal.clr <- microbiome::transform(physeq.fecal, "clr") # the function adds pseudocounts itself




#############
# FUNCTIONS #
#############

#________________________________________________________
# Function to agglomerate to given taxonomic level
aggTable <- function(physeq, tax_rank){
  
  cat("\n++ AGGLOMERATE TO", tax_rank, "LEVEL++\n")
  long.agg <- physeq %>%
    tax_glom(taxrank = tax_rank) %>%
    psmelt()
  
  # Get a matrix Samples (rows) x TaxRank (columns)
  if(tax_rank=="Phylum"){agglomeratedTable <- acast(long.agg, Sample ~ Phylum, value.var = 'Abundance')}
  else if(tax_rank=="Class"){agglomeratedTable <- acast(long.agg, Sample ~ Class, value.var = 'Abundance')}
  else if(tax_rank=="Order"){agglomeratedTable <- acast(long.agg, Sample ~ Order, value.var = 'Abundance')}
  else if(tax_rank=="Family"){agglomeratedTable <- acast(long.agg, Sample ~ Family, value.var = 'Abundance')}
  else if(tax_rank=="Genus"){agglomeratedTable <- acast(long.agg, Sample ~ Genus, value.var = 'Abundance')}
  cat("-> Dimensions matrix:", dim(agglomeratedTable), "\n")
  cat("-> Any 0 counts?", table(agglomeratedTable == 0), "\n")
  
  # Save
  filepath <- paste0("~/IBS/UMAP/data_clr/", paste0(lowercase(tax_rank), "_agg.rds", sep=""), sep="")
  saveRDS(object=long.agg, file=filepath)
  
  return(agglomeratedTable)
}


#________________________________________________________
# Function to get the table to feed into the UMAP
getTable <- function(physeq, tax_rank){
  
  # Get the agglomerated taxa
  agglomeratedTable <- aggTable(physeq, tax_rank)
  
  return(agglomeratedTable)
}


#________________________________________________________
# Function to run the UMAP
runUMAP <- function(physeq, tax_rank){
  
  # Get the data table
  agglomeratedTable <- getTable(physeq, tax_rank)
  
  cat("\n++RUN UMAP...++\n")
  # Run UMAP
  set.seed(123)
  umap <- uwot::umap(agglomeratedTable, # umap on samples (rows) and taxa abundances CLR transformed (columns)
                     n_neighbors=50, n_components=3, n_threads=16)
  
  # Save
  filepath <- paste0("~/IBS/UMAP/data_clr/umap", paste0(tax_rank, ".rds", sep=""), sep="")
  saveRDS(object=umap, file=filepath)
  
  # Get the (x,y) coordinates from the UMAP
  dims.umap <- umap %>% as.data.frame()
  colnames(dims.umap) <- c("UMAP_1", "UMAP_2", "UMAP_3")
  rownames(dims.umap) <- rownames(agglomeratedTable)
  
  # Add covariates
  covariates <- data.frame(disease = sample_data(physeq)[,'host_disease'],
                           subtype = sample_data(physeq)[,'host_subtype'],
                           seq_tech = sample_data(physeq)[,'sequencing_tech'],
                           author = sample_data(physeq)[,'author'],
                           variable_region = sample_data(physeq)[,'variable_region'],
                           age = sample_data(physeq)[,'host_age'],
                           bmi = sample_data(physeq)[,'host_bmi'])
  dims.umap <- merge(as.data.frame(dims.umap), covariates, by="row.names") # umap+covariates
  
  # Save
  filepath <- paste0("~/IBS/UMAP/data_clr/dims_umap", paste0(tax_rank, ".rds", sep=""), sep="")
  saveRDS(object=dims.umap, file=filepath)
  
  return(dims.umap)
}


#________________________________________________________
# Function to plot the UMAP
plotUMAP <- function(physeq, tax_rank){
  
  # Get the data table
  dims.umap <- runUMAP(physeq, tax_rank)
  
  cat("\n++PLOT UMAP...++\n")
  
  # PLOT host_disease
  ggplot(dims.umap, aes(x = UMAP_1, y = UMAP_2, color = host_disease))+
    geom_point(size = 2, alpha = 0.7)+
    scale_color_manual(values = c('blue', 'red', 'black'))+
    labs(title = paste0('Agglomeration at taxonomic level:', tax_rank))+
    theme_bw()
  ggsave(paste0("~/IBS/UMAP/data_clr/umap", paste0(tax_rank, "_disease.jpg", sep=""), sep=""), width=6, height=4, type="cairo")
  
  # PLOT host_subtype
  ggplot(dims.umap, aes(x = UMAP_1, y = UMAP_2, color = host_subtype))+
    geom_point(size = 1, alpha = 0.7)+
    scale_color_manual(values = c('#99CCFF', '#FF3300', '#990000', '#FF66CC','#FFFF66', '#CCCCCC'))+
    labs(title = paste0('Agglomeration at taxonomic level:', tax_rank))+
    theme_bw()
  ggsave(paste0("~/IBS/UMAP/data_clr/umap", paste0(tax_rank, "_subtype.jpg", sep=""), sep=""), width=6, height=4, type="cairo")
  
  # PLOT author
  ggplot(dims.umap, aes(x = UMAP_1, y = UMAP_2, color = author))+
    geom_point(size = 2, alpha = 0.7)+
    labs(title = paste0('Agglomeration at taxonomic level:', tax_rank))+
    theme_bw()
  ggsave(paste0("~/IBS/UMAP/data_clr/umap", paste0(tax_rank, "_author.jpg", sep=""), sep=""), width=6, height=4, type="cairo")
  
  # PLOT seqtech
  ggplot(dims.umap, aes(x = UMAP_1, y = UMAP_2, color = sequencing_tech))+
    geom_point(size = 2, alpha = 0.7)+
    scale_color_manual(values = c('#6600FF', '#33CC33', '#006600', '#FF6633'))+
    labs(title = paste0('Agglomeration at taxonomic level:', tax_rank))+
    theme_bw()
  ggsave(paste0("~/IBS/UMAP/data_clr/umap", paste0(tax_rank, "_seqtech.jpg", sep=""), sep=""), width=6, height=4, type="cairo")
  # plotly::plot_ly(dims.umapFam_fecal_clr, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, color=~sequencing_tech, type="scatter3d", mode="markers")
  
  # PLOT variable region
  ggplot(dims.umap, aes(x = UMAP_1, y = UMAP_2, color = variable_region))+
    geom_point(size = 2, alpha = 0.7)+
    labs(title = paste0('Agglomeration at taxonomic level:', tax_rank))+
    theme_bw()
  ggsave(paste0("~/IBS/UMAP/data_clr/umap", paste0(tax_rank, "_vregion.jpg", sep=""), sep=""), width=6, height=4, type="cairo")
}




############
# RUN UMAP #
############

# taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")
taxranks <- c("Genus")

for(taxa in taxranks){
  plotUMAP(physeq=physeq_fecal.clr, tax_rank = taxa)
}





