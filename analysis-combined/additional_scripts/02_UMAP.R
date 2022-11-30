##########################
# Purpose: UMAP plotting at all taxonomic levels
# Date: March 2021
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




###################
# PREPROCESS DATA #
###################

# Merge phyloseq objects
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

# Separate fecal & sigmoid samples
physeq.fecal <- subset_samples(physeq.all, sample_type == 'stool') # 2,220 samples
physeq.fecal <- prune_taxa(taxa_sums(physeq.fecal)>0, physeq.fecal) # remove ASVs that are not present anymore
print(nsamples(physeq.fecal))

# physeq.sigmoid <- subset_samples(physeq.all, sample_type == 'sigmoid') # 431 samples
# physeq.sigmoid <- prune_taxa(taxa_sums(physeq.sigmoid)>0, physeq.sigmoid)
# nsamples(physeq.sigmoid)

# Have pseudocounts
# physeq_fecal.pseudocts <- physeq.fecal
# otu_table(physeq_fecal.pseudocts)[otu_table(physeq_fecal.pseudocts) == 0] <- 0.5

# physeq_sigmoid.pseudocts <- physeq.sigmoid
# otu_table(physeq_sigmoid.pseudocts)[otu_table(physeq_sigmoid.pseudocts) == 0] <- 0.5




#############
# FUNCTIONS #
#############

#________________________________________________________
# Function to obtain table with log-ratios between taxa
LogRatios <- function(abundanceTable, tax_rank="Family"){
  
  cat("\n++ GET LOG-RATIOS BETWEEN", tax_rank, "++\n")
  
  # Get combinations between all taxa: 1st column numerator, 2nd column denominator, 3rd column taxa1/taxa2
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




###################################
# GET LOG-RATIOS BETWEEN FAMILIES #
###################################

# Agglomerate data to family level
family.agg <- physeq.fecal %>%
  tax_glom(taxrank = "Family") %>%
  psmelt()
# head(family.agg)

# Get a matrix Families (rows) x Samples (columns)
familyTable <- acast(family.agg, Family ~ Sample, fun.aggregate=sum, value.var = 'Abundance')

# Sanity checks
# familyTable[1:10,1:10]
# min(familyTable)
# max(familyTable)
# randomTaxa   <- rownames(familyTable)[7]
# randomSample <- colnames(familyTable)[7]
# valueAgg <- sum(family.agg[family.agg[,"Family"] == randomTaxa & family.agg$Sample == randomSample, "Abundance"])
# wideAgg  <- familyTable[randomTaxa, randomSample]
# cat("-> Sanity check: value in long shaped df is", valueAgg, "while in wide shaped df is", wideAgg)

# Add pseudocounts (before computing log-ratios)
table(familyTable==0)
familyTable[familyTable == 0] <- 0.5

# Compute log-ratios between families
ratiosFamily <- LogRatios(familyTable, "Family")

# Mean-center
# ratios.scaled <- scale(ratiosFamily, center = TRUE, scale = FALSE)




########
# UMAP #
########

# Run UMAP
set.seed(123)
umap <- uwot::umap(ratiosFamily, # umap on samples (rows) and taxa ratios (columns)
                   n_neighbors=50, n_components=3, n_threads=8)

# Get the (x,y) coordinates from the UMAP
dims.umap <- umap %>% as.data.frame()
colnames(dims.umap) <- c("UMAP_1", "UMAP_2", "UMAP_3")
rownames(dims.umap) <- rownames(ratiosFamily)

# Add covariates
covariates <- data.frame(disease = sample_data(physeq)[,'host_disease'],
                         subtype = sample_data(physeq)[,'host_subtype'],
                         seq_tech = sample_data(physeq)[,'sequencing_tech'],
                         author = sample_data(physeq)[,'author'],
                         variable_region = sample_data(physeq)[,'variable_region'],
                         age = sample_data(physeq)[,'host_age'],
                         bmi = sample_data(physeq)[,'host_bmi'])
dims.umap <- merge(as.data.frame(dims.umap), covariates, by="row.names") # umap+covariates
# saveRDS(object=dims.umap, file=filepath)




########
# PLOT #
########




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
  # ggsave(paste0("~/IBS/UMAP/data_logratio/umap", paste0(tax_rank, "_disease.jpg", sep=""), sep=""), width=6, height=4, type="cairo")
  
  # PLOT host_subtype
  ggplot(dims.umap, aes(x = UMAP_1, y = UMAP_2, color = host_subtype))+
    geom_point(size = 1, alpha = 0.7)+
    scale_color_manual(values = c('#99CCFF', '#FF3300', '#990000', '#FF66CC','#FFFF66', '#CCCCCC'))+
    labs(title = paste0('Agglomeration at taxonomic level:', tax_rank))+
    theme_bw()
  # ggsave(paste0("~/IBS/UMAP/data_logratio/umap", paste0(tax_rank, "_subtype.jpg", sep=""), sep=""), width=6, height=4, type="cairo")
  
  # PLOT author
  ggplot(dims.umap, aes(x = UMAP_1, y = UMAP_2, color = author))+
    geom_point(size = 2, alpha = 0.7)+
    labs(title = paste0('Agglomeration at taxonomic level:', tax_rank))+
    theme_bw()
  # ggsave(paste0("~/IBS/UMAP/data_logratio/umap", paste0(tax_rank, "_author.jpg", sep=""), sep=""), width=6, height=4, type="cairo")
  
  # PLOT seqtech
  ggplot(dims.umap, aes(x = UMAP_1, y = UMAP_2, color = sequencing_tech))+
    geom_point(size = 2, alpha = 0.7)+
    scale_color_manual(values = c('#6600FF', '#33CC33', '#006600', '#FF6633'))+
    labs(title = paste0('Agglomeration at taxonomic level:', tax_rank))+
    theme_bw()
  # ggsave(paste0("~/IBS/UMAP/data_logratio/umap", paste0(tax_rank, "_seqtech.jpg", sep=""), sep=""), width=6, height=4, type="cairo")
  # plotly::plot_ly(dims.umapFam_fecal_clr, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, color=~sequencing_tech, type="scatter3d", mode="markers")
  
  # PLOT variable region
  ggplot(dims.umap, aes(x = UMAP_1, y = UMAP_2, color = variable_region))+
    geom_point(size = 2, alpha = 0.7)+
    labs(title = paste0('Agglomeration at taxonomic level:', tax_rank))+
    theme_bw()
  # ggsave(paste0("~/IBS/UMAP/data_logratio/umap", paste0(tax_rank, "_vregion.jpg", sep=""), sep=""), width=6, height=4, type="cairo")
}








