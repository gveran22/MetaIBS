##########################
# Purpose: Obtain log-ratios between taxa
# Date: March 2022
# Author: Salomé Carcy
##########################




##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(tidyverse)
library(reshape2)
library(gtools)

# Data
path.phy <- "~/IBS/PhyloTree/input"
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
cat("\n++ MERGE PHYLOSEQ OBJECTS ++\n")
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
cat("Nb of taxa:", ntaxa(physeq.all), "\n\n\n") # should be 79,943




#############
# FUNCTIONS #
#############

#________________________________________________________
# Function to agglomerate to given taxonomic level
aggTable <- function(physeq, tax_rank, add_pseudocounts="after"){
  
  if(add_pseudocounts=="before"){
    cat("\n-> Add pseudocounts", add_pseudocounts, "agglomerating to", tax_rank, "level\n")
    otu_table(physeq)[otu_table(physeq) == 0] <- 0.5
  }
  
  # sanity check
  cat("\nSanity check ASV table:\n")
  test <- otu_table(physeq)[1:5,1:5]
  colnames(test) <- NULL
  print(test)
  
  # Agglomerate
  cat("\n++ AGGLOMERATE TO", tax_rank, "LEVEL++\n")
  long.agg <- physeq %>%
    tax_glom(taxrank = tax_rank) %>%
    psmelt()
  
  # Get a matrix TaxRank (rows) x Samples (columns)
  if(tax_rank=="Phylum")     {agglomeratedTable <- acast(long.agg, Phylum ~ Sample, fun.aggregate=sum, value.var = 'Abundance')}
  else if(tax_rank=="Class") {agglomeratedTable <- acast(long.agg, Class ~ Sample,  fun.aggregate=sum, value.var = 'Abundance')}
  else if(tax_rank=="Order") {agglomeratedTable <- acast(long.agg, Order ~ Sample,  fun.aggregate=sum, value.var = 'Abundance')}
  else if(tax_rank=="Family"){agglomeratedTable <- acast(long.agg, Family ~ Sample, fun.aggregate=sum, value.var = 'Abundance')}
  else if(tax_rank=="Genus") {agglomeratedTable <- acast(long.agg, Genus ~ Sample,  fun.aggregate=sum, value.var = 'Abundance')}
  cat("-> Dimensions matrix:", dim(agglomeratedTable), "\n")
  cat("-> Any 0 counts?", table(agglomeratedTable == 0), "\n")
  
  # Sanity check reshape df
  randomTaxa   <- rownames(agglomeratedTable)[7]
  randomSample <- colnames(agglomeratedTable)[7]
  valueAgg <- sum(long.agg[long.agg[,tax_rank] == randomTaxa & long.agg$Sample == randomSample, "Abundance"])
  wideAgg  <- agglomeratedTable[randomTaxa, randomSample]
  cat("-> Sanity check: value in long shaped df is", valueAgg, "while in wide shaped df is", wideAgg, "\n")
  
  cat("-> Head of agglomerated matrix:\n")
  print(agglomeratedTable[1:10,1:5])
  
  if(add_pseudocounts=="after"){
    cat("\n-> Add pseudocounts", add_pseudocounts, "agglomerating to", tax_rank, "level\n")
    agglomeratedTable[agglomeratedTable == 0] <- 0.5
  }
  
  cat("\n->Head of agglomerated matrix:\n")
  print(agglomeratedTable[1:10,1:5])
  
  return(agglomeratedTable)
}



#________________________________________________________
# Function to obtain table with log-ratios between taxa
LogRatios <- function(abundanceTable, tax_rank, add_pseudocounts="after"){
  
  cat("\n++ GET LOG-RATIOS BETWEEN", tax_rank, "++\n")
  
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
  
  # Save
  if(add_pseudocounts=="before"){filepath <- paste0("~/IBS/LogRatios-taxa/pseudocounts_bef-agg/ratios", paste0(tax_rank, ".rds", sep=""), sep="")}
  else if(add_pseudocounts=="after"){filepath <- paste0("~/IBS/LogRatios-taxa/pseudocounts_aft-agg/ratios", paste0(tax_rank, ".rds", sep=""), sep="")}
  
  # Save table with samples as rows and logratios as columns
  saveRDS(object=t(ratios), file=filepath)
}



#________________________________________________________
# Function to get the table to feed into the UMAP
getTable <- function(physeq, tax_rank, add_pseudocounts){
  
  # Get the agglomerated taxa
  agglomeratedTable <- aggTable(physeq, tax_rank, add_pseudocounts)
  
  # Get log-ratios between taxa
  ratiosTable <- LogRatios(agglomeratedTable, tax_rank, add_pseudocounts)
  
  cat("\n=> Finished computing", tax_rank, "log-ratios\n")
}






######################
# COMPUTE LOG RATIOS #
######################

taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")

# Get log-ratios while adding pseudocounts before agglomeration
cat("*****************************************\n")
cat("*** PSEUDOCOUNTS BEFORE AGGLOMERATION *** \n")
cat("*****************************************\n")
for(taxa in taxranks){
  cat("\n >>>>", taxa, "level <<<<\n")
  getTable(physeq=physeq.all, tax_rank=taxa, add_pseudocounts="before")
}

# Get log-ratios while adding pseudocounts after agglomeration
cat("\n\n\n*****************************************\n")
cat("*** PSEUDOCOUNTS AFTER AGGLOMERATION *** \n")
cat("*****************************************\n")
for(taxa in taxranks){
  cat("\n >>>>", taxa, "level <<<<\n")
  getTable(physeq=physeq.all, tax_rank=taxa, add_pseudocounts="after")
}



