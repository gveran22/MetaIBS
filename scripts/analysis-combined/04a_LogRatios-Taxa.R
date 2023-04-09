# ***************************************
# Purpose: Obtain log-ratios between taxa
# Date: March 2022
# Author: Salomé Carcy
# Note: this script is highly recommended to run on a cluster!!!
# ***************************************




# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(phyloseq)
library(tidyverse)
library(reshape2)
library(gtools)

## 1.2. Data ####
path.root   <- "~/Projects/MetaIBS" # CHANGE THIS ROOT DIRECTORY ON YOUR CLUSTER
# path.output <- file.path(path.root, "data/analysis-combined/04a_LogRatios-Taxa/pseudocounts_aft-agg") # for authors of MetaIBS paper
path.output <- file.path(path.root, "data/analysis-combined/04a_LogRatios-Taxa/output") # CHANGE NAME OF THE OUTPUT FOLDER WHERE YOU WANT TO SAVE THE MATRICES

path.phylobj    <- file.path(path.root, "data/phyloseq-objects/phyloseq-without-phylotree")
datasets        <- list.files(path.phylobj)
phyloseqobjects <- sapply(datasets, function(x) readRDS(file.path(path.phylobj, x)), USE.NAMES=T, simplify=F)
# names(phyloseqobjects) # sanity check




# ***********************
# 2. PREPROCESS DATA ####
# ***********************

# Merge phyloseq objects
physeq.all <- merge_phyloseq(phyloseqobjects[[1]], phyloseqobjects[[2]]) # Merge first two phyloseq objects in the list
# if there are more than 2 phyloseq objects, merge the rest of them
if(length(phyloseqobjects)>2){
  for (i in 3:length(phyloseqobjects)){
    print(paste0("merging with phyloseq object #", i))
    physeq.all <- merge_phyloseq(physeq.all, phyloseqobjects[[i]])
  }
}

cat("Nb of samples:", nsamples(physeq.all), "\n") # should be 2,651
cat("Nb of taxa:", ntaxa(physeq.all), "\n\n\n") # should be 79,943





# *****************
# 3. FUNCTIONS ####
# *****************

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
  
  print("\n->Head of agglomerated matrix:\n")
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
  # if(add_pseudocounts=="before"){filepath <- paste0(path.output, "/pseudocounts_bef-agg/ratios", tax_rank, ".rds")}
  if(add_pseudocounts=="after"){filepath <- paste0(path.output, "/ratios", tax_rank, ".rds")}
  
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






# **************************
# 3. COMPUTE LOG RATIOS ####
# **************************

taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")

# # Get log-ratios while adding pseudocounts before agglomeration (bad idea because then the pseudocounts add up during agglomeration!!!)
# cat("*****************************************\n")
# cat("*** PSEUDOCOUNTS BEFORE AGGLOMERATION *** \n")
# cat("*****************************************\n")
# for(taxa in taxranks){
#   cat("\n >>>>", taxa, "level <<<<\n")
#   getTable(physeq=physeq.all, tax_rank=taxa, add_pseudocounts="before")
# }

# Get log-ratios while adding pseudocounts after agglomeration
cat("\n\n\n*****************************************\n")
cat("*** PSEUDOCOUNTS AFTER AGGLOMERATION *** \n")
cat("*****************************************\n")
for(taxa in taxranks){
  cat("\n >>>>", taxa, "level <<<<\n")
  getTable(physeq=physeq.all, tax_rank=taxa, add_pseudocounts="after")
}



