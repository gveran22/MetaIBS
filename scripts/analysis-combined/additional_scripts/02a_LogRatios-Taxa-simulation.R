#
# Purpose: Verify that log ratios are correctly computed
# Date: July 2022
# Author: Salomé Carcy
#




# 1. IMPORT ----------------------------------------------------------------------------------------

# Libraries
library(phyloseq)
library(tidyverse)
library(reshape2)
library(gtools)


# Create mock phyloseq object
otumat <- data.frame(cbind(c(1,2,4), c(2,4,6), c(10,2,5)))
colnames(otumat) <- c("spl1", "spl2", "spl3")
rownames(otumat) <- c("OTU1", "OTU2", "OTU3")
otumat

taxmat = matrix(sample(letters, 9, replace = F), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat

phytest <- phyloseq(otu_table(as.matrix(otumat), taxa_are_rows = T),
                    tax_table(as.matrix(taxmat)))




# 2. FUNCTIONS ----------------------------------------------------------------------------------------
# Note: functions were copied-pasted from the script 02a_LogRatios-Taxa.R
# The functions were slightly adjusted for the print(), as the mock phyloseq is only (3,3) dimensions compared to real data phyloseq objects


#_________________________________________________
# Function to agglomerate to given taxonomic level
aggTable <- function(physeq, tax_rank){
  
  # sanity check
  cat("\nSanity check ASV table:\n")
  # test <- otu_table(physeq)[1:5,1:5]
  # colnames(test) <- NULL
  # print(test)
  print(otu_table(physeq))
  
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
  randomTaxa   <- rownames(agglomeratedTable)[2]
  randomSample <- colnames(agglomeratedTable)[2]
  valueAgg <- sum(long.agg[long.agg[,tax_rank] == randomTaxa & long.agg$Sample == randomSample, "Abundance"])
  wideAgg  <- agglomeratedTable[randomTaxa, randomSample]
  cat("-> Sanity check: value in long shaped df is", valueAgg, "while in wide shaped df is", wideAgg, "\n")
  
  cat("-> Head of agglomerated matrix:\n")
  # print(agglomeratedTable[1:10,1:5])
  print(agglomeratedTable)
  
  cat("\n-> Add pseudocounts after agglomerating to", tax_rank, "level\n")
  agglomeratedTable[agglomeratedTable == 0] <- 0.5
  
  print("\n->Head of agglomerated matrix:\n")
  # print(agglomeratedTable[1:10,1:5])
  print(agglomeratedTable)
  
  return(agglomeratedTable)
}


#______________________________________________________
# Function to obtain table with log-ratios between taxa
LogRatios <- function(abundanceTable, tax_rank){
  
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
  
  # Return ratios df
  return(t(ratios))
}


#________________________________________________
# Function to get the table to feed into the UMAP
getTable <- function(physeq, tax_rank){
  
  # Get the agglomerated taxa
  agglomeratedTable <- aggTable(physeq, tax_rank)
  
  # Get log-ratios between taxa
  ratiosTable <- LogRatios(agglomeratedTable, tax_rank)
  
  cat("\n=> Finished computing", tax_rank, "log-ratios\n")
  
  return(ratiosTable)
}



# 3. TEST THE LOG-RATIOS COMPUTATION ----------------------------------------------------------------------------------------

# On phyloseq object with no 0s and no family agglomeration (all families are distinct)
table(as.data.frame(otu_table(phytest)) == 0)
table(as.data.frame(tax_table(phytest))$Family)
getTable(physeq = phytest, tax_rank = 'Family') # run


# On phyloseq object with some 0s and no family agglomeration (all families are distinct)
# Create phyloseq object with zeroes
phytest_with0 <- phytest
otu_table(phytest_with0)[1,1] <- 0
otu_table(phytest_with0)[3,2] <- 0

table(as.data.frame(otu_table(phytest_with0)) == 0)
table(as.data.frame(tax_table(phytest_with0))$Family)
getTable(physeq = phytest_with0, tax_rank = 'Family') # run


# On phyloseq object with some 0s and family agglomeration (some OTUs belong to the same family)
phytest_withagglom <- phytest_with0
tax_table(phytest_withagglom)[1,"Family"] <- tax_table(phytest_withagglom)[2,"Family"]

table(as.data.frame(otu_table(phytest_withagglom)) == 0)
table(as.data.frame(tax_table(phytest_withagglom))$Family)
getTable(physeq = phytest_withagglom, tax_rank = 'Family') # run
