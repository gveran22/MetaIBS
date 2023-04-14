# *********************************
# Purpose: Exporting ASV sequences in fasta format for annotIEM
# Date: November 2021
# Author: Salomé Carcy
# *********************************




# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(Biostrings)


## 1.2. Data ####
path.root <- "~/Projects/MetaIBS" # CHANGE THIS ROOT DIRECTORY ON YOUR COMPUTER

path.asvtables <- file.path(path.root, "data/analysis-individual/CLUSTER/taxonomy/input") # these are the asv tables used for taxonomic alignment
datasets       <- list.files(path.asvtables, pattern=".rds", include.dirs=FALSE)
asvtables      <- sapply(datasets, function(x) readRDS(file.path(path.asvtables, x)), USE.NAMES=T, simplify=F)

# Change name of phyloseq objects to make it easier later on
names(asvtables)# sanity check
names(asvtables) <- gsub("seqtable_", "", names(asvtables))
names(asvtables) <- gsub(".rds", "", names(asvtables))
names(asvtables)# sanity check




# ***************************************
# 2. EXPORT TO FASTA FORMAT THE ASVs ####
# ***************************************

path.output <- file.path(path.root, "data/asv-sequences/")

for (author in names(asvtables)){
  print(author)
  print(dim(asvtables[[author]]))
  
  # Store ASVs sequences as DNAStringSet
  asv <- DNAStringSet(colnames(asvtables[[author]]))
  
  # Save ASVs sequences in FASTA file
  file.name <- paste0(path.output, "asv_", author, ".fasta")
  writeXStringSet(x=asv, filepath=file.name, compress=FALSE, format="fasta")
  
}

