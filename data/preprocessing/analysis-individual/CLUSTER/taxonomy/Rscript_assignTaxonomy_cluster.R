##########################
# Purpose: Assign taxonomy to ASVs in ASV table (on computer cluster)
# Date: January 2023
# Author: Salomé Carcy
##########################


# FILE PATHS (to change depending on computer cluster)
lib.loc <- "~/R/x86_64-pc-linux-gnu-library/3.6/" # R libraries location
cluster.wd <- "~/MetaIBS/Dada2/taxonomy"
input.data <- file.path(cluster.wd, "input")
output.data <- file.path(cluster.wd, "output")

input.data <- "~/Projects/MetaIBS/data_local/analysis-individual/CLUSTER/Taxonomy/input"
list.files(input.data)


#____________________________________________________________________
# IMPORT LIBRARIES
library(phyloseq, lib.loc=lib.loc)
library(dada2, lib.loc=lib.loc)
library(stats)
library(Biostrings)
library(base)


#____________________________________________________________________
# IMPORT DATA

# Import ASV tables into a list (name of each element is the file name)
print("... Importing ASV tables into list...")
datasets  <- list.files(input.data)
seqtables <- sapply(datasets, function(x) readRDS(file.path(input.data, x)), USE.NAMES=T, simplify=F)

# seqtable_ringel <- readRDS("~/IBS/Dada2/seqtable_ringel.rds")
# seqtable_labus <- readRDS("~/IBS/Dada2/seqtable_labus.rds")
# seqtable_lopresti <- readRDS("~/IBS/Dada2/seqtable_lopresti.rds")
# seqtable_zhuang <- readRDS("~/IBS/Dada2/seqtable_zhuang.rds")
# seqtable_zhu <- readRDS("~/IBS/Dada2/seqtable_zhu.rds")
# seqtable_nagel <- readRDS("~/IBS/Dada2/seqtable_nagel.rds")
# 
# # Create a vector
# print("... Putting ASV tables into vector ...")
# datasets <- c("ringel", "labus", "lopresti", "zhuang", "zhu", "nagel")
# 
# seqtables <- vector("list", length(datasets))
# names(seqtables) <- datasets
# 
# seqtables[[1]] <- seqtable_ringel
# seqtables[[2]] <- seqtable_labus
# seqtables[[3]] <- seqtable_lopresti
# seqtables[[4]] <- seqtable_zhuang
# seqtables[[5]] <- seqtable_zhu
# seqtables[[6]] <- seqtable_nagel


#____________________________________________________________________
# ASSIGN TAXONOMY WITH SILVA v138

# Save path to save taxonomic table

for(i in 1:length(seqtables)){
  # Print current dataset
  print(paste("++ Taxa assignment", datasets[i], "++\n"))
  
  # Initialize
  taxa <- NULL
  taxa.print <- NULL
  file_name <- NULL
  
  # Sanity check
  print(dim(seqtables[[i]]))
  
  # Assign taxonomy
  taxa <- assignTaxonomy(seqtables[[i]], file.path(cluster.wd, "silva-taxonomic-ref/silva_nr_v138_train_set.fa.gz"),
                         tryRC = TRUE, # try reverse complement of the sequences
                         multithread=TRUE, verbose = TRUE)
  taxa <- addSpecies(taxa, file.path(cluster.wd, "silva-taxonomic-ref/silva_species_assignment_v138.fa.gz"))
  
  # Sanity check
  taxa.print <- taxa
  rownames(taxa.print) <- NULL # removing sequence rownames for display only
  head(taxa.print)
  table(taxa.print[,1]) # Show the different kingdoms (should be only bacteria)
  table(taxa.print[,2]) # Show the different phyla
  
  # Save taxa table as [path]/taxonomy/output/taxa_NameDataset.rds
  file_name  <- paste0(paste0(paste0(output.data, "/taxa_"), # get [path]/taxonomy/output/taxa_
                              str_match(datasets[i], "seqtable_(.*?).rds")[,2]), # get NameDataset (in "seqtable_NameDataset.rds")
                       ".rds")
  saveRDS(taxa, file_name)
}