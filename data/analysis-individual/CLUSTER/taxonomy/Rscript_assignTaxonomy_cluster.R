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


#____________________________________________________________________
# IMPORT LIBRARIES
library(phyloseq)
library(dada2)
library(stats)
library(Biostrings)
library(base)


#____________________________________________________________________
# IMPORT DATA

# Import ASV tables into a list (name of each element is the file name)
print("... Importing ASV tables into list...")
datasets  <- list.files(input.data, pattern=".rds")
seqtables <- sapply(datasets, function(x) readRDS(file.path(input.data, x)), USE.NAMES=T, simplify=F)


#____________________________________________________________________
# ASSIGN TAXONOMY WITH SILVA v138

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
  set.seed(123)
  taxa <- assignTaxonomy(seqtables[[i]], file.path(cluster.wd, "silva-taxonomic-ref/silva_nr99_v138.1_train_set.fa.gz"),
                         tryRC = TRUE, # try reverse complement of the sequences
                         multithread=TRUE, verbose = TRUE)
  set.seed(123)
  taxa <- addSpecies(taxa, file.path(cluster.wd, "silva-taxonomic-ref/silva_species_assignment_v138.1.fa.gz"))
  
  # Sanity check
  taxa.print <- taxa
  rownames(taxa.print) <- NULL # removing sequence rownames for display only
  print(head(taxa.print))
  print(table(taxa.print[,1], useNA="ifany")) # Show the different kingdoms (should be only bacteria)
  print(table(taxa.print[,2], useNA="ifany")) # Show the different phyla
  
  # Save taxa table as [path]/taxonomy/output/taxa_NameDataset.rds
  file_name  <- paste0(paste0(paste0(output.data, "/taxa_"), # get [path]/taxonomy/output/taxa_
                              stringr::str_match(datasets[i], "seqtablenochim_(.*?).rds")[,2]), # get NameDataset (in "seqtable_NameDataset.rds")
                       ".rds")
  saveRDS(taxa, file_name)
}