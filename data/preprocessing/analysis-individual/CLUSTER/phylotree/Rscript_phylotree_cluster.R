##########################
# Purpose: Add phylogenetic trees to phyloseq objects
# Date: January 2023
# Author: Salomé Carcy
##########################


# FILE PATHS (to change depending on computer cluster)
lib.loc <- "~/R/x86_64-pc-linux-gnu-library/3.6/" # R libraries location on cluster
cluster.wd <- "~/MetaIBS/Dada2/phylotree"
input.data <- file.path(cluster.wd, "input")
output.data <- file.path(cluster.wd, "output")


#____________________________________________________________________
# IMPORT LIBRARIES
library(phyloseq)
library(dada2)
library(DECIPHER)
library(phangorn)
library(stats)
library(Biostrings)
library(base)


#____________________________________________________________________
# IMPORT DATA

# Import phyloseq objects into a list (name of each element is the file name)
print("... Importing phyloseq objects into list...")
datasets  <- list.files(input.data)
phyloseqobjects <- sapply(datasets, function(x) readRDS(file.path(input.data, x)), USE.NAMES=T, simplify=F)


#____________________________________________________________________
# FUNCTION TO COMPUTE PHYLOGENETIC TREE IN EACH DATASET/PHYLOSEQ OBJECT

phylotree <- function(physeq, namedataset){
  # Print current dataset/phyloseq object
  print("______________________________________________")
  print(paste("++ Phyloseq object", namedataset, "++\n"))
  
  # Initialize
  seqtable.nochim <- NULL
  seqs            <- NULL
  alignment       <- NULL
  phang.align     <- NULL
  dm              <- NULL
  treeNJ          <- NULL
  fit             <- NULL
  fitGTR          <- NULL
  
  
  # Get ASV table from phyloseq object
  seqtable.nochim <- as.matrix(as.data.frame(otu_table(physeq)))
  print(dim(seqtable.nochim)) # sanity check
  
  
  # Multiple-sequence alignment (know which regions are conserved/different to be able to do the phylogeny)
  print("-- Multiple-sequence alignment --")
  seqs <- getSequences(seqtable.nochim) # get the ASVs
  names(seqs) <- seqs
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA) # by default 70% of sequences must share a common region to anchor the alignment space
  
  # Construct a neighbor-joining tree
  print("-- Neighbor-joining tree --")
  phang.align <- phyDat(as(alignment, "matrix"), type="DNA") # transform the aligned DNA sequences into a phyDat (phangorn) object
  dm <- dist.ml(phang.align) # compute pairwise distance between DNA sequences
  treeNJ <- NJ(dm) # create a neighbor-joining tree estimation based on the distance matrix
  fit <- pml(treeNJ, data=phang.align) # get the likelihood of the phylogenetic tree given the sequence alignment (then we'll optimize it)
  
  ## negative edges length changed to 0!
  
  # Fit a generalized time-reversible with gamma rate variation (GTR+G+I)
  print("-- Optimize --")
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, # gamma rate and proportion of variable size get optimized
                      rearrangement = "stochastic", control = pml.control(trace = 0)) # (trace = 0) don't show output during optimization
  # stochastic tree rearrangement
  
  # Add tree to physeq
  print("-- Add to phyloseq object --")
  physeq <- merge_phyloseq(physeq, phy_tree(fitGTR$tree))

  
  # Give surnames to sequence variants & store the sequence variants in "refseq" in the phyloseq object
  dna <- DNAStringSet(taxa_names(physeq)) # get the sequence variants (ASVs)
  names(dna) <- taxa_names(physeq)
  physeq <- merge_phyloseq(physeq, dna) # store the dna sequences in the refseq of the phyloseq object
  taxa_names(physeq) <- paste0("ASV", seq(ntaxa(physeq))) # replace the whole dna sequences in the taxa_names by a surname ASV1, ASV2, etc.
  
  
  # Save phyloseq object as [path]/phylotree/output/physeq_NameDataset.rds
  file_name  <- paste0(paste0(paste0(output.data, "/physeq_"), # get [path]/phylotree/output/physeq_
                              namedataset), # NameDataset (in "physeq_NameDataset.rds")
                       ".rds")
  saveRDS(physeq, file_name)
}


#____________________________________________________________________
# RUN PHYLOGENETIC TREE ON SELECTED DATASETS

# Call function for whichever dataset you want (by default, all present in input/ subdirectory)
for(i in 1:length(phyloseqobjects)){
  
  # Initialize
  physeq <- NULL
  dtset <- NULL
  
  # get phyloseq object and name of the dataset
  physeq <- phyloseqobjects[[i]]
  dtset <- str_match(datasets[i], "physeq_(.*?).rds")[,2] # get NameDataset (in "physeq_NameDataset.rds")
  
  # compute phylogenetic tree
  phylotree(physeq=physeq, namedataset=dtset)
}

