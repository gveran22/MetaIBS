##########################
# Purpose: Add phylogenetic tree to Zhu phyloseq object
# Date: February 2023
# Author: Salomé Carcy
##########################


#____________________________________________________________________
# IMPORT LIBRARIES

library(phyloseq)
library(dada2)
library(DECIPHER)
library(phangorn)
library(stats)
library(Biostrings)
library(base)

# CHANGE THIS PATH (for your cluster or your local computer)
path.root <- "/grid/wsbs/home_norepl/scarcy/IBS/PhyloTree"


#____________________________________________________________________
# IMPORT DATA
physeq <- readRDS(file.path(path.root, "input/physeq_zhu.rds"))
seqtable.nochim <- as.matrix(as.data.frame(otu_table(physeq)))

# Sanity check
print(dim(seqtable.nochim))


#____________________________________________________________________
# PHYLOGENETIC TREE

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


#____________________________________________________________________
# GIVE SURNAMES TO OTUs

# Give surnames to sequence variants & store the sequence variants in "refseq" in the phyloseq object
dna <- DNAStringSet(taxa_names(physeq)) # get the sequence variants (ASVs)
names(dna) <- taxa_names(physeq)
physeq <- merge_phyloseq(physeq, dna) # store the dna sequences in the refseq of the phyloseq object
taxa_names(physeq) <- paste0("ASV_zhu_", seq(ntaxa(physeq))) # replace the whole dna sequences in the taxa_names by a surname ASV_zhu_1, ASV_zhu_2, etc.

# Save physeq object
saveRDS(physeq, file.path(path.root, "output/physeq_zhu.rds"))
