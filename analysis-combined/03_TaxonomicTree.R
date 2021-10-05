##########################
# Purpose: Taxonomic tree plotting
# Date: August 2021
# Author: Salomé Carcy
##########################




##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(trac)
library(treeio)
library(ggtree)
library(RColorBrewer)

# Data
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
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
physeq.fecal <- subset_samples(physeq.all, sample_type == 'stool') # 2,220 samples
physeq.sigmoid <- subset_samples(physeq.all, sample_type == 'sigmoid') # 431 samples
cat("Nb of fecal samples:", nsamples(physeq.fecal))
cat("\nNb of sigmoid samples:", nsamples(physeq.sigmoid))




####################
# GET PHYLO OBJECT #
####################

physeq.glom <- tax_glom(physeq.all, "Genus")
# saveRDS(physeq.glom, "~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/physeq_all_glomGenus.rds")
# physeq.glom <- readRDS("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/physeq_all_glomGenus.rds")

## ***********************
## 1 - GET TAXONOMIC TABLE
tax <- physeq.glom@tax_table@.Data[,1:6]
# Add a "Life" column (rank0)
tax <- cbind("Life", tax)
colnames(tax)[1] <- "Life"
# add an ASV column
# tax <- cbind(tax, rownames(tax))
# colnames(tax)[ncol(tax)] <- "ASV"


#_____________________________________________________
# FOR UNAGGLOMERATED DATA (if we have some taxa with unassigned taxonomy)
# replace NA taxa by a character
# tax[is.na(tax[,"Class"]) == TRUE, "Class"] <- "c__"
# tax[is.na(tax[,"Order"]) == TRUE, "Order"] <- "o__"
# tax[is.na(tax[,"Family"]) == TRUE, "Family"] <- "f__"
# tax[is.na(tax[,"Genus"]) == TRUE, "Genus"] <- "g__"
# # tax[is.na(tax[,"Species"]) == TRUE, "Species"] <- "s__"
# 
# 
# # Make it so labels are unique, from Class (4th column) to Genus (7th column)
# for (i in 4:7) {
#   # add a number when the type is unknown... e.g. "g__"
#   NAtaxa <- grep("__", tax[,i], value = TRUE)
#   # 
#   if(length(NAtaxa) > 0){
#     tax[names(NAtaxa), i] <- paste0(tax[names(NAtaxa), i], 1:length(NAtaxa))
#   }
# }
#_____________________________________________________

# See how taxonomic table looks like (matrix)
tax[1:5,1:7]


## ***********************
## 2 - GET TAXONOMIC TABLE WITH CUMULATIVE LABELS
full_tax <- tax
for (i in 2:7) {
  full_tax[, i] <- paste(full_tax[, i-1], full_tax[, i], sep = "::")
}
full_tax <- as.data.frame(full_tax, stringsAsFactors = TRUE)
full_tax[1:5,1:7]


## ***********************
## 3 - GET PHYLO OBJECT
tree.phylo <- tax_table_to_phylo(~Life/Kingdom/Phylum/Class/Order/Family/Genus,
                                 data = full_tax, collapse = TRUE)




#######################
# PLOT TAXONOMIC TREE #
#######################

## ***********************
## 1 - GET TREEDATA OBJECT

# We will define new labels for the tips
d <- data.frame(label=full_tax$Genus,
                new_label=tax[,"Genus"],
                phylum=tax[,"Phylum"])
d[!d$phylum %in% c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota"), "phylum"] <- "Other"

# Create a treedata object, that will contain the new_label
tree.td <- full_join(tree.phylo, d, by="label")

# Get the list of tip labels belonging to each phylum (phyla are the names of the elements in the list)
l <- list()
for (phylum in unique(d$phylum)){
  print(phylum)
  l[[phylum]] <- d[d$phylum==phylum,"label"]
}

# Group the nodes and edges by phylum (to color the branches by phylum)
tree.td <- groupOTU(tree.td, l)


## ***********************
## 2 - PLOT <3
# Find the last 2 characters for alpha
rgb(0,255,0, max=255, alpha=0.6*255)

colors <- paste0(brewer.pal(6, "Set1"), "B2", sep="")
names(colors) <- c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota", "Other")

jpeg("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/taxTree_all.jpg", width=2000, height=2000, res=400)
ggtree(tree.td, layout="circular", aes(color=group))+
  # geom_tiplab(aes(label=new_label), size=0.5, color="black")+
  geom_nodepoint(size=0.1, color="black")+
  scale_color_manual(values=colors)+
  labs(color="Phylum")
dev.off()






