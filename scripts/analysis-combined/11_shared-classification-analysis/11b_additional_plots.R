# *****************************************************************************
# Purpose: Generate Figures for the paper
# Date: April 2023
# Author: (Minh) Viet Tran, based on Salomés Code
# *****************************************************************************


# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(phyloseq)
library(tidyverse)
library(trac)
library(treeio)
library(ggtree)
library(ggnewscale)
library(RColorBrewer)
library(reshape2)
library(ggupset)
library(cowplot)

## 1.2. Data ####
path_root <- "~/MetaIBS" # CHANGE THIS ROOT DIRECTORY ON YOUR COMPUTER
path_intermed <- paste0(path_root,
                        "/data/analysis-combined/11_shared_classification/")
path_phylobj <- 
  file.path(path_root,
            "/data/phyloseq-objects")
path_phylobj_nagel_pozuelo <- 
  file.path(path_root,
            "/data/analysis-combined/05_Common-ASVs/physeq_commonASV_Nagel-Pozuelo.rds")

physeq.ringel <- readRDS(file.path(path_phylobj, "physeq_ringel.rds"))
physeq.labus <- readRDS(file.path(path_phylobj, "physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path_phylobj, "physeq_lopresti.rds"))
physeq.pozuelo <- readRDS(file.path(path_phylobj, "physeq_pozuelo.rds"))
physeq.zhuang <- readRDS(file.path(path_phylobj, "physeq_zhuang.rds"))
physeq.zhu <- readRDS(file.path(path_phylobj, "physeq_zhu.rds"))
physeq.hugerth <- readRDS(file.path(path_phylobj, "physeq_hugerth.rds"))
physeq.fukui <- readRDS(file.path(path_phylobj, "physeq_fukui.rds"))
physeq.mars <- readRDS(file.path(path_phylobj, "physeq_mars.rds"))
physeq.liu <- readRDS(file.path(path_phylobj, "physeq_liu.rds"))
physeq.agp <- readRDS(file.path(path_phylobj, "physeq_agp.rds"))
physeq.nagel <- readRDS(file.path(path_phylobj, "physeq_nagel.rds"))
physeq.zeber <- readRDS(file.path(path_phylobj, "physeq_zeber.rds"))

## 1.4. Preprocessing ####
#Rename the ASV names to the actual string so its possible to to merge them
taxa_names(physeq.ringel) <- refseq(physeq.ringel)
taxa_names(physeq.labus) <- refseq(physeq.labus)
taxa_names(physeq.lopresti) <- refseq(physeq.lopresti)
taxa_names(physeq.pozuelo) <- refseq(physeq.pozuelo)
taxa_names(physeq.zhuang) <- refseq(physeq.zhuang)
taxa_names(physeq.zhu) <- refseq(physeq.zhu)
taxa_names(physeq.hugerth) <- refseq(physeq.hugerth)
taxa_names(physeq.fukui) <- refseq(physeq.fukui)
taxa_names(physeq.mars) <- refseq(physeq.mars)
#Build new phyloseq objects without the phylogenetic tree
physeq.ringel <- phyloseq(otu_table(physeq.ringel),
                          tax_table(physeq.ringel),
                          sample_data(physeq.ringel))
physeq.labus <- phyloseq(otu_table(physeq.labus),
                         tax_table(physeq.labus),
                         sample_data(physeq.labus))
physeq.lopresti <- phyloseq(otu_table(physeq.lopresti),
                            tax_table(physeq.lopresti),
                            sample_data(physeq.lopresti))
physeq.pozuelo <- phyloseq(otu_table(physeq.pozuelo),
                           tax_table(physeq.pozuelo),
                           sample_data(physeq.pozuelo))
physeq.zhuang <- phyloseq(otu_table(physeq.zhuang),
                          tax_table(physeq.zhuang),
                          sample_data(physeq.zhuang))
physeq.zhu <- phyloseq(otu_table(physeq.zhu),
                       tax_table(physeq.zhu),
                       sample_data(physeq.zhu))
physeq.hugerth <- phyloseq(otu_table(physeq.hugerth),
                           tax_table(physeq.hugerth),
                           sample_data(physeq.hugerth))
physeq.fukui <- phyloseq(otu_table(physeq.fukui),
                         tax_table(physeq.fukui),
                         sample_data(physeq.fukui))
physeq.mars <- phyloseq(otu_table(physeq.mars),
                        tax_table(physeq.mars),
                        sample_data(physeq.mars))
physeq.liu <- phyloseq(otu_table(physeq.liu),
                       tax_table(physeq.liu),
                       sample_data(physeq.liu))
physeq.agp <- phyloseq(otu_table(physeq.agp),
                       tax_table(physeq.agp),
                       sample_data(physeq.agp))
physeq.nagel <- phyloseq(otu_table(physeq.nagel),
                         tax_table(physeq.nagel),
                         sample_data(physeq.nagel))
physeq.zeber <- phyloseq(otu_table(physeq.zeber),
                         tax_table(physeq.zeber),
                         sample_data(physeq.zeber))
# merge all datasets
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

datasets <- list("Labus"    = physeq.labus,
                 "Lopresti" = physeq.lopresti,
                 "Ringel"   = physeq.ringel,
                 "AGP"      = physeq.agp,
                 "Liu"      = physeq.liu,
                 "Pozuelo"  = physeq.pozuelo,
                 "Fukui"    = physeq.fukui,
                 "Mars"     = physeq.mars,
                 "Hugerth"  = physeq.hugerth,
                 "Zhu"      = physeq.zhu,
                 "Zhuang"   = physeq.zhuang,
                 "Nagel"    = physeq.nagel,
                 "Zeber"    = physeq.zeber)

asv.df <- melt(lapply(datasets, function(x) taxa_names(x)))
colnames(asv.df) <- c("asv", "author")
length(unique(asv.df$asv)) # we do find 79,943 unique sequences

## 1.5. Plot number of shared ASVs ####


# Let's see which datasets share the exact same ASV sequence
common.asv <- asv.df %>%
  group_by(asv) %>%
  summarize(Datasets = list(unique(author))) %>%
  filter(lengths(Datasets)>1)


common_asvs_plot <- ggplot(common.asv, aes(x = Datasets))+
  geom_bar(fill = "#737373") +
  scale_x_upset() +
  theme_cowplot() +
  labs(y = "# common ASVs")

common_asvs_plot

ggsave(filename = paste0(path_intermed, "common_asvs.jpg"), 
       plot = common_asvs_plot, width = 6, 
       height = 5, device = 'jpeg', dpi = 700)
saveRDS(common_asvs_plot, file = paste0(path_intermed, "common_asvs.RData"))

## 1.6. Plot previously shown tree ####

data_nagel_pozuelo <- 
  readRDS(path_phylobj_nagel_pozuelo)
data_nagel_pozuelo <- subset_taxa(data_nagel_pozuelo , Kingdom == "Bacteria")
data_nagel_pozuelo <- subset_samples(data_nagel_pozuelo, Collection == "1st")
data_nagel_pozuelo <- subset_samples(data_nagel_pozuelo, 
                                     host_subtype %in% c("HC", "IBS-D"))

# ASV have appear at least in 10% of the observations
freq <- colSums(sign(data_nagel_pozuelo@otu_table@.Data))
data_nagel_pozuelo <- 
  phyloseq::prune_taxa(freq >= 0.1 * phyloseq::nsamples(data_nagel_pozuelo),
                       data_nagel_pozuelo)
names_nagel_pozuelo <- taxa_names(data_nagel_pozuelo)

data_nagel_pozuelo <- tax_glom(data_nagel_pozuelo, taxrank = "Genus")
names_nagel_pozuelo <- tax_table(data_nagel_pozuelo)@.Data

colors <- paste0(brewer.pal(6, "Dark2"), "99", sep = "")
names(colors) <- c("Actinobacteriota", "Bacteroidota", 
                   "Firmicutes", "Proteobacteria", "Verrucomicrobiota", "Other")

physeq.glom <- tax_glom(physeq.all, "Genus")

phyloseq_to_treedata_highlight <- function(physeq, comparison){
  ## ***********************
  ## 1 - GET TAXONOMIC TABLE
  tax <- physeq@tax_table@.Data[,1:6]
  # Add a "Life" column (rank0)
  tax <- cbind("Life", tax)
  colnames(tax)[1] <- "Life"
  # tax[1:5,1:7]
  
  
  ## ***********************
  ## 2 - GET TAXONOMIC TABLE WITH CUMULATIVE LABELS
  full_tax <- tax
  for (i in 2:7) {
    full_tax[, i] <- paste(full_tax[, i-1], full_tax[, i], sep = "::")
  }
  full_tax <- as.data.frame(full_tax, stringsAsFactors = TRUE)
  # full_tax[1:5,1:7]
  
  
  ## ***********************
  ## 3 - GET PHYLO OBJECT
  tree.phylo <- tax_table_to_phylo(~Life/Kingdom/Phylum/Class/Order/Family/Genus,
                                   data = full_tax, collapse = TRUE)
  
  
  ## ***********************
  ## 4 - GET TREEDATA OBJECT
  
  # We will define new labels for the tips
  d <- data.frame(label=full_tax$Genus,
                  new_label=tax[,"Genus"],
                  phylum=tax[,"Phylum"])
  d[!d$phylum %in% c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota"), "phylum"] <- "Other"
  d$comparision <- d$new_label %in% as.data.frame(comparison)$Genus
  # Create a treedata object, that will contain the new_label
  tree.td <- full_join(tree.phylo, d, by = "label")
  
  # Highlight selected nodes by
  # checking if the node is part of the selected tips for coloring the 
  # selected nodes
  # get the names of the nodes
  nodes_names <- str_replace(tree.td@phylo$node.label, "'", "")
  nodes_names <- str_replace(nodes_names, "'", "")
  # iterate over all nodes and check if they contain at least one selected tip
  comp_node <- c()
  for (i in seq_len(length(nodes_names))) {
    comp_node[i] <- 
      (sum(str_detect(as.character(d$label[d$comparision == TRUE]), 
                      nodes_names[i])) > 0)
  }
  tree.td@data$comparision[is.na(tree.td@data$comparision)] <- comp_node
  # Get the list of tip labels belonging to each phylum (phyla are the names of the elements in the list)
  l <- list()
  for (phylum in unique(d$phylum)){
    print(phylum)
    l[[phylum]] <- d[d$phylum == phylum, "label"]
  }
  present_yn <- list()
  for (compa in unique(d$comparision)){
    print(compa)
    present_yn[[paste0(compa)]] <- d[d$comparision==compa,"label"]
  }
  
  
  # Group the nodes and edges by phylum (to color the branches by phylum)
  tree.td <- groupOTU(tree.td, l)
  tree.td <- groupOTU(tree.td, present_yn, "selected")
  
  ## ***********************
  ## 5 - RETURN TREEDATA OBJECT
  return(tree.td)
}

tree.all <- phyloseq_to_treedata_highlight(physeq = physeq.glom, 
                                           comparison = names_nagel_pozuelo)

p1 <- ggtree(tree.all, layout = "circular", aes(color = group, 
                                                alpha = comparision), 
             lwd = 0.35) +
  scale_color_manual(values = colors) +
  labs(color = "Phylum",
       alpha = "Present in subset?")
p1 <- p1 %>%
  rotate_tree(-270)

ggsave(filename = paste0(path_intermed, "highlight_tree.jpg"), 
       plot = p1, width = 6, 
       height = 5, device = 'jpeg', dpi = 700)
saveRDS(p1, file = paste0(path_intermed, "highlight_tree.RData"))

        