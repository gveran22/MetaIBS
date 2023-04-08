# *********************************
# Purpose: Taxonomic tree plotting
# Date: August 2021
# Author: Salomé Carcy
# *********************************




# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(trac)
library(treeio)
library(ggtree)
library(ggnewscale)
library(RColorBrewer)
library(scales)


## 1.2. Data ####
path.root <- "~/Projects/MetaIBS" # CHANGE THIS ROOT DIRECTORY ON YOUR COMPUTER

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


# Separate fecal & sigmoid samples
# physeq.fecal <- subset_samples(physeq.all, sample_type == 'stool') # 2,220 samples
# physeq.fecal <- prune_taxa(taxa_sums(physeq.fecal)>0, physeq.fecal) # remove ASVs that are not present anymore
# 
# physeq.sigmoid <- subset_samples(physeq.all, sample_type == 'sigmoid') # 431 samples
# physeq.sigmoid <- prune_taxa(taxa_sums(physeq.sigmoid)>0, physeq.sigmoid) # remove ASVs that are not present anymore
# 
# cat("Nb of fecal samples:", nsamples(physeq.fecal))
# cat("\nNb of sigmoid samples:", nsamples(physeq.sigmoid))




# ************************
# 3. DEFINE FUNCTIONS ####
# ************************

#___________________________
# 3.1. Get a treedata object
phyloseq_to_treedata <- function(physeq){
  
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
  ## 5 - RETURN TREEDATA OBJECT
  return(tree.td)
}


#______________________
# 3.2. Function to get dataframe with information on number of ASVs detected for each genus
nb_ASV_per_genus <- function(physeq){
  
  asvGenus.df <- as.data.frame(tax_table(physeq)) %>%
    rownames_to_column("ASV") %>%
    # get a column with Life::Kingdom::Phylum::Class::Order::Family::Genus
    mutate(leaf = paste("Life", Kingdom, Phylum, Class, Order, Family, Genus, sep = "::")) %>%
    select(ASV, Genus, leaf) %>%
    # get nb of ASV per genus
    group_by(leaf, Genus) %>%
    summarize(n_ASV = n_distinct(ASV)) %>%
    # about 30,851 ASVs belong to an unknown genus
    filter(!is.na(Genus))
  
  asvGenus.df <- data.frame(n_ASV=asvGenus.df$n_ASV,
                            row.names = asvGenus.df$leaf)
  
  # Return df
  return(asvGenus.df)
}


#______________________
# Function to get a dataframe with information on average count per sample for each genus
# count_per_genus <- function(physeq){
#   # Get relative counts (common-scale normalization: sum per sample==100)
#   cat("Common-scale normalization...\n")
#   physeq.CSN <- transform_sample_counts(physeq, function(x) (x*100) / sum(x) )
#   print(table(sample_sums(physeq.CSN))) # sanity check total of 100 counts per sample
# 
#   # Get the mean count per sample for each ASV
#   cat("\nAverage count of each genus per sample....\n")
#   meanCountASV <- t(data.frame(as.list(colMeans(otu_table(physeq.CSN)))))
#   # sum(meanCountASV[,1]) # sanity check should sum to 100
# 
#   # Get a dataframe with the average count of each Genus per sample
#   countGenus.df <- meanCountASV %>%
#     as.data.frame()%>%
#     rownames_to_column("ASV") %>%
#     rename(Count=V1) %>%
#     # Add taxonomy
#     left_join(as.data.frame(tax_table(physeq.CSN)) %>%
#                 rownames_to_column("ASV"),
#               by="ASV") %>%
#     # get a column with Life::Kingdom::Phylum::Class::Order::Family::Genus
#     mutate(leaf = paste("Life", Kingdom, Phylum, Class, Order, Family, Genus, sep = "::")) %>%
#     # get average count per sample for each Genus
#     group_by(leaf, Genus) %>%
#     summarize(count=sum(Count)) %>% # if several ASVs belong to same Genus, need to sum their counts
#     # about 30,851 ASVs belong to an unknown genus (which represents ~10% of the counts)
#     filter(!is.na(Genus))
# 
#   countGenus.df <- data.frame(count=countGenus.df$count,
#                               row.names = countGenus.df$leaf)
# 
#   # Return dataframe
#   return(countGenus.df)
# }


#______________________
# 3.3. Function to get a dataframe with information on number of datasets where this genus is present
dataset_per_genus <- function(physeq){

  print("CAUTION: give the agglomerated phyloseq object for shorter running time")
  # Get a dataframe with the number of dataset(s) each Genus is found in
  datasetGenus.df <- physeq %>%
    # melt to long format
    psmelt() %>%
    # keep only ASVs with a count != 0 in this specific sample
    filter(Abundance != 0) %>%
    # get a column with Life::Kingdom::Phylum::Class::Order::Family::Genus
    mutate(leaf = paste("Life", Kingdom, Phylum, Class, Order, Family, Genus, sep = "::")) %>%
    # get the number of datasets with non-0 count for each genus
    group_by(leaf, Genus) %>%
    summarize(n_dataset = n_distinct(author))
  
  datasetGenus.df <- data.frame(n_dataset=datasetGenus.df$n_dataset,
                                row.names = datasetGenus.df$leaf)

  # Return dataframe
  return(datasetGenus.df)
}




# ******************
# 4. PLOT TREES ####
# ******************

path.data <- file.path(path.root, "data/analysis-combined/01_TaxonomicTree")
colors <- paste0(brewer.pal(6, "Dark2"), "99", sep="")
names(colors) <- c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota", "Other")

## *********************
## 4.1. All samples ####

# Agglomerate to Genus level
physeq.glom <- tax_glom(physeq.all, "Genus")
# saveRDS(physeq.glom, file.path(path.data, "physeq_all_glomGenus.rds")) # recommend to save it as takes a while to compute
# physeq.glom <- readRDS(file.path(path.data, "physeq_all_glomGenus.rds"))

# Get the treedata object, and info on nb of ASV per genus, average count of each genus per sample, sample types containing the genus.
tree.all <- phyloseq_to_treedata(physeq=physeq.glom)
asvDF.all <- nb_ASV_per_genus(physeq=physeq.all)
# countDF.all <- count_per_genus(physeq=physeq.all)
datasetDF.all <- dataset_per_genus(physeq=physeq.glom)
sampleType.df <- physeq.glom %>%
  # melt to long format
  psmelt() %>%
  # keep only ASVs with a count != 0
  filter(Abundance != 0) %>%
  # get a column with Life::Kingdom::Phylum::Class::Order::Family::Genus
  mutate(leaf = paste("Life", Kingdom, Phylum, Class, Order, Family, Genus, sep = "::")) %>%
  # get the sample_type(s) each genus can be found in
  group_by(leaf, Genus) %>%
  summarize(spltype = paste(sort(unique(sample_type)), collapse=","))
sampleType.df <- data.frame(spltype=sampleType.df$spltype,
                            row.names=sampleType.df$leaf)


# Plot tree (all genera)
p1 <- ggtree(tree.all, layout="circular", aes(color=group), lwd=0.4)+
  # geom_tiplab(aes(label=new_label), size=2, color="black")+
  scale_color_manual(values=colors)+
  labs(color="Phylum")


# Add outer layers
p2 <- gheatmap(p1, asvDF.all, offset=0.5, width=.08, colnames=FALSE)+
  scale_fill_gradient2(low = muted("red"),
                       mid = "white",
                       high = "#2b8cbe",
                       trans="log",
                       breaks=c(1,10,100,1000),
                       name="Number of ASVs")

# p3 <- gheatmap(p2+new_scale_fill(),
#                countDF.all, offset=1.1, width=.08, colnames = FALSE)+
#   scale_fill_gradient(low = "white",
#                       high = muted("red"),
#                       trans="log",
#                       breaks=c(1e-5,1e-3,1e-1,10),
#                       labels = parse(text=c("10^-5","10^-3","10^-1", "10")),
#                       name="Mean relative abundance (%)")+
#   labs(title="All samples")

p3 <- gheatmap(p2+new_scale_fill(),
         datasetDF.all, offset=1.1, width=.08, colnames = FALSE)+
  scale_fill_gradient(low = "#ffffe5",
                      high = "#ec7014",
                      breaks = c(1,5,10),
                      name="Number of datasets")

gheatmap(p3+new_scale_fill(),
         sampleType.df, offset=1.7, width=.08, colnames = FALSE)+
  scale_fill_manual(values=c("#f0f0f0","#636363","#bdbdbd"), name="Sample type")+
  labs(title="All samples")+
  scale_y_reverse() # flip tree up/down
ggsave(file.path(path.data, "taxTree_all.jpg"), width=12, height=8)



## ***********************
## 4.2. Fecal samples ####

# Agglomerate to Genus level
physeq.fecal.glom <- tax_glom(physeq.fecal, "Genus")
# saveRDS(physeq.fecal.glom, file.path(path.data, "physeq_fecal_glomGenus.rds")) # recommend to save it as takes a while to compute
# physeq.fecal.glom <- readRDS(file.path(path.data, "physeq_fecal_glomGenus.rds"))

# Get the treedata object, and info on nb of ASV per genus, average count of each genus per sample.
tree.fecal <- phyloseq_to_treedata(physeq=physeq.fecal.glom)
asvDF.fecal <- nb_ASV_per_genus(physeq=physeq.fecal)
countDF.fecal <- count_per_genus(physeq=physeq.fecal)
datasetDF.fecal <- dataset_per_genus(physeq=physeq.fecal.glom)

p1_fecal <- ggtree(tree.fecal, layout="circular", aes(color=group), lwd=0.4)+
  scale_color_manual(values=colors)+
  labs(color="Phylum")

p2_fecal <- gheatmap(p1_fecal, asvDF.fecal, offset=0.5, width=.08, colnames=FALSE)+
  scale_fill_gradient2(low = muted("red"),
                       mid = "white",
                       high = "#2b8cbe",
                       trans="log",
                       breaks=c(1,10,100,1000),
                       limits=c(1,4200),
                       name="Number of ASVs")

p3_fecal <- gheatmap(p2_fecal+new_scale_fill(), countDF.fecal, offset=1.1, width=.08, colnames = FALSE)+
  scale_fill_gradient(low = "white",
                      high = muted("red"),
                      trans="log",
                      breaks=c(1e-5,1e-3,1e-1,10),
                      limits=c(1e-6, 25),
                      name="Average count per sample (%)")+
  labs(title="Fecal samples")

gheatmap(p3_fecal+new_scale_fill(),
         datasetDF.fecal, offset=1.7, width=.08, colnames = FALSE)+
  scale_fill_gradient(low = "#ffffe5",
                      high = "#ec7014",
                      limits=c(1,13),
                      breaks = c(1,5,10),
                      name="Number of datasets")+
  scale_y_reverse() +# flip tree up/down
  labs(title="Fecal samples")
ggsave(file.path(path.data, "taxTree_fecal.jpg"), width=12, height=8)



## *************************
## 4.3. Sigmoid samples ####

# Agglomerate to Genus level
physeq.sigmoid.glom <- tax_glom(physeq.sigmoid, "Genus")
# saveRDS(physeq.sigmoid.glom, file.path(path.data, "physeq_sigmoid_glomGenus.rds")) # recommend to save it as takes a while to compute
# physeq.sigmoid.glom <- readRDS(file.path(path.data, "physeq_sigmoid_glomGenus.rds"))

# Get the treedata object, and info on nb of ASV per genus, average count of each genus per sample.
tree.sigmoid <- phyloseq_to_treedata(physeq=physeq.sigmoid.glom)
asvDF.sigmoid <- nb_ASV_per_genus(physeq=physeq.sigmoid)
countDF.sigmoid <- count_per_genus(physeq=physeq.sigmoid)
datasetDF.sigmoid <- dataset_per_genus(physeq=physeq.sigmoid.glom)

p1_sigmoid <- ggtree(tree.sigmoid, layout="circular", aes(color=group), lwd=0.4)+
  scale_color_manual(values=colors)+
  labs(color="Phylum")

p2_sigmoid <- gheatmap(p1_sigmoid, asvDF.sigmoid, offset=0.5, width=.08, colnames=FALSE)+
  scale_fill_gradient2(low = muted("red"),
                       mid = "white",
                       high = "#2b8cbe",
                       trans="log",
                       breaks=c(1,10,100,1000),
                       limits=c(1,4200),
                       name="Number of ASVs")

p3_sigmoid <- gheatmap(p2_sigmoid+new_scale_fill(), countDF.sigmoid, offset=1.1, width=.08, colnames = FALSE)+
  scale_fill_gradient(low = "white",
                      high = muted("red"),
                      trans="log",
                      breaks=c(1e-5,1e-3,1e-1,10),
                      limits=c(1e-6, 25),
                      name="Average count per sample (%)")+
  labs(title="Sigmoid samples")

gheatmap(p3_sigmoid+new_scale_fill(),
         datasetDF.sigmoid, offset=1.7, width=.08, colnames = FALSE)+
  scale_fill_gradient(low = "#ffffe5",
                      high = "#ec7014",
                      limits=c(1,13),
                      breaks = c(1,5,10),
                      name="Number of datasets")+
  scale_y_reverse() +
  labs(title="Sigmoid samples")
ggsave(file.path(path.data, "taxTree_sigmoid.jpg", width=12, height=8))

