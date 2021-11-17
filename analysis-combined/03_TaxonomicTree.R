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
library(ggnewscale)
library(RColorBrewer)
library(scales)

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
physeq.fecal <- prune_taxa(taxa_sums(physeq.fecal)>0, physeq.fecal) # remove ASVs that are not present anymore

physeq.sigmoid <- subset_samples(physeq.all, sample_type == 'sigmoid') # 431 samples
physeq.sigmoid <- prune_taxa(taxa_sums(physeq.sigmoid)>0, physeq.sigmoid) # remove ASVs that are not present anymore

cat("Nb of fecal samples:", nsamples(physeq.fecal))
cat("\nNb of sigmoid samples:", nsamples(physeq.sigmoid))




####################
# DEFINE FUNCTIONS #
####################

#______________________
# Get a treedata object
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
# Function to get dataframe with information on number of ASVs detected for each genus
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
count_per_genus <- function(physeq){
  # Get relative counts (common-scale normalization: sum per sample==100)
  cat("Common-scale normalization...\n")
  physeq.CSN <- transform_sample_counts(physeq, function(x) (x*100) / sum(x) )
  print(table(sample_sums(physeq.CSN))) # sanity check total of 100 counts per sample
  
  # Get the mean count per sample for each ASV
  cat("\nAverage count of each genus per sample....\n")
  meanCountASV <- t(data.frame(as.list(colMeans(otu_table(physeq.CSN)))))
  # sum(meanCountASV[,1]) # sanity check should sum to 100
  
  # Get a dataframe with the average count of each Genus per sample
  countGenus.df <- meanCountASV %>%
    as.data.frame()%>%
    rownames_to_column("ASV") %>%
    rename(Count=V1) %>%
    # Add taxonomy
    left_join(as.data.frame(tax_table(physeq.CSN)) %>%
                rownames_to_column("ASV"),
              by="ASV") %>%
    # get a column with Life::Kingdom::Phylum::Class::Order::Family::Genus
    mutate(leaf = paste("Life", Kingdom, Phylum, Class, Order, Family, Genus, sep = "::")) %>%
    # get average count per sample for each Genus
    group_by(leaf, Genus) %>%
    summarize(count=sum(Count)) %>% # if several ASVs belong to same Genus, need to sum their counts
    # about 30,851 ASVs belong to an unknown genus (which represents ~10% of the counts)
    filter(!is.na(Genus))
  
  countGenus.df <- data.frame(count=countGenus.df$count,
                              row.names = countGenus.df$leaf)
  
  # Return dataframe
  return(countGenus.df)
}


#______________________
# Function to get a dataframe with information on number of datasets where this genus is present
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




##############
# PLOT TREES #
##############

colors <- paste0(brewer.pal(6, "Dark2"), "99", sep="")
names(colors) <- c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota", "Other")

## ***********************
## 1 - ALL SAMPLES

# Agglomerate to Genus level
physeq.glom <- tax_glom(physeq.all, "Genus")
# saveRDS(physeq.glom, "~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/physeq_all_glomGenus.rds")
# physeq.glom <- readRDS("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/physeq_all_glomGenus.rds")

# Get the treedata object, and info on nb of ASV per genus, average count of each genus per sample.
tree.all <- phyloseq_to_treedata(physeq=physeq.glom)
asvDF.all <- nb_ASV_per_genus(physeq=physeq.all)
countDF.all <- count_per_genus(physeq=physeq.all)
datasetDF.all <- dataset_per_genus(physeq=physeq.glom)
sampleType.df <- physeq.glom %>%
  # melt to long format
  psmelt() %>%
  # keep only ASVs with a count != 0 in this specific sample
  filter(Abundance != 0) %>%
  # get a column with Life::Kingdom::Phylum::Class::Order::Family::Genus
  mutate(leaf = paste("Life", Kingdom, Phylum, Class, Order, Family, Genus, sep = "::")) %>%
  # get the number of datasets with non-0 count for each genus, and also the sample_type
  group_by(leaf, Genus) %>%
  summarize(spltype = paste(sort(unique(sample_type)), collapse=","))
sampleType.df <- data.frame(spltype=df$spltype,
                            row.names = df$leaf)


p1 <- ggtree(tree.all, layout="circular", aes(color=group), lwd=0.4)+
  # geom_tiplab(aes(label=new_label), size=2, color="black")+
  # geom_nodepoint(size=0.1, color="black")+
  scale_color_manual(values=colors)+
  labs(color="Phylum")

p2 <- gheatmap(p1, asvDF.all, offset=0.5, width=.08, colnames=FALSE)+
  scale_fill_gradient2(low = muted("red"),
                       mid = "white",
                       high = "#2b8cbe",
                       trans="log",
                       breaks=c(1,10,100,1000),
                       name="Number of ASVs")

p3 <- gheatmap(p2+new_scale_fill(),
               countDF.all, offset=1.1, width=.08, colnames = FALSE)+
  scale_fill_gradient(low = "white",
                      high = muted("red"),
                      trans="log",
                      breaks=c(1e-5,1e-3,1e-1,10),
                      labels = parse(text=c("10^-5","10^-3","10^-1", "10")),
                      name="Mean relative abundance (%)")+
  labs(title="All samples")

p4 <- gheatmap(p3+new_scale_fill(),
         datasetDF.all, offset=1.7, width=.08, colnames = FALSE)+
  scale_fill_gradient(low = "#ffffe5",
                      high = "#ec7014",
                      breaks = c(1,5,10),
                      name="Number of datasets")+
  labs(title="All samples")

jpeg("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/taxTree_all.jpg", width=6000, height=4000, res=400)
gheatmap(p4+new_scale_fill(),
         sampleType.df, offset=2.3, width=.08, colnames = FALSE)+
  scale_fill_manual(values=c("#f0f0f0","#636363","#bdbdbd"), name="Sample type")+
  labs(title="All samples")
dev.off()



## ***********************
## 2 - FECAL SAMPLES

# Agglomerate to Genus level
physeq.fecal.glom <- tax_glom(physeq.fecal, "Genus")
# saveRDS(physeq.fecal.glom, "~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/physeq_fecal_glomGenus.rds")
# physeq.fecal.glom <- readRDS("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/physeq_fecal_glomGenus.rds")

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
                       name="Number of ASVs")

p3_fecal <- gheatmap(p2_fecal+new_scale_fill(),
         countDF.fecal, offset=1.1, width=.08, colnames = FALSE)+
  scale_fill_gradient(low = "white",
                      high = muted("red"),
                      trans="log",
                      breaks=c(1e-5,1e-3,1e-1,10),
                      name="Average count per sample (%)")+
  labs(title="Fecal samples")

jpeg("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/taxTree_fecal.jpg", width=6000, height=3000, res=400)
gheatmap(p3_fecal+new_scale_fill(),
         datasetDF.fecal, offset=1.7, width=.08, colnames = FALSE)+
  scale_fill_gradient(low = "#ffffe5",
                      high = "#ec7014",
                      limits=c(1,13),
                      breaks = c(1,5,10),
                      name="Number of datasets")+
  labs(title="Fecal samples")
dev.off()



## ***********************
## 3 - SIGMOID SAMPLES

# Agglomerate to Genus level
physeq.sigmoid.glom <- tax_glom(physeq.sigmoid, "Genus")
# saveRDS(physeq.sigmoid.glom, "~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/physeq_sigmoid_glomGenus.rds")
# physeq.sigmoid.glom <- readRDS("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/physeq_sigmoid_glomGenus.rds")

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
                       name="Number of ASVs")

p3_sigmoid <- gheatmap(p2_sigmoid+new_scale_fill(),
         countDF.sigmoid, offset=1.1, width=.08, colnames = FALSE)+
  scale_fill_gradient(low = "white",
                      high = muted("red"),
                      trans="log",
                      breaks=c(1e-5,1e-3,1e-1,10),
                      name="Average count per sample (%)")+
  labs(title="Sigmoid samples")

jpeg("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/taxTree_sigmoid.jpg", width=6000, height=3000, res=400)
gheatmap(p3_sigmoid+new_scale_fill(),
         datasetDF.sigmoid, offset=1.7, width=.08, colnames = FALSE)+
  scale_fill_gradient(low = "#ffffe5",
                      high = "#ec7014",
                      limits=c(1,13),
                      breaks = c(1,5,10),
                      name="Number of datasets")+
  labs(title="Sigmoid samples")
dev.off()






















