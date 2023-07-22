# ************************************
# Purpose: Testing annotated heatmaps
# Date: July 2023
# Author: VT
# ************************************
# The script is mostly based on Salomes script 03_Heatmaps.R. The script aims
# to do a compositional mean test


# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(pals)


## 1.2. Data ####
path.root <- "~/Projects/MetaIBS" # CHANGE THIS ROOT DIRECTORY ON YOUR COMPUTER
path.data <- file.path(path.root, "data/analysis-combined/12_compositional_mean_test/")

path.phylobj    <- file.path(path.root, "data/phyloseq-objects/phyloseq-without-phylotree")
datasets        <- list.files(path.phylobj, pattern=".rds")
phyloseqobjects <- sapply(datasets, function(x) readRDS(file.path(path.phylobj, x)), USE.NAMES=T, simplify=F)
# names(phyloseqobjects) # sanity check



# ***********************
# 2. PREPROCESS DATA ####
# ***********************

## 2.1. Merge phyloseq objects ####
physeq.all <- merge_phyloseq(phyloseqobjects[[1]], phyloseqobjects[[2]]) # Merge first two phyloseq objects in the list
# if there are more than 2 phyloseq objects, merge the rest of them
if(length(phyloseqobjects)>2){
  for (i in 3:length(phyloseqobjects)){
    print(paste0("merging with phyloseq object #", i))
    physeq.all <- merge_phyloseq(physeq.all, phyloseqobjects[[i]])
  }
}


## 2.2. Separate fecal & sigmoid samples ####
physeq.fecal <- subset_samples(physeq.all, sample_type == 'stool') # 2,220 samples
physeq.fecal <- prune_taxa(taxa_sums(physeq.fecal)>0, physeq.fecal) # remove ASVs that are not present anymore
cat("Nb of fecal samples:", nsamples(physeq.fecal))

physeq.sigmoid <- subset_samples(physeq.all, sample_type == 'sigmoid') # 431 samples
physeq.sigmoid <- prune_taxa(taxa_sums(physeq.sigmoid)>0, physeq.sigmoid) # remove ASVs that are not present anymore
cat("\nNb of sigmoid samples:", nsamples(physeq.sigmoid))



# ********************************
# 4. HEATMAP FAMILIES AS ROWS ####
# ********************************

##___________________________
## 4.1. Agglomerate to Family level ####
# Download speedyseq here: https://github.com/mikemc/speedyseq
family.agg_phylo <- physeq.fecal %>%
  speedyseq::tax_glom(taxrank = "Family")

family.agg <- family.agg_phylo %>%
  psmelt()

# Identify families present in at least 3 datasets
list.family <- family.agg %>%
  filter(Abundance > 0) %>%
  group_by(Family) %>%
  summarise(nb_datasets = n_distinct(author)) %>%
  filter(nb_datasets>2) %>%
  ungroup()
list.family <- list.family$Family # 119 Families


# Agglomerate again at family level, but keeping only families present in at least 3 datasets
# family.agg <- subset_taxa(physeq.fecal, Family %in% list.family) %>%
#   tax_glom(taxrank = "Family") %>%
#   transform_sample_counts(function(x) {x/sum(x)} ) %>%
#   psmelt()

# Get dataframe family x samples
familyTable <- acast(family.agg %>% filter(Family %in% list.family),
                     Family ~ Sample, fun.aggregate=sum, value.var = 'Abundance')
family_phylo <- subset_taxa(family.agg_phylo, Family %in% list.family)


# Add 1 on zeros
otu_table(family_phylo)[otu_table(family_phylo)==0] <- 1

## Sanity checks

# test_taxa <- as.data.frame(family_phylo@tax_table@.Data) %>%
#   rownames_to_column("value") %>%
#   select(value, Family)
# test <- as.data.frame(as.matrix(family_phylo@otu_table@.Data))
# test_names <- enframe(colnames(test)) %>%
#   left_join(test_taxa)
# colnames(test) <- test_names$Family


# ******************************************************
# 5. Compositional mean test per author and disease ####
# ******************************************************
unique_authors <- unique(sample_data(family_phylo)$author)
unique_authors <- unique_authors[unique_authors != "Ringel"]

df_results <- data.frame(author = c(), p_value = c(), test_stat = c())
# iterate over all datasets and do a compositional mean test for each dataset individual
for (i in seq_len(length(unique_authors))) {
  # Filter the data
  
  authors <- unique_authors[i]
  family_phylo_subset <- subset_samples(
    family_phylo, sample_data(family_phylo)$author == authors)
  
# Code from https://github.com/bio-datascience/Causal_Microbiome_Tutorial/blob/main/3_mean_diff_test_AG/3.1_high_dim_test_AG.R
  

  
  
  # Separate the sample into two groups
  ps_control = subset_samples(family_phylo_subset, 
                              sample_data(family_phylo_subset)$host_disease == "Healthy")
  ps_disease = subset_samples(family_phylo_subset, 
                              sample_data(family_phylo_subset)$host_disease != "Healthy")
  
  n <- dim(otu_table(family_phylo_subset))[1]
  p <- dim(otu_table(family_phylo_subset))[2]
  
  X_c = otu_table(ps_control); dim(X_c)
  X_c = as(X_c, "matrix")
  
  X_t = otu_table(ps_disease); dim(X_t)
  X_t = as(X_t, "matrix")
  
  nx <- dim(X_c)[1]
  ny <- dim(X_t)[1]
  
  x <- X_c/(rowSums(X_c)%*%matrix(1,1,p))
  y <- X_t/(rowSums(X_t)%*%matrix(1,1,p))
  
  log_x<-log(x)
  log_y<-log(y)
  
  clog_TX <- log_x-1/p*rowSums(log_x)%*%matrix(1,1,p)
  clog_TY <- log_y-1/p*rowSums(log_y)%*%matrix(1,1,p)
  
  clr_x_mean <- colSums(clog_TX)/nx
  clr_y_mean <- colSums(clog_TY)/ny
  
  
  clr_x_var <- diag(var(clog_TX))*(nx-1)/nx
  clr_y_var <- diag(var(clog_TY))*(ny-1)/ny
  clr_x_stat_var <- (clr_x_var*nx + clr_y_var*ny)/(nx*ny)
  
  # p-values
  clr_x_stat <- max((clr_x_mean - clr_y_mean)^2/clr_x_stat_var)
  p_clrx <- 1-exp(-1/sqrt(pi)*exp(-(clr_x_stat-(2*log(p)-log(log(p))))/2))
  df_results <- bind_rows(df_results, c(author = authors, 
                          p_value = p_clrx, 
                          test_stat = clr_x_stat)) 
}

df_results$p_value_adj <- p.adjust(df_results$p_value, method = "BH")

write.csv(df_results, file = paste0(path.data, "comp_mean_test_IBS_Healthy.csv"))


# *****************************************
# 6. Compositional mean across authors ####
# *****************************************

unique_authors <- unique(sample_data(family.agg_phylo)$author)
# unique_authors <- unique_authors[unique_authors != "Ringel"]

author_comparision <- t(combn(unique_authors, 2))
otu_table(family.agg_phylo)[otu_table(family.agg_phylo)==0] <- 1

df_results <- data.frame(author = c(), p_value = c(), test_stat = c())
# iterate over all datasets and do a compositional mean test for each dataset individual
for (i in seq_len(nrow(author_comparision))) {
  # Filter the data
  
  author_1 <- author_comparision[i, 1]
  author_2 <- author_comparision[i, 2]
  family_phylo_subset <- subset_samples(
    family.agg_phylo, sample_data(family.agg_phylo)$author == author_1 |
      sample_data(family.agg_phylo)$author == author_2)
  
  # Code from https://github.com/bio-datascience/Causal_Microbiome_Tutorial/blob/main/3_mean_diff_test_AG/3.1_high_dim_test_AG.R
  
  
  
  
  # Separate the sample into two groups
  ps_control = subset_samples(family_phylo_subset, 
                              sample_data(family_phylo_subset)$author == author_1)
  ps_disease = subset_samples(family_phylo_subset, 
                              sample_data(family_phylo_subset)$author == author_2)
  
  n <- dim(otu_table(family_phylo_subset))[1]
  p <- dim(otu_table(family_phylo_subset))[2]
  
  X_c = otu_table(ps_control); dim(X_c)
  X_c = as(X_c, "matrix")
  
  X_t = otu_table(ps_disease); dim(X_t)
  X_t = as(X_t, "matrix")
  
  nx <- dim(X_c)[1]
  ny <- dim(X_t)[1]
  
  x <- X_c/(rowSums(X_c)%*%matrix(1,1,p))
  y <- X_t/(rowSums(X_t)%*%matrix(1,1,p))
  
  log_x<-log(x)
  log_y<-log(y)
  
  clog_TX <- log_x-1/p*rowSums(log_x)%*%matrix(1,1,p)
  clog_TY <- log_y-1/p*rowSums(log_y)%*%matrix(1,1,p)
  
  clr_x_mean <- colSums(clog_TX)/nx
  clr_y_mean <- colSums(clog_TY)/ny
  
  
  clr_x_var <- diag(var(clog_TX))*(nx-1)/nx
  clr_y_var <- diag(var(clog_TY))*(ny-1)/ny
  clr_x_stat_var <- (clr_x_var*nx + clr_y_var*ny)/(nx*ny)
  
  # p-values
  clr_x_stat <- max((clr_x_mean - clr_y_mean)^2/clr_x_stat_var)
  p_clrx <- 1-exp(-1/sqrt(pi)*exp(-(clr_x_stat-(2*log(p)-log(log(p))))/2))
  df_results <- bind_rows(df_results, c(author = paste0(author_1, "_", author_2), 
                                        p_value = p_clrx, 
                                        test_stat = clr_x_stat)) 
}

write.csv(df_results, file = paste0(path.data, "comp_mean_test_across_datasets.csv"))


