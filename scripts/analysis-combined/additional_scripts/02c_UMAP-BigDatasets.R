#
# Purpose: UMAP & BrayCurtis plotting on big datasets (Family agglom)
# Date: July 2022
# Author: Salomé Carcy
#




# 1. IMPORT ----------------------------------------------------------------------------------------

# Libraries
library(phyloseq)
library(ggplot2)
library(cowplot)
library(umap)
library(tidyverse)
library(reshape2)
library(gtools)


# Phyloseq objects
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
physeq.agp      <- readRDS(file.path(path.phy, "physeq_agp.rds"))
physeq.pozuelo  <- readRDS(file.path(path.phy, "physeq_pozuelo.rds"))
physeq.hugerth  <- readRDS(file.path(path.phy, "physeq_hugerth.rds"))




# 2. UMAP ----------------------------------------------------------------------------------------

# 2.1. Functions ===================================================

#______________________________________________________
# Function to agglomerate to given taxonomic level
aggTable <- function(physeq, tax_rank){
  
  # sanity check
  cat("\nSanity check ASV table:\n")
  test <- otu_table(physeq)[1:5,1:5]
  colnames(test) <- NULL
  print(test)
  
  cat("\nSanity check 2 ASV table: is there any ASV with 0 total count?\n")
  print(table(colSums(otu_table(physeq)) == 0))
  
  # Agglomerate
  cat("\n++ AGGLOMERATE TO", tax_rank, "LEVEL++\n")
  long.agg <- physeq %>%
    tax_glom(taxrank = tax_rank) %>%
    psmelt()
  
  # Get a matrix TaxRank (rows) x Samples (columns)
  if(tax_rank=="Phylum")     {agglomeratedTable <- acast(long.agg, Phylum ~ Sample, fun.aggregate=sum, value.var = 'Abundance')}
  else if(tax_rank=="Class") {agglomeratedTable <- acast(long.agg, Class ~ Sample,  fun.aggregate=sum, value.var = 'Abundance')}
  else if(tax_rank=="Order") {agglomeratedTable <- acast(long.agg, Order ~ Sample,  fun.aggregate=sum, value.var = 'Abundance')}
  else if(tax_rank=="Family"){agglomeratedTable <- acast(long.agg, Family ~ Sample, fun.aggregate=sum, value.var = 'Abundance')}
  else if(tax_rank=="Genus") {agglomeratedTable <- acast(long.agg, Genus ~ Sample,  fun.aggregate=sum, value.var = 'Abundance')}
  cat("-> Dimensions matrix:", dim(agglomeratedTable), "\n")
  cat("-> Any 0 counts?", table(agglomeratedTable == 0), "\n")
  
  # Sanity check reshape df
  randomTaxa   <- rownames(agglomeratedTable)[7]
  randomSample <- colnames(agglomeratedTable)[7]
  valueAgg <- sum(long.agg[long.agg[,tax_rank] == randomTaxa & long.agg$Sample == randomSample, "Abundance"])
  wideAgg  <- agglomeratedTable[randomTaxa, randomSample]
  cat("-> Sanity check: value in long shaped df is", valueAgg, "while in wide shaped df is", wideAgg, "\n")
  
  cat("-> Head of agglomerated matrix:\n")
  print(agglomeratedTable[1:10,1:5])
  
  cat("\n-> Add pseudocounts after agglomerating to", tax_rank, "level\n")
  agglomeratedTable[agglomeratedTable == 0] <- 0.5
  
  cat("\n->Head of agglomerated matrix:\n")
  print(agglomeratedTable[1:10,1:5])
  
  return(agglomeratedTable)
}


#______________________________________________________
# Function to obtain table with log-ratios between taxa
LogRatios <- function(abundanceTable, tax_rank){
  
  cat("\n++ GET LOG-RATIOS BETWEEN", tax_rank, "++\n")
  
  # Get combinations between all taxa, 1st column numerator, 2nd column denominator, 3rd column taxa1/taxa2
  comb <- as.data.frame(combinations(nrow(abundanceTable), 2, rownames(abundanceTable), repeats.allowed = FALSE))
  comb[,3] <- paste0(comb[,1], "/", comb[,2])
  
  cat("Number of combinations:", nrow(comb), "\n")
  
  # Compute the ratios
  ratios <- as.data.frame(abundanceTable) # taxa as rows, samples as columns
  test <- ratios[1:3,1:3] # for later sanity check
  
  ratios[ comb[,3] ,] <- mapply(function(x, y) log2(x/y), # compute the log ratio
                                abundanceTable[ as.character(comb[,1]) ,], # between first column
                                abundanceTable[ as.character(comb[,2]) ,]) # and second column of all combinations
  
  ratios <- ratios[-c(1:nrow(abundanceTable)),] # remove first x rows that are not ratios (it is simply the taxa)
  
  # sanity check
  cat("Sanity check: log-ratio is", log2(test[1,1]/test[2,1]),
      "and the corresponding value in the log-ratio table is", ratios[1,1], "\n")
  cat("Sanity check2: log-ratio is", log2(test[1,2]/test[3,2]),
      "and the corresponding value in the log-ratio table is", ratios[2,2], "\n")
  
  cat("We have", dim(ratios)[2], "samples and", dim(ratios)[1], "predictors (log-ratios) \n")
  
  # Return ratios df
  return(t(ratios))
  
}


#______________________________________________________
# Function to get the table to feed into the UMAP
getTable <- function(physeq, tax_rank){
  
  # Get the agglomerated taxa
  agglomeratedTable <- aggTable(physeq, tax_rank)
  
  # Get log-ratios between taxa
  ratiosTable <- LogRatios(agglomeratedTable, tax_rank)
  
  cat("\n=> Finished computing", tax_rank, "log-ratios\n")
  
  return(ratiosTable)
}


#______________________________________________________
# Function to compute the UMAP
getUMAP <- function(physeq, tax_rank){
  
  # Get the table with log-ratios between taxa
  cat("** 1 - GET LOG-RATIOS BETWEEN", tax_rank, "** \n")
  ratiosTable <- getTable(physeq, tax_rank)
  
  # Compute UMAP
  cat("\n\n** 2 - COMPUTE UMAP ** \n")
  cat("Nb of samples", nrow(ratiosTable), "and nb of predictors", ncol(ratiosTable), "\n")
  set.seed(123)
  umap <- uwot::umap(ratiosTable, # umap on samples (rows) and taxa ratios (columns)
                     n_neighbors=20, n_components=3, scale=F, n_threads=8)
  cat("\n=> Finished computing the UMAP \n")
  
  # Get the (x,y) coordinates from the UMAP
  dims.umap <- umap %>% as.data.frame()
  colnames(dims.umap) <- c("UMAP_1", "UMAP_2", "UMAP_3")
  cat("Verify that the rows in the UMAP coordinates df correspond to the rows in the logratios table (same samples): \n")
  print(table(rownames(dims.umap) == rownames(ratiosTable))) # sanity check
  
  # Add covariates
  covariates <- data.frame(disease         = sample_data(physeq)[,'host_disease'],
                           subtype         = sample_data(physeq)[,'host_subtype'],
                           sample_type     = sample_data(physeq)[,'sample_type'],
                           collection      = sample_data(physeq)[,'Collection'],
                           seq_tech        = sample_data(physeq)[,'sequencing_tech'],
                           author          = sample_data(physeq)[,'author'])
  cat("Verify that the rows in the UMAP coordinates df correspond to the rows in the covariates table: \n")
  print(table(rownames(dims.umap) == rownames(covariates))) # sanity check
  dims.umap <- merge(as.data.frame(dims.umap), covariates, by="row.names") # umap+covariates
  colnames(dims.umap)[1] <- "Sample"
  cat("\nNb of samples", nrow(dims.umap), "\n\n")
  
  return(dims.umap)
}



# 2.2. Compute UMAP ===============================================
umap.agp <- getUMAP(physeq=physeq.agp, tax_rank="Family")
umap.pozuelo <- getUMAP(physeq=physeq.pozuelo, tax_rank="Family")
# umap.pozuelo <- getUMAP(physeq=subset_samples(physeq.pozuelo, Collection=="1st"), tax_rank="Family")
umap.hugerth <- getUMAP(physeq=physeq.hugerth, tax_rank="Family")
# umap.hugerth <- getUMAP(physeq=subset_samples(physeq.hugerth, sample_type=="stool"), tax_rank="Family")



# 2.3. Plot UMAP ==================================================

# 2.3.1. AGP ####
a <- ggplot(umap.agp, aes(x = UMAP_1, y = UMAP_2, color = host_disease))+
  geom_point(size = 2)+
  scale_color_manual(values = c('blue', 'red', 'black'), name="")+ # disease
  labs(x="UMAP 1", y="UMAP 2", title = "AGP")+
  theme_cowplot()+
  theme(line = element_blank(),
        plot.title = element_text(hjust = 0.5))


# 2.3.2. Pozuelo ####
b <- ggplot(umap.pozuelo, aes(x = UMAP_1, y = UMAP_2, color = host_disease))+
  geom_point(size = 2)+
  scale_color_manual(values = c('blue', 'red', 'black'), name="")+ # disease
  labs(x="UMAP 1", y="UMAP 2", title = "Pozuelo")+
  theme_cowplot()+
  theme(line = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggplot(umap.pozuelo, aes(x = UMAP_1, y = UMAP_2, color = Collection))+
  geom_point(size = 2)+
  scale_color_manual(values = c('#f1a340', '#998ec3'), name="Time point")+
  labs(x="UMAP 1", y="UMAP 2", title = "Pozuelo")+
  theme_cowplot()+
  theme(line = element_blank(),
        plot.title = element_text(hjust = 0.5))


# 2.3.3. Hugerth ####
c <- ggplot(umap.hugerth, aes(x = UMAP_1, y = UMAP_2, color = host_disease))+
  geom_point(size = 2)+
  scale_color_manual(values = c('blue', 'red', 'black'), name="")+ # disease
  labs(x="UMAP 1", y="UMAP 2", title = "Hugerth")+
  theme_cowplot()+
  theme(line = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggplot(umap.hugerth, aes(x = UMAP_1, y = UMAP_2, color = sample_type))+
  geom_point(size = 2)+
  scale_color_manual(values = c('#7fbf7b', '#af8dc3'))+
  labs(x="UMAP 1", y="UMAP 2", title = "Hugerth")+
  theme_cowplot()+
  theme(line = element_blank(),
        plot.title = element_text(hjust = 0.5))


# 2.3.4. Combine plots ####
ggdraw() +
  draw_plot(a, x = 0,     y = 0, width = .32, height = 1) +
  draw_plot(b, x = .33,   y = 0, width = .32, height = 1) +
  draw_plot(c, x = .66,   y = 0, width = .32, height = 1)
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/02_DimensionReduction/plots_dimension-reduc/AGP-Poz-Hug_fam-agglom_UMAP.jpeg", width=15, height=5)




# 3. Bray-Curtis on PCoA -------------------------------------------------------------------------------

# 3.1. Functions ===================================================

getBray <- function(physeq, tax_rank){
  
  # If want to agglomerate to a certain taxonomic level, do it before normalization
  if(tax_rank != "none"){
    cat("\n++ Agglomerate to", tax_rank, "level ++\n")
    physeq <- tax_glom(physeq, taxrank=tax_rank)
  }
  
  # Normalize
  cat("++ Normalize ++\n")
  physeq.CSN <- transform_sample_counts(physeq, function(x) (x*min(sample_sums(physeq))) / sum(x) )
  print(table(rowSums(otu_table(physeq.CSN)))) # check that all rows are summing to the same total
  
  # Compute BrayCurtis
  cat("\n++ Compute Bray-Curtis distance ++\n")
  bray.dist <- phyloseq::distance(physeq.CSN, method = "bray")
  
  # Perform PCoA
  cat("++ Perform PCoA ++\n")
  set.seed(123)
  iMDS.Bray <- ordinate(physeq=physeq.CSN, method="PCoA", distance=bray.dist)
  plt.data <- plot_ordination(physeq.CSN, iMDS.Bray)$data
  
  return(plt.data)
}



# 3.2. Compute Bray-Curtis ==========================================

bray.agp <- getBray(physeq=physeq.agp, tax_rank="Family")
bray.pozuelo <- getBray(physeq=physeq.pozuelo, tax_rank="Family")
# bray.pozuelo <- getBray(physeq=subset_samples(physeq.pozuelo, Collection=="1st"), tax_rank="Family")
bray.hugerth <- getBray(physeq=physeq.hugerth, tax_rank="Family")
# bray.hugerth <- getBray(physeq=subset_samples(physeq.hugerth, sample_type=="stool"), tax_rank="Family")



# 3.3. Plot PCoA ==================================================

# 3.3.1. AGP ####
d <- ggplot(bray.agp, aes(x = Axis.1, y = Axis.2, color = host_disease))+
  geom_point(size = 1.5)+
  scale_color_manual(values = c('blue', 'red', 'black'), name="")+ # disease
  labs(title = "AGP")+
  theme_cowplot()+
  theme(line = element_blank(),
        plot.title = element_text(hjust = 0.5))


# 3.3.2. Pozuelo ####
e1 <- ggplot(bray.pozuelo, aes(x = Axis.1, y = Axis.2, color = host_disease))+
  geom_line(aes(group=host_ID), color="black", lwd=0.1) +
  geom_point(size = 1)+
  scale_color_manual(values = c('blue', 'red', 'black'), name="")+ # disease
  labs(title = "Pozuelo")+
  theme_cowplot()+
  theme(line = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.1,0.95),
        plot.title = element_text(hjust = 0.5))

e2 <- ggplot(bray.pozuelo, aes(x = Axis.1, y = Axis.2, color = Collection))+
  geom_line(aes(group=host_ID), color="black", lwd=0.1) +
  geom_point(size = 1)+
  scale_color_manual(values = c('#f1a340', '#998ec3'), name="Time point")+
  theme_cowplot()+
  theme(line = element_blank(),
        legend.position = c(0.1,0.95),
        plot.title = element_text(hjust = 0.5))

# 3.3.3. Hugerth ####
f1 <- ggplot(bray.hugerth, aes(x = Axis.1, y = Axis.2, color = host_disease))+
  geom_line(aes(group=host_ID), color="black", lwd=0.1) +
  geom_point(size = 1)+
  scale_color_manual(values = c('blue', 'red', 'black'), name="")+ # disease
  labs(title = "Hugerth")+
  theme_cowplot()+
  theme(line = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.1,0.95),
        plot.title = element_text(hjust = 0.5))

f2 <- ggplot(bray.hugerth, aes(x = Axis.1, y = Axis.2, color = sample_type))+
  geom_line(aes(group=host_ID), color="black", lwd=0.1) +
  geom_point(size = 1)+
  scale_color_manual(values = c('#7fbf7b', '#af8dc3'), name="Sample type")+
  theme_cowplot()+
  theme(line = element_blank(),
        axis.text = element_text(size=8),
        legend.position = c(0.1,0.95),
        plot.title = element_text(hjust = 0.5))


# 3.3.4. Combine plots ####
ggdraw() +
  draw_plot(d,  x = 0,    y = .1, width = .39, height = .8) +
  draw_plot(e1, x = .4,   y = .5, width = .29, height = .5) +
  draw_plot(e2, x = .4,   y = 0,  width = .29, height = .5) +
  draw_plot(f1, x = .7,   y = .5, width = .29, height = .5) +
  draw_plot(f2, x = .7,   y = 0,  width = .29, height = .5)
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/02_DimensionReduction/plots_dimension-reduc/AGP-Poz-Hug_fam-agglom_bray-PCoA.jpeg", width=15, height=5)




