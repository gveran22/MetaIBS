##########################
# Purpose: Dimension reduction of all samples
# Date: January 2022
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(ggplot2)
#library(umap)
library(tidyverse)
#library(reshape2)
#library(gtools)




###############
# IMPORT DATA #
###############

# table(sample_sums(physeq.all)<500) # sanity check
path <- "~/IBS/dimension-reduction"
physeq.NZcomp <- readRDS(file.path(path, "physeq_all_NZcomp.rds"))
physeq.CSN <- readRDS(file.path(path, "physeq_all_CSN.rds"))
physeq.clr <- readRDS(file.path(path, "physeq_all_clr.rds"))


# Sanity checks
print(table(rowSums(otu_table(physeq.NZcomp)))) # check if there is any row not summing to 1
print(table(rowSums(otu_table(physeq.CSN)))) # check that all rows are summing to the same total
print(otu_table(physeq.clr)[1:5, 1:5]) # should all be relative




#######################
# DIMENSION REDUCTION #
#######################

# Compute distances
getDistances <- function(){
  print("Aitchison")
  glom.ait <- phyloseq::distance(physeq.clr, method = 'euclidean') # aitchison
  print("Bray")
  glom.bray <- phyloseq::distance(physeq.CSN, method = "bray") # bray-curtis
  print("Canberra")
  glom.can <- phyloseq::distance(physeq.NZcomp, method = "canberra") # canberra
  dist.list <- list("Ait" = glom.ait, "Bray" = glom.bray, "Canb" = glom.can)
  return(dist.list)
}

# Function to plot distances
plotDistances2D <- function(dlist, ordination="MDS", coloring="host_disease"){
  plist <- NULL
  plist <- vector("list", 3)
  names(plist) <- c("Aitchison", "Bray-Curtis", "Canberra")
  
  print("Aitchison")
  # Aitchison
  set.seed(123)
  iMDS.Ait <- ordinate(physeq=physeq.clr, method=ordination, distance=dlist$Ait)
  plist[[1]] <- plot_ordination(physeq.clr, iMDS.Ait, color=coloring)
  
  print("Bray")
  # Bray-Curtis
  set.seed(123)
  iMDS.Bray <- ordinate(physeq=physeq.CSN, method=ordination, distance=dlist$Bray)
  plist[[2]] <- plot_ordination(physeq.CSN, iMDS.Bray, color=coloring)
  
  print("Canberra")
  # Canberra
  set.seed(123)
  iMDS.Can <- ordinate(physeq=physeq.NZcomp, method=ordination, distance=dlist$Can)
  plist[[3]] <- plot_ordination(physeq.NZcomp, iMDS.Can, color=coloring)
  
  # Creating a dataframe to plot everything
  plot.df = plyr::ldply(plist, function(x) x$data)
  names(plot.df)[1] <- "distance"
  
  return(plot.df)
}




########
# PLOT #
########

dist.all <- getDistances()
saveRDS(dist.all, "~/IBS/dimension-reduction/output_ait-bray-can-distances.rds")
plot.df <- plotDistances2D(dlist=dist.all, ordination="MDS", coloring="host_disease")
saveRDS(plot.df, "~/IBS/dimension-reduction/output_plot-df-distances.rds")

# Plot
ggplot(plot.df, aes(Axis.1, Axis.2, color=host_disease))+
  geom_point(size=6, alpha=0.5)  + scale_color_manual(values = c('blue', 'red'))+
  facet_wrap(distance~., scales='free', nrow=1)+
  theme_bw()+
  theme(strip.text.x = element_text(size=20))+
  labs(color="Disease")
ggsave("~/IBS/dimension-reduction/output_plot.jpg", width=10, height=4, type="cairo")
