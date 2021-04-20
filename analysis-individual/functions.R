##################################################################################
# This R script stores common functions used systematically for all datasets




#############
# DISTANCES #
#############


#____________________________________________________________________________________
# This function computes a number of distances provided by the phyloseq package
# and plots them in 2D by multi-dimensional scaling
#____________________________________________________________________________________

plotAllDistances <- function(physeq_obj, ait_dist, ordination="MDS"){
  
  library(phyloseq)
  
  # Get the distance methods
  dist_methods <- unlist(distanceMethodList)
  dist_methods <- dist_methods[-which(dist_methods=="ANY")] # Remove the user-defined distance
  dist_methods <- dist_methods[-c(1,15,18)] #remove unweighted unifrac, mountford, chao (not gonna work on pseudocounts)
  dist_methods <- dist_methods[1:16] #keep only the 16 first methods (remove the beta diversity distances)
  
  # Create a list where plots will be saved
  plist <- vector("list", length(dist_methods)+1)
  names(plist) <- dist_methods
  names(plist)[17] <- "aitchison" # add aitchison
  
  # Loop through all distance methods
  for(i in dist_methods){
    # Calculate distance matrix
    set.seed(123)
    iDist <- phyloseq::distance(physeq_obj, method=i)
    # Calculate ordination
    set.seed(123)
    iMDS  <- ordinate(physeq_obj, ordination, distance=iDist)
    ## Make plot
    p <- NULL
    p <- plot_ordination(physeq_obj, iMDS, color="host_disease") +
      ggtitle(paste("MDS using distance method ", i, sep=""))
    # Save the graphic to the plot list
    plist[[i]] = p
  }
  
  # Add Aitchison
  iMDS  <- ordinate(physeq_obj, ordination, distance=ait_dist)
  p <- NULL
  p <- plot_ordination(physeq_obj, iMDS, color="host_disease")
  p <- p + ggtitle("MDS using distance method Aitchison")
  plist[[17]] = p
  
  # Create a dataframe to plot everything
  plot.df = ldply(plist, function(x) x$data)
  names(plot.df)[1] <- "distance"
  
  # Plot
  p.alldist <-  ggplot(plot.df, aes(Axis.1, Axis.2, color=host_disease))+
    geom_point(size=5, alpha=0.5)  + scale_color_manual(values = c('blue', 'red'))+
    facet_wrap(~distance, scales='free')+
    theme(strip.text = element_text(size = 40),
          legend.text = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20))
  
  return(p.alldist)
}












