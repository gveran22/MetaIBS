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

# Data
path.phy <- "~/Projects/IBS_Meta-analysis_16S/phyloseq-objects"
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



#############################
# TESTING STUFF WITH PHYLOG #
#############################

# With MicrobiotaProcess
library(MicrobiotaProcess)
library(ggtree)

# data("kostic2012crc")
# taxa <- phyloseq::tax_table(kostic2012crc)
# tree <- as.treedata(taxa)
# ggtree(tree, layout="circular", size=0.2) +
#   geom_tiplab(size=1)

taxa <- tax_table(physeq.labus)
tree <- as.treedata(taxa)

ggtree(tree, layout='circular') + 
  aes(color="Phylum")


# With taxo2phylog
library(ade4)
data(taxo.eg)
taxo.eg[[1]]
as.taxo(taxo.eg[[1]])


tax <- as.data.frame(physeq.labus@tax_table@.Data, stringsAsFactors = TRUE)
tax$Kingdom <- NULL
tax$Species <- NULL
as.taxo(tax[,1:2]) # problem

tax.phy <- taxo2phylog(tax)




# With as.phylo
library(treeio)
library(ggtree)
data("GNI2014", package="treemap")
n <- GNI2014[, c(3,1)]
n[,1] <- as.character(n[,1])
n[,1] <- gsub("\\s\\(.*\\)", "", n[,1]) # n is a dataframe with continent and country code

w <- cbind("World", as.character(unique(n[,1])))
colnames(w) <- colnames(n) # w is a dataframe with world & continent (but same column name?...)

edgelist <- rbind(n, w)

y <- as.phylo(edgelist)
ggtree(y, layout='circular') %<+% GNI2014 + 
  aes(color=continent) + geom_tippoint(aes(size=population), alpha=.6) + 
  geom_tiplab(aes(label=country), offset=.1) +
  theme(plot.margin=margin(60,60,60,60))

# try with labus
tax <- as.data.frame(physeq.labus@tax_table@.Data, stringsAsFactors = FALSE, row.names=NULL)

test1 <- unique( tax[,2:3] )
rownames(test1) <- NULL
test1[is.na(test1$Class),"Class"] <- paste0("c__", 1:2)
test2 <- cbind("Bacteria", as.character(unique(tax[,2])))
colnames(test2) <- colnames(test1)

test <- rbind(test1, test2)
testplot <- as.phylo(test)

ggtree(testplot, layout='circular') %<+% tax + 
  aes(color=Phylum) +
  geom_tippoint(size=1, alpha=.6) +
  geom_tiplab(aes(label="Class"), offset=.1)


###########################
# TESTING STUFF WITH TRAC #
###########################

library(trac)

tax <- physeq.labus@tax_table@.Data
# Add a "Life" column (rank0)
tax <- cbind("Life", tax); colnames(tax)[1] <- "Rank0"
# add an ASV column
tax <- cbind(tax, rownames(tax))
colnames(tax)[ncol(tax)] <- "ASV"

# replace NA taxa by a character
tax[is.na(tax[,"Class"]) == TRUE, "Class"] <- "c__"
tax[is.na(tax[,"Order"]) == TRUE, "Order"] <- "o__"
tax[is.na(tax[,"Family"]) == TRUE, "Family"] <- "f__"
tax[is.na(tax[,"Genus"]) == TRUE, "Genus"] <- "g__"
tax[is.na(tax[,"Species"]) == TRUE, "Species"] <- "s__"


# make it so labels are unique
for (i in seq(2, 8)) {
  # add a number when the type is unknown... e.g. "g__"
  # ii <- nchar(tax[, i]) == 3
  ii <- grep("__", tax[,i], value = TRUE)
  # if (sum(ii) > 0){
  if(length(ii) > 0){
    # tax[ii, i] <- paste0(tax[ii, i], 1:sum(ii))
    tax[names(ii), i] <- paste0(tax[names(ii), i], 1:length(ii))
  }
}


# cumulative labels are harder to read but easier to work with:
for (i in 2:9) {
  tax[, i] <- paste(tax[, i-1], tax[, i], sep = "::")
}
tax <- as.data.frame(tax, stringsAsFactors = TRUE)


###### DOESN'T WORK (R SESSION ABORTS) ##########
# form phylo object:
# tree1 <- tax_table_to_phylo(~Rank0/Rank1/Rank2/Rank3/Rank4/Rank5/Rank6/Rank7/ASV,
#                             data = tax, collapse = TRUE)
tree1 <- tax_table_to_phylo(~Rank0/Kingdom/Phylum/Class/Order/Family/Genus/Species/ASV,
                            data = tax, collapse = TRUE)


# convert this to an A matrix to be used for aggregation:
A <- phylo_to_A(tree1)

dat <- list(y = y[keep],
            x = t(agp@otu_table@.Data[, keep]),
            tree = tree1,
            tax = tax,
            A = A,
            sample_data = as_tibble(sample_data(agp))[keep, ])
# rows of A correspond to OTUs as do columns of x
# rearrange columns of x to be in the order of rows of A:
dat$x <- dat$x[, match(str_match(rownames(A), "::([^:]+)$")[, 2],
                       colnames(dat$x))]
identical(str_match(rownames(A), "::([^:]+)$")[,2],
          colnames(dat$x))
saveRDS(dat, file = "AGP_processed.RDS")


# aggregate to higher levels (for comparison to log-contrast regression at
# fixed aggregation levels)
dat <- readRDS("AGP_processed.RDS")
level_names <- c("Phylum",
                 "Class",
                 "Order",
                 "Family",
                 "Genus",
                 "Species",
                 "OTU")
dat_agg <- list()
for (i in seq(1,6)) {
  dat_agg[[level_names[i]]] <- aggregate_to_level(x = dat$x,
                                                  y = dat$y, 
                                                  A = dat$A,
                                                  tax = dat$tax, 
                                                  level = i + 2,
                                                  collapse = TRUE)
}
dat_agg[["OTU"]] <- dat
saveRDS(dat_agg, file = "AGP_aggregated.RDS")