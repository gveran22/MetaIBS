##########################
# Purpose: Exporting ASV sequences in fasta format for annotIEM
# Date: November 2021
# Author: Salomé Carcy
##########################



##########
# IMPORT #
##########

# Libraries
library(Biostrings)


# Data
path <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/Taxonomy/input"
seqtable.agp      <- readRDS(file.path(path, "seqtable_agp.rds"))
seqtable.fukui    <- readRDS(file.path(path, "seqtable_fukui.rds"))
seqtable.hugerth  <- readRDS(file.path(path, "seqtable_hugerth.rds"))
seqtable.labus    <- readRDS(file.path(path, "seqtable_labus.rds"))
seqtable.liu      <- readRDS(file.path(path, "seqtable_liu.rds"))
seqtable.lopresti <- readRDS(file.path(path, "seqtable_lopresti.rds"))
seqtable.mars     <- readRDS(file.path(path, "seqtable_mars.rds"))
seqtable.nagel    <- readRDS(file.path(path, "seqtable_nagel.rds"))
seqtable.pozuelo  <- readRDS(file.path(path, "seqtable_pozuelo.rds"))
seqtable.ringel   <- readRDS(file.path(path, "seqtable_ringel.rds"))
seqtable.zeber    <- readRDS(file.path(path, "seqtable_zeber.rds"))
seqtable.zhu      <- readRDS(file.path(path, "seqtable_zhu.rds"))
seqtable.zhu.small<- readRDS(file.path(path, "seqtable_zhu_cutFWDprimer.rds"))
seqtable.zhuang   <- readRDS(file.path(path, "seqtable_zhuang.rds"))

# Get a list
datasets <- list("agp"              = seqtable.agp,
                 "fukui"            = seqtable.fukui,
                 "hugerth"          = seqtable.hugerth,
                 "labus"            = seqtable.labus,
                 "liu"              = seqtable.liu,
                 "lopresti"         = seqtable.lopresti,
                 "mars"             = seqtable.mars,
                 "nagel"            = seqtable.nagel,
                 "pozuelo"          = seqtable.pozuelo,
                 "ringel"           = seqtable.ringel,
                 "zeber"            = seqtable.zeber,
                 "zhu"              = seqtable.zhu,
                 "zhu_cutFWDprimer" = seqtable.zhu.small,
                 "zhuang"           = seqtable.zhuang)




###################################
# EXPORT TO FASTA FORMAT THE ASVs #
###################################

for (author in names(datasets)){
  print(author)
  print(dim(datasets[[author]]))
  
  # Store ASVs sequences as DNAStringSet
  asv <- DNAStringSet(colnames(datasets[[author]]))
  
  # Save ASVs sequences in FASTA file
  file.name <- paste("~/Projects/IBS_Meta-analysis_16S/asv-sequences/asv_", author, ".fasta", sep="")
  writeXStringSet(x=asv, filepath=file.name, compress=FALSE, format="fasta")
  
}











