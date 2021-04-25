##########################
# Purpose: merge phyloseq object and merge ASVs
# Date: April 2021
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(tidyverse)
library(phyloseq)

# Data
path <- "~/Projects/IBS_Meta-analysis_16S"
physeq.fukui <- readRDS(file.path(path, "phyloseq-objects/physeq_fukui.rds"))
physeq.hugerth <- readRDS(file.path(path, "phyloseq-objects/physeq_hugerth.rds"))
physeq.labus <- readRDS(file.path(path, "phyloseq-objects/physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path, "phyloseq-objects/physeq_lopresti.rds"))
physeq.nagel <- readRDS(file.path(path, "phyloseq-objects/physeq_nagel.rds"))
physeq.pozuelo <- readRDS(file.path(path, "phyloseq-objects/physeq_pozuelo.rds"))
physeq.zeber <- readRDS(file.path(path, "phyloseq-objects/physeq_zeber.rds"))
physeq.zhu <- readRDS(file.path(path, "phyloseq-objects/physeq_zhu.rds"))
physeq.zhuang <- readRDS(file.path(path, "phyloseq-objects/physeq_zhuang.rds"))

datasets <- list("fukui" = physeq.fukui,
                 "hugerth" = physeq.hugerth,
                 "labus" = physeq.labus,
                 "lopresti" = physeq.lopresti,
                 "nagel" = physeq.nagel,
                 "pozuelo" = physeq.pozuelo,
                 "zeber" = physeq.zeber,
                 "zhu" = physeq.zhu,
                 "zhuang" = physeq.zhuang)

#####################
# CHECK COMMON ASVs #
#####################

# Get ASV sequences for each dataset
ASVs <- vector("list", length(datasets))
names(ASVs) <- names(datasets)

i=0
for (physeq.cohort in datasets){
  i=i+1
  taxa_names(physeq.cohort) <- refseq(physeq.cohort)
  ASVs[[i]] <- colnames(otu_table(physeq.cohort))
}


# Quick peek if we can find Nagel's ASVs in Pozuelo's ASVs (V4 variable region)
mean(nchar(ASVs$nagel)) # ~254 bp per ASV
mean(nchar(ASVs$pozuelo)) # ~246 bp per ASV

for (asv in ASVs$pozuelo){
  
  print(table(str_detect(ASVs$nagel, asv)))
  
}


# Previous code
for (i in 1:420){
  table(substr(taxa_names(physeq.zhu), start = i, stop = i+248) %in% taxa_names(physeq.pozuelo))
}

table(substr(taxa_names(physeq.zhu), start = 160, stop = 413) %in% taxa_names(physeq.pozuelo))

# 9 OTUs of length 246bp are common
# pozuelo: 200 - 248 bp