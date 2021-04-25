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

# MOCK DATA
# Create mock data
df <- data.frame("ASV" = paste("ASV", rep(1:5), sep=""),
                 "Sequence" = c("ATTTGAC", "GCATTAGCTTTA", "TTTGA", "GCTATCG", "AGCTTTA"),
                 "Phylum" = c("Firmicutes", "Bacteroidota", "Firmicutes", "Actinobacteria", "Bacteroidota"),
                 "Class" = c("Clostridia", "Bacteria", "Clostridia", "Proteomachin", "Bacteria")) %>%
  arrange(nchar(Sequence))

# Init
new.df <- data.frame(matrix(ncol=5, nrow=0))
colnames(new.df) <- c(colnames(df), "ASV_new")

allseq <- df$Sequence
i=1

# Identify ASVs with same sequence and give them same ASV reference
while(length(allseq) !=0 ){
  
  seq <- allseq[1]
  
  cat("\n++ Sequence:", seq, "++\n")
  
  # Subset df by sequence
  tempdf <- df %>%
    filter(str_detect(Sequence, seq)) %>%
    mutate(ASV_new = paste("ASV", i, sep=""))
  
  # Remove sequences already identified
  allseq <- allseq[!allseq %in% tempdf$Sequence]
  cat("\n")
  cat("Seq vector:", allseq, "\n\n")
  
  # Keep only 1 sequence (the shortest one)
  tempdf <- tempdf %>%
    filter(Sequence == Sequence[which.min(nchar(Sequence))])
  
  # Add to new.df
  new.df <- rbind(new.df, tempdf)
  print(new.df)
  
  i=i+1
}




#_______________________________________

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