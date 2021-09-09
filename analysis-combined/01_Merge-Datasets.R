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
# MERGE COMMON ASVs #
#####################

#_______________________________________
# MOCK DATA
# Create mock data
df <- data.frame("author" = c("fukui", "fukui", "zhu", "pozuelo", "fukui", "pozuelo"),
                 "ASV" = paste("ASV", rep(1:6), sep=""),
                 "Sequence" = c("ATTTGAC", "GCATTAGCTTTA", "TTTGA", "GCTATCG", "AGCTTTA", "TTTGAT"),
                 "Phylum" = c("Firmicutes", "Bacteroidota", "Firmicutes", "Actinobacteria", "Bacteroidota", "Proteobacteria"),
                 "Class" = c("Clostridia", "Bacteria", NA, "Proteomachin", "Bacteria", "Gammaproteo")) %>%
  arrange(nchar(Sequence))


# Init
new.df <- data.frame(matrix(ncol=6, nrow=0))
colnames(new.df) <- c(colnames(df), "ASV_new")

allseq <- df$Sequence
i=1

# Identify ASVs with same sequence and give them same ASV reference
while(length(allseq) !=0 ){
  
  seq <- allseq[1]
  
  cat("\n++ Sequence:", seq, "++\n")
  
  # Subset df by sequence
  tempdf <- df %>%
    filter(str_detect(string=Sequence, pattern=seq)) %>%
    filter(n_distinct(Phylum) == 1) %>%
    mutate(ASV_new = paste("ASV", i, sep=""))
  
  # Add to new.df
  new.df <- rbind(new.df, tempdf)
  print(new.df)
  
  # Remove sequences already identified
  allseq <- allseq[!allseq %in% tempdf$Sequence]
  cat("\n")
  cat("Seq vector:", allseq, "\n\n")
  
  i=i+1
}


# Keep only 1 sequence per ASV (the shortest one)
new.df %>%
  group_by(ASV_new) %>%
  filter(Sequence == Sequence[which.min(nchar(Sequence))]) %>%
  select(-ASV)

# See common ASVs between datasets (authors)
new.df %>%
  select(author, ASV_new)


#_______________________________________

# BUILD TAXONOMIC TABLE WITH ALL DATASSETS

taxtable <- lapply(datasets, function(physeq){
  
  taxa_names(physeq) <- refseq(physeq)
  tax_table(physeq) %>%
    as.data.frame() %>%
    mutate(author = unique(sample_data(physeq)$author)) %>%
    mutate(vregion = unique(sample_data(physeq)$variable_region)) %>%
    mutate(sequence = rownames(tax_table(physeq)))
  
}) %>%
  bind_rows() %>%
  arrange(nchar(sequence))

print(nrow(taxtable)) # 49,619 ASVs


# MERGE COMMON ASVs

taxtable.new <- data.frame(matrix(ncol=ncol(taxtable)+1, nrow=0))
colnames(taxtable.new) <- c(colnames(taxtable), "ASV")

allseq <- taxtable$sequence
i=1

# Identify ASVs with same sequence and give them same ASV reference
while(length(allseq) !=0 ){
  
  seq <- allseq[1]
  
  # Subset df by sequence
  tempdf <- taxtable %>%
    filter(str_detect(string=sequence, pattern=seq)) %>%
    mutate(ASV = paste("ASV", i, sep=""))
  
  # Add to new.df
  taxtable.new <- rbind(taxtable.new, tempdf)

  # Remove sequences already identified
  allseq <- allseq[!allseq %in% tempdf$sequence]
  
  # Follow progress of computation
  if(length(allseq)%%100==0){ print(length(allseq)) }
  
  i=i+1
}

# saveRDS(taxtable.new, "~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/01_Merge-Datasets/taxtable_merged.rds")

# Sanity check
taxtable.new <- readRDS(file.path(path, "data/analysis-combined/01_Merge-Datasets/taxtable_merged.rds"))
# check that the shortest sequence of each ASV can be found in the other sequences
taxtable.new %>%
  group_by(ASV) %>%
  mutate(verif = str_detect(string=sequence, pattern=sequence[which.min(nchar(sequence))])) %>%
  ungroup() %>%
  summarise(n_true = sum(verif), n_false=sum(!verif))
# Only TRUEs!! =)

# check that for each ASV, the "middle-length" sequences can be found in the bigger length sequences
test <- taxtable.new %>%
  group_by(ASV) %>%
  filter(sequence != sequence[which.min(nchar(sequence))]) %>% # remove shortest sequence
  mutate(verif = str_detect(string=sequence, pattern=sequence[which.min(nchar(sequence))]),
         length = nchar(sequence)) %>%
  ungroup() %>%
  summarise(n_true = sum(verif), n_false=sum(!verif))


# check that ASVs with several sequences have the same taxa assigned to them
test <- taxtable.new %>%
  group_by(ASV) %>%
  # if several sequences for the ASV, replace NA values by the known Genus of other sequences
  fill(Genus, .direction="downup") %>%
  filter(n_distinct(Genus)>1) %>%
  arrange(nchar(sequence)) %>%
  arrange(ASV) %>%
  mutate(seqlength=nchar(sequence))

# Keep only 1 sequence per ASV (the shortest one)
taxtable.new %>%
  group_by(ASV) %>%
  filter(sequence == sequence[which.min(nchar(sequence))])


# See common ASVs between datasets (authors)
test <- taxtable.new %>%
  select(author, ASV)








