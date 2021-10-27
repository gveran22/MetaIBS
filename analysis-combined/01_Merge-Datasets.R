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
library(reshape2)
library(ggupset)

# Data
path <- "~/Projects/IBS_Meta-analysis_16S"
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




##############################################
# MERGE PHYLOSEQ OBJECTS WITH MERGE_PHYLOSEQ #
##############################################

# Merge phyloseq objects
physeq.all <- merge_phyloseq(physeq.labus,
                             physeq.lopresti,
                             physeq.ringel,
                             physeq.agp,
                             physeq.liu,
                             physeq.pozuelo,
                             physeq.fukui,
                             physeq.mars,
                             physeq.hugerth,
                             physeq.zhu,
                             physeq.zhuang,
                             physeq.nagel,
                             physeq.zeber)

# Compare the number of ASVs before/after merging
sum(ntaxa(physeq.labus)+
    ntaxa(physeq.lopresti)+
    ntaxa(physeq.ringel)+
    ntaxa(physeq.agp)+
    ntaxa(physeq.liu)+
    ntaxa(physeq.pozuelo)+
    ntaxa(physeq.fukui)+
    ntaxa(physeq.mars)+
    ntaxa(physeq.hugerth)+
    ntaxa(physeq.zhu)+
    ntaxa(physeq.zhuang)+
    ntaxa(physeq.nagel)+
    ntaxa(physeq.zeber)) # 81,452 before

ntaxa(physeq.all) # 79,918 after
min(nchar(taxa_names(physeq.all)))

# Put datasets in a list
datasets <- list("labus"    = physeq.labus,
                 "lopresti" = physeq.lopresti,
                 "ringel"   = physeq.ringel,
                 "agp"      = physeq.agp,
                 "liu"      = physeq.liu,
                 "pozuelo"  = physeq.pozuelo,
                 "fukui"    = physeq.fukui,
                 "mars"     = physeq.mars,
                 "hugerth"  = physeq.hugerth,
                 "zhu"      = physeq.zhu,
                 "zhuang"   = physeq.zhuang,
                 "nagel"    = physeq.nagel,
                 "zeber"    = physeq.zeber)


# Let's get all ASVs in a dataframe and check if we can find common ones
asv.df <- melt(lapply(datasets, function(x) taxa_names(x)))
colnames(asv.df) <- c("asv", "author")
length(unique(asv.df$asv)) # we do find 79,918 unique sequences

# Let's see which datasets share the exact same ASV sequence
common.asv <- asv.df %>%
  group_by(asv) %>%
  summarize(Datasets = list(unique(author))) %>%
  filter(lengths(Datasets)>1)

# jpeg("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/01_Merge-Datasets/commonASV_merge-phyloseq-funct.jpg", width=2000, height=2000, res=400)
ggplot(common.asv, aes(x=Datasets))+
  geom_bar() +
  scale_x_upset()+
  labs(y="# of common ASVs")
# dev.off()






##############################
# MERGE COMMON ASVs MANUALLY #
##############################

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
  # fill(Genus, .direction="downup") %>%
  filter(n_distinct(Genus, na.rm=TRUE)>1) %>%
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








