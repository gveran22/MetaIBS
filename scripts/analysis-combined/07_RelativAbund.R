# *********************************************
# Purpose: Plotting relative abundance of phyla
# Date: August 2021
# Author: Salomé Carcy
# *********************************************




# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(phyloseq)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(RColorBrewer)


## 1.2. Data ####
path.root <- "~/Projects/MetaIBS" # CHANGE THIS ROOT DIRECTORY ON YOUR COMPUTER
path.plots <- file.path(path.root, "data/analysis-combined/07_RelativAbund")

path.phylobj    <- file.path(path.root, "data/phyloseq-objects/phyloseq-without-phylotree")
datasets        <- list.files(path.phylobj)
phyloseqobjects <- sapply(datasets, function(x) readRDS(file.path(path.phylobj, x)), USE.NAMES=T, simplify=F)
# names(phyloseqobjects) # sanity check

# Remove Ringel dataset from phyloseqobjects list (we don't have any healthy/IBS information in that dataset)
ringel_name <- names(phyloseqobjects)[grepl("ringel", names(phyloseqobjects))]
phyloseqobjects <- phyloseqobjects[names(phyloseqobjects) != ringel_name]




# ***********************
# 2. PREPROCESS DATA ####
# ***********************

# Merge phyloseq objects
physeq.all <- merge_phyloseq(phyloseqobjects[[1]], phyloseqobjects[[2]]) # Merge first two phyloseq objects in the list
# if there are more than 2 phyloseq objects, merge the rest of them
if(length(phyloseqobjects)>2){
  for (i in 3:length(phyloseqobjects)){
    print(paste0("merging with phyloseq object #", i))
    physeq.all <- merge_phyloseq(physeq.all, phyloseqobjects[[i]])
  }
}
print(physeq.all) # 2,576 samples

# Obtain relative abundances
phylum.table <- physeq.all %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()

# Look at all the phyla & the main ones (ordered by mean abundance per sample)
table(phylum.table$Phylum)
phylum.table %>%
  group_by(Phylum) %>%
  summarize(Mean=mean(Abundance)) %>%
  arrange(-Mean)

# Separate the main phyla (top 5) from the rest
main_phyla <- c("Firmicutes", "Bacteroidota", "Proteobacteria", "Actinobacteriota", "Verrucomicrobiota")
phylum.table.main <- phylum.table %>%
  mutate(Phylum=replace(Phylum, !(Phylum %in% main_phyla), "Other")) %>%
  mutate(Phylum=factor(Phylum, levels=c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota", "Other")))
table(phylum.table.main$Phylum) # sanity check




# *************************************************
# 3. PLOT PHYLA RELATIVE ABUNDANCES PER SAMPLE ####
# *************************************************

# Set the color theme
# rgb(0,1,0, alpha=0.55) # find code for alpha
colors <- paste0(brewer.pal(6, "Dark2"), "80", sep="")
names(colors) <- c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota", "Other")

## 3.1. Plot main phyla by sample (fecal samples) ####
a <- ggplot(phylum.table.main %>% filter(sample_type == "stool"),
       aes(x = reorder(Sample, Sample, function(x) mean(phylum.table.main[Sample == x & Phylum == 'Bacteroidota', 'Abundance'])),
           y = Abundance, fill = Phylum))+
  facet_wrap(~ host_disease, scales = "free_x") +
  geom_bar(stat = "identity", width=1) +
  scale_fill_manual(values=colors, guide=guide_legend(nrow=6))+
  scale_y_continuous(expand = c(0, 0))+ # remove empty space between axis and plot
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15),
        legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "Fecal samples", y = "Relative abundance")

# Save figure
# ggsave(a, file.path(path.plots, "phyla_relabund_fecal_01.jpg"), width=8, height=5) # high
# ggsave(a, file.path(path.plots, "phyla_relabund_fecal_02.jpg"), width=10, height=5) # large


## 3.2. Plot main phyla by sample (sigmoid samples) ####
b <- ggplot(phylum.table.main %>% filter(sample_type == "sigmoid"),
       aes(x = reorder(Sample, Sample, function(x) mean(phylum.table.main[Sample == x & Phylum == 'Bacteroidota', 'Abundance'])),
           y = Abundance, fill = Phylum))+
  facet_wrap(~ host_disease, scales = "free_x") +
  geom_bar(stat = "identity", width=1) +
  scale_fill_manual(values=colors, guide=guide_legend(nrow=6))+
  scale_y_continuous(expand = c(0, 0))+ # remove empty space between axis and plot
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15),
        # legend.position = "None",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "Biopsy samples", y = "Relative abundance")

# Save figure
# ggsave(b, file.path(path.plots, "phyla_relabund_sigmoid_01.jpg"), width=8, height=5) # high
# ggsave(b, file.path(path.plots, "phyla_relabund_sigmoid_02.jpg"), width=10, height=5) # large




# **************************************************
# 4. PLOT PHYLA RELATIVE ABUNDANCES PER DISEASE ####
# **************************************************

# Get a df with abundance of each Phylum per sample
phylum.sample <- phylum.table.main %>%
  # for each sample, sum the abundance of all ASVs belonging to the same phylum
  group_by(Sample, host_disease, host_subtype, Collection, sample_type, author, Phylum) %>%
  dplyr::summarize(Abundance=sum(Abundance)) %>%
  ungroup()
# Sanity check (sum per sample == 1)
# test <- phylum.sample %>%
#   group_by(Sample) %>%
#   summarize(sum_per_sample=sum(Abundance))
# table(test$sum_per_sample)


#___________________________
## 4.1. By host_disease ####

# Get the mean abundance of each phylum per disease status (across samples)
phylum.disease <- phylum.sample %>%
  # keep only stool samples, 1st collection time point
  filter(Collection=="1st") %>% # only 2,450 samples afterwards
  # get mean of each phylum per disease phenotype
  group_by(sample_type, host_disease, Phylum) %>%
  summarize(Mean=mean(Abundance)) %>%
  ungroup() %>%
  mutate(sample_type=factor(sample_type, levels=c("stool", "sigmoid")))
# Sanity check (sum per disease == 1)
# phylum.disease %>%
#   group_by(sample_type, host_disease) %>%
#   summarize(total=sum(Mean))

# Get info on nb of samples per disease
nsamples.disease <- phylum.table.main %>%
  filter(Collection=="1st") %>%
  group_by(sample_type, host_disease) %>%
  summarize(number=paste0("n=", n_distinct(Sample))) %>% # should sum to 2,450 samples
  mutate(sample_type=factor(sample_type, levels=c("stool", "sigmoid")))

# PLOT samples by host_disease
d <- ggplot(phylum.disease,
       aes(x = host_disease, y = Mean, fill = Phylum))+
  facet_wrap(~sample_type)+
  geom_bar(stat = "identity", aes(color=host_disease)) +
  geom_text(data=nsamples.disease, aes(x = host_disease, y=1.03, label=number, fill=NULL))+
  scale_color_manual(values=c("#3182bd", "#de2d26"))+
  scale_fill_manual(values=colors, guide=guide_legend(nrow=6))+
  scale_y_continuous(expand = c(0, 0), limits=c(0,1.05))+ # remove empty space between axis and plot
  theme_cowplot()+
  theme(axis.text.x = element_text(size=10, color="black"),
        legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank(),
        panel.background=element_blank(),
        axis.line.x = element_line(size=0.5, color="black"))+
  labs(x = "", y = "Proportion")
# ggsave(d, file.path(path.plots, "phyla_relabund_1.jpg"), width=7, height=5)


# ___________________________
## 4.2. By host_subtype ####

# Get the mean abundance of each phylum per disease status (across samples)
phylum.subtype <- phylum.sample %>%
  # keep only stool samples, 1st collection time point
  filter(sample_type == "stool" & Collection=="1st") %>% # only 2,051 samples afterwards
  # get mean of each phylum per IBS subtype
  group_by(host_subtype, Phylum) %>%
  summarize(Mean=mean(Abundance)) %>%
  ungroup() %>%
  mutate(host_subtype=recode(host_subtype, "HC"="Healthy"))
# Sanity check (sum per subtype == 1)
# phylum.subtype %>%
#   group_by(host_subtype) %>%
#   summarize(total=sum(Mean))

# Get info on nb of fecal samples per subtype
nsamples.subtype <- phylum.table.main %>%
  filter(sample_type == "stool" & Collection=="1st") %>%
  group_by(host_subtype) %>%
  summarize(number=paste0("n=", n_distinct(Sample))) %>% # should sum to 2,051 samples
  mutate(host_subtype=recode(host_subtype, "HC"="Healthy"))

# PLOT fecal samples by host_subtype
c <- ggplot(phylum.subtype %>% filter(host_subtype != "IBS-unspecified"),
       aes(x = host_subtype, y = Mean, fill = Phylum))+
  geom_bar(stat = "identity", aes(color=host_subtype)) +
  geom_text(data=nsamples.subtype %>% filter(host_subtype != "IBS-unspecified"), aes(x = host_subtype, y=1.03, label=number, fill=NULL))+
  scale_color_manual(values=c("#3182bd", "#de2d26", "#de2d26", "#de2d26", "#de2d26"))+
  scale_fill_manual(values=colors, guide=guide_legend(nrow=6))+
  scale_y_continuous(expand = c(0, 0), limits=c(0,1.05))+ # remove empty space between axis and plot
  theme_cowplot()+
  theme(axis.text.x = element_text(size=10, color="black"),
        legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank(),
        panel.background=element_blank(),
        axis.line.x = element_line(size=0.5, color="black"))+
  labs(x = "", y = "Proportion", title="Fecal samples")
# ggsave(c, file.path(path.plots, "phyla_relabund_2.jpg"), width=6, height=5)


# ___________________________
## 4.3. By collection time point ####

# Get the mean abundance of each phylum per disease status (across samples)
phylum.collection <- phylum.sample %>%
  # keep only stool samples, 1st collection time point
  filter(author %in% c("Pozuelo", "Mars")) %>%
  # get mean of each phylum per IBS subtype
  group_by(host_disease, author, Collection, Phylum) %>%
  summarize(Mean=mean(Abundance)) %>%
  ungroup() %>%
  mutate(Collection=replace(Collection, Collection=="1st", "t0"),
         Collection=replace(Collection, Collection=="2nd" & author=="Mars", "t0+6mo"),
         Collection=replace(Collection, Collection=="2nd" & author=="Pozuelo", "t0+1mo")) %>%
  mutate(author=factor(author, levels=c("Pozuelo", "Mars")),
         author=recode(author, "Mars"="Mars (sigmoid)"),
         author=recode(author, "Pozuelo"="Pozuelo (stool)"))
# Sanity check (sum per subtype == 1)
# phylum.collection %>%
#   group_by(host_disease, author, Collection) %>%
#   summarize(total=sum(Mean))

# Get info on nb of fecal samples per subtype
nsamples.collection <- phylum.table.main %>%
  filter(author %in% c("Pozuelo", "Mars")) %>%
  group_by(host_disease, author, Collection) %>%
  summarize(number=paste0("n=", n_distinct(Sample))) %>% # should sum to 342 samples
  mutate(Collection=replace(Collection, Collection=="1st", "t0"),
         Collection=replace(Collection, Collection=="2nd" & author=="Mars", "t0+6mo"),
         Collection=replace(Collection, Collection=="2nd" & author=="Pozuelo", "t0+1mo")) %>%
  mutate(author=factor(author, levels=c("Pozuelo", "Mars")),
         author=recode(author, "Mars"="Mars (sigmoid)"),
         author=recode(author, "Pozuelo"="Pozuelo (stool)"))

# PLOT fecal samples by host_subtype (will manually change the facet strips on adobe illustrator)
e <- ggplot(phylum.collection,
       aes(x = host_disease, y = Mean, fill = Phylum))+
  facet_wrap(~author+Collection, nrow=1)+
  geom_bar(aes(color=host_disease), stat = "identity") +
  geom_text(data=nsamples.collection, aes(x = host_disease, y=1.03, label=number, fill=NULL))+
  scale_fill_manual(values=colors, guide=guide_legend(nrow=6))+
  scale_color_manual(values=c("#3182bd", "#de2d26"), name="Disease phenotype")+
  scale_y_continuous(expand = c(0, 0), limits=c(0,1.05))+ # remove empty space between axis and plot
  theme_cowplot()+
  theme(#legend.position = "None",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid = element_blank(),
        panel.background=element_blank(),
        axis.line.x = element_line(size=0.5, color="black"))+
  labs(x = "", y = "Proportion", title="")
# ggsave(e, file.path(path.plots, "phyla_relabund_3.jpg"), width=7, height=5)


# ___________________________
## 4.4. Combine plots ####
ggdraw() +
  draw_plot(a, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(b, x = .55, y = .5, width = .45, height = .5) +
  draw_plot(c, x = 0, y = 0, width = .21, height = .45) +
  draw_plot(d, x = 0.27, y = 0, width = 0.23, height = .45) +
  draw_plot(e, x = 0.55, y = 0, width = 0.45, height = .5) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 15,
                  x = c(0, .55, 0, 0.27, 0.55), y = c(1, 1, 0.45, 0.45, 0.45))
# ggsave(file.path(path.plots, "phyla_relabund_all.jpg"), width=15, height=10)




