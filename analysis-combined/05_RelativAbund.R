##########################
# Purpose: Plotting relative abundance of phyla
# Date: August 2021
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(RColorBrewer)

# Data
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
# physeq.ringel <- readRDS(file.path(path.phy, "physeq_ringel.rds"))
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




###################
# PREPROCESS DATA #
###################

# Merge phyloseq objects
physeq <- merge_phyloseq(physeq.labus,
                         physeq.lopresti,
                         physeq.pozuelo,
                         physeq.zhuang,
                         physeq.zhu,
                         physeq.hugerth,
                         physeq.fukui,
                         physeq.mars,
                         physeq.liu,
                         physeq.agp,
                         physeq.nagel,
                         physeq.zeber)

# Obtain relative abundances
phylum.table <- physeq %>%
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


#############################################
# PLOT PHYLA RELATIVE ABUNDANCES PER SAMPLE #
#############################################

# Set the color theme
# rgb(0,1,0, alpha=0.55) # find code for alpha
colors <- paste0(brewer.pal(6, "Dark2"), "80", sep="")
names(colors) <- c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota", "Other")

# Plot main phyla by sample (fecal samples)
ggplot(phylum.table.main %>% filter(sample_type == "stool"),
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
        #legend.position = "None",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "Fecal samples", y = "Relative abundance")

# Save figure
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/05_Relative-Abund/phyla_relabund_fecal_02.jpg", width=8, height=5) # high
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/05_Relative-Abund/phyla_relabund_fecal_03.jpg", width=10, height=5) # large


# Plot main phyla by sample (sigmoid samples)
ggplot(phylum.table.main %>% filter(sample_type == "sigmoid"),
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
        #legend.position = "None",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "Biopsy samples", y = "Relative abundance")

# Save figure
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/05_Relative-Abund/phyla_relabund_sigmoid_02.jpg", width=8, height=5) # high
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/05_Relative-Abund/phyla_relabund_sigmoid_03.jpg", width=10, height=5) # large




##############################################
# PLOT PHYLA RELATIVE ABUNDANCES PER DISEASE #
##############################################

# *********
# BY HOST_DISEASE
# *********
# Get a df with average abundance of each Phylum per sample
test <- phylum.table.main %>%
  # keep only stool samples, 1st collection time point
  filter(sample_type == "stool" & Collection=="1st") %>%
  select(Sample, Abundance, Phylum, host_disease, host_subtype) %>%
  # for each sample, sum the abundance of all ASVs belonging to the same phylum
  group_by(Sample, Phylum) %>%
  mutate(abundance=sum(Abundance)) %>%
  ungroup() %>%
  group_by(Sample) %>%
  summarize(sum(abundance))
  # get the mean abundance of each phylum per disease status (across samples)
  group_by(host_disease, Phylum) %>%
  summarize(Mean=mean(abundance))

  # Plot fecal samples by host_disease
ggplot(test,
       aes(x = host_disease, y = Mean, fill = Phylum))+
  geom_bar(stat = "identity", position='fill') +
  scale_fill_manual(values=colors, guide=guide_legend(nrow=6))+
  scale_y_continuous(expand = c(0, 0))+ # remove empty space between axis and plot
  theme_cowplot()+
  theme(axis.text.x = element_text(size=10, color="black"),
        #legend.position = "None",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.2, 'cm'),
        panel.grid = element_blank(),
        panel.background=element_blank(),
        axis.line.x = element_line(size=0.5, color="black"))+
  labs(x = "", y = "Proportion")


test <- phylum.table.main %>%
  filter(sample_type == "stool" & Collection=="1st") %>%
  group_by(host_subtype) %>%
  summarize(number=paste0("n=", n_distinct(Sample)))

ggplot(phylum.table.main %>% filter(sample_type == "stool" & Collection=="1st"),
       aes(x = host_subtype, y = Abundance, fill = Phylum))+
  geom_bar(stat = "identity", position = "fill") +
  geom_text(data=test, aes(x = host_subtype, y=1.03, label=number, fill=NULL))+
  scale_fill_manual(values=colors, guide=guide_legend(nrow=6))+
  scale_y_continuous(expand = c(0, 0), limits=c(0,1.1))+ # remove empty space between axis and plot
  theme_cowplot()+
  theme(axis.text.x = element_text(size=10, color="black"),
        #legend.position = "None",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.2, 'cm'),
        panel.grid = element_blank(),
        panel.background=element_blank(),
        axis.line.x = element_line(size=0.5, color="black"))+
  labs(x = "", y = "Proportion")

# Save figure
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/05_Relative-Abund/phyla_relabund_2.jpg", width=5, height=5)
