##########################
# Purpose: UMAP plotting
# Date: April 2022
# Author: Salomé Carcy
##########################




##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(ggplot2)
library(cowplot)
library(umap)
library(tidyverse)


# Phyloseq objects
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
physeq.labus    <- readRDS(file.path(path.phy, "physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path.phy, "physeq_lopresti.rds"))
physeq.ringel   <- readRDS(file.path(path.phy, "physeq_ringel.rds"))
physeq.agp      <- readRDS(file.path(path.phy, "physeq_agp.rds"))
physeq.liu      <- readRDS(file.path(path.phy, "physeq_liu.rds"))
physeq.pozuelo  <- readRDS(file.path(path.phy, "physeq_pozuelo.rds"))
physeq.fukui    <- readRDS(file.path(path.phy, "physeq_fukui.rds"))
physeq.hugerth  <- readRDS(file.path(path.phy, "physeq_hugerth.rds"))
physeq.mars     <- readRDS(file.path(path.phy, "physeq_mars.rds"))
physeq.zhu      <- readRDS(file.path(path.phy, "physeq_zhu.rds"))
physeq.zhuang   <- readRDS(file.path(path.phy, "physeq_zhuang.rds"))
physeq.nagel    <- readRDS(file.path(path.phy, "physeq_nagel.rds"))
physeq.zeber    <- readRDS(file.path(path.phy, "physeq_zeber.rds"))

physeq.all <- merge_phyloseq(physeq.labus,
                             physeq.lopresti,
                             physeq.ringel,
                             physeq.agp,
                             physeq.liu,
                             physeq.pozuelo,
                             physeq.fukui,
                             physeq.hugerth,
                             physeq.mars,
                             physeq.zhu,
                             physeq.zhuang,
                             physeq.nagel,
                             physeq.zeber)
cat("Nb of samples:", nsamples(physeq.all), "\n") # should be 2651
cat("Nb of taxa:", ntaxa(physeq.all), "\n") # should be 79,943

covariates <- data.frame(disease         = sample_data(physeq.all)[,'host_disease'],
                         subtype         = sample_data(physeq.all)[,'host_subtype'],
                         sample_type     = sample_data(physeq.all)[,'sample_type'],
                         seq_tech        = sample_data(physeq.all)[,'sequencing_tech'],
                         author          = sample_data(physeq.all)[,'author'],
                         variable_region = sample_data(physeq.all)[,'variable_region'],
                         age             = sample_data(physeq.all)[,'host_age'],
                         bmi             = sample_data(physeq.all)[,'host_bmi'])
covariates <- covariates %>%
  mutate(author = factor(covariates$author, levels = c('Labus', 'LoPresti', 'Ringel', # 454 pyrosequencing
                                                       'AGP', 'Liu', 'Pozuelo', # Illumina single end
                                                       'Fukui', 'Hugerth', 'Mars', 'Zhu', 'Zhuang', # Illumina paired end
                                                       'Nagel', 'Zeber-Lubecka'))) %>%
  mutate(sequencing_tech = factor(covariates$sequencing_tech, levels= c("454 pyrosequencing", "Illumina single-end",
                                                                        "Illumina paired-end", "Ion Torrent")))


# Log-ratios between families (computed on the cluster)
path <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/02_UMAP"
ratiosFamily <- readRDS(file.path(path, "pseudocounts_aft-agg/ratiosFamily.rds"))
# Sanity checks
ratiosFamily[1:10,1:10]
dim(ratiosFamily) # should have 2651 samples (rows) and 64,620 predictors (family log ratios)




#################
# STOOL SAMPLES #
#################

# *********
# RUN UMAP 
# *********

# Subset samples to only stool samples
stool.samples <- sample_names(subset_samples(physeq.all, sample_type=="stool"))
# head(stool.samples)
length(stool.samples) # 2220 stool samples
ratiosFamily.stool <- ratiosFamily[rownames(ratiosFamily) %in% stool.samples,]
dim(ratiosFamily.stool)


# Run UMAP
cat("Nb of samples", nrow(ratiosFamily.stool), "and nb of predictors", ncol(ratiosFamily.stool), "\n")
set.seed(123)
umap <- uwot::umap(ratiosFamily.stool, # umap on samples (rows) and taxa ratios (columns)
                   n_neighbors=20, n_components=3, scale=F, n_threads=8)

# Get the (x,y) coordinates from the UMAP
dims.umap <- umap %>% as.data.frame()
colnames(dims.umap) <- c("UMAP_1", "UMAP_2", "UMAP_3")
# table(rownames(dims.umap) == rownames(ratiosFamily.stool)) # sanity check

# Add covariates
# setdiff(rownames(dims.umap), rownames(covariates)) # sanity check
dims.umap <- merge(as.data.frame(dims.umap), covariates, by="row.names") # umap+covariates
colnames(dims.umap)[1] <- "Sample"
cat("Nb of samples", nrow(dims.umap), "\n")


# **********
# PLOT UMAP
# **********

# Per author/dataset
a <- ggplot(dims.umap,
            aes(x = UMAP_1, y = UMAP_2, fill = author))+
  geom_point(size = 2, pch=21, colour="black", stroke=.1)+
  scale_fill_manual(values=pals::brewer.paired(12), name="")+ # author
  # xlim(c(-2.5,2.5))+
  # ylim(c(-2.5,3.5))+
  labs(x="UMAP 1", y="UMAP 2", title = "Dataset")+
  theme_cowplot()+
  theme(line = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Per host_disease
b <- ggplot(dims.umap,
       aes(x = UMAP_1, y = UMAP_2, color = host_disease))+
  geom_point(size = .5)+
  scale_color_manual(values = c('blue', 'red', 'black'), name="")+ # disease
  # xlim(c(-2.5,2.5))+
  # ylim(c(-2.5,3.5))+
  labs(x="UMAP 1", y="UMAP 2", title = "Disease phenotype")+
  theme_cowplot()+
  theme(line = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Per seqtech
c <- ggplot(dims.umap,
       aes(x = UMAP_1, y = UMAP_2, color = sequencing_tech))+
  geom_point(size = .5)+
  scale_color_manual(values=c("#6a51a3", "#a1d99b", "#238b45", "#f16913"), name="")+ # seqtech
  # xlim(c(-2.5,2.5))+
  # ylim(c(-2.5,3.5))+
  labs(x="UMAP 1", y="UMAP 2", title = "Seq. technology")+
  theme_cowplot()+
  theme(line = element_blank(),
        plot.title = element_text(hjust = 0.5))


# Combine three plots
ggdraw() +
  draw_plot(a, x = .25,  y = .4, width = .6, height = .6) +
  draw_plot(b, x = .1,   y = 0,  width = .4, height = .4) +
  draw_plot(c, x = .5,   y = 0,  width = .4, height = .4)
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/02_UMAP/plots/umap_family_legends.jpeg", width=10, height=10)




###################
# SIGMOID SAMPLES #
###################

# *********
# RUN UMAP 
# *********

# Subset samples to only sigmoid biopsy samples
sigm.samples <- sample_names(subset_samples(physeq.all, sample_type=="sigmoid"))
# head(sigm.samples)
length(sigm.samples) # 431 sigmoid samples
ratiosFamily.sigm <- ratiosFamily[rownames(ratiosFamily) %in% sigm.samples,]
dim(ratiosFamily.sigm)


# Run UMAP
cat("Nb of samples", nrow(ratiosFamily.sigm), "and nb of predictors", ncol(ratiosFamily.sigm), "\n")
set.seed(123)
umap.sigm <- uwot::umap(ratiosFamily.sigm, # umap on samples (rows) and taxa ratios (columns)
                        n_neighbors=20, n_components=3, scale=F, n_threads=8)

# Get the (x,y) coordinates from the UMAP
dims.umap.sigm <- umap.sigm %>% as.data.frame()
colnames(dims.umap.sigm) <- c("UMAP_1", "UMAP_2", "UMAP_3")
# table(rownames(dims.umap.sigm) == rownames(ratiosFamily.sigm)) # sanity check

# Add covariates
# setdiff(rownames(dims.umap.sigm), rownames(covariates)) # sanity check
dims.umap.sigm <- merge(as.data.frame(dims.umap.sigm), covariates, by="row.names") # umap+covariates
colnames(dims.umap.sigm)[1] <- "Sample"
cat("Nb of samples", nrow(dims.umap.sigm), "\n")


# **********
# PLOT UMAP
# **********

# Per author/dataset
d <- ggplot(dims.umap.sigm,
            aes(x = UMAP_1, y = UMAP_2, fill = author))+
  geom_point(size = 2, pch=21, colour="black", stroke=.1)+
  scale_fill_manual(values=c("#1F78B4", "#FF7F00", "#dfc27d"), name="")+ # author
  labs(x="UMAP 1", y="UMAP 2", title = "Dataset")+
  theme_cowplot()+
  theme(line = element_blank(),
        # legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# Per host_disease
e <- ggplot(dims.umap.sigm,
            aes(x = UMAP_1, y = UMAP_2, color = host_disease))+
  geom_point(size = 2)+
  scale_color_manual(values = c('blue', 'red'), name="")+ # disease
  labs(x="UMAP 1", y="UMAP 2", title = "Disease phenotype")+
  theme_cowplot()+
  theme(line = element_blank(),
        # legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# Per seqtech
f <- ggplot(dims.umap.sigm,
            aes(x = UMAP_1, y = UMAP_2, color = sequencing_tech))+
  geom_point(size = 2)+
  scale_color_manual(values=c("#6a51a3", "#238b45"), name="")+ # seqtech
  labs(x="UMAP 1", y="UMAP 2", title = "Seq. technology")+
  theme_cowplot()+
  theme(line = element_blank(),
        # legend.position = "none",
        plot.title = element_text(hjust = 0.5))


# Combine three plots
ggdraw() +
  draw_plot(d, x = 0,     y = 0, width = .32, height = 1) +
  draw_plot(e, x = .33,   y = 0, width = .32, height = 1) +
  draw_plot(f, x = .66,   y = 0, width = .32, height = 1)
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/02_UMAP/plots/sigmoid_umap_family.jpeg", width=15, height=5)




################
# HOST SUBTYPE #
################

# Stool
g <- ggplot(dims.umap %>% arrange(match(host_subtype, c("IBS-unspecified", "NA", "HC", "IBS-C", "IBS-D", "IBS-M"))),
       aes(x = UMAP_1, y = UMAP_2, color = host_subtype))+
  geom_point(size = 2)+
  scale_color_manual(values = c("#3182bd", "#a50f15", "#fcae91", "#fb6a4a", "#f0f0f0", "#d9d9d9"), name="")+ # disease
  labs(x="UMAP 1", y="UMAP 2", title = "Stool samples")+
  theme_cowplot()+
  theme(line = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# Sigmoid
h <- ggplot(dims.umap.sigm %>% arrange(match(host_subtype, c("IBS-unspecified", "HC", "IBS-C", "IBS-D", "IBS-M"))),
       aes(x = UMAP_1, y = UMAP_2, color = host_subtype))+
  geom_point(size = 2)+
  scale_color_manual(values = c("#3182bd", "#a50f15", "#fcae91", "#fb6a4a", "#f0f0f0"), name="")+ # disease
  labs(x="UMAP 1", y="UMAP 2", title = "Sigmoid samples")+
  theme_cowplot()+
  theme(line = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# Combine two plots
ggdraw() +
  draw_plot(g, x = 0,  y = 0, width = .5, height = 1) +
  draw_plot(h, x = .5, y = 0, width = .5, height = 1)
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/02_UMAP/plots/host-subtype_umap_family.jpeg", width=11, height=6)

