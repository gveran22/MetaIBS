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
ratiosFamily[1:5,1:5]
dim(ratiosFamily) # should have 2651 samples (rows) and 64,620 predictors (family log ratios)




############
# RUN UMAP #
############

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
# table(rownames(dims.umap) == rownames(ratiosFamily)) # sanity check

# Add covariates
# setdiff(rownames(dims.umap), rownames(covariates)) # sanity check
dims.umap <- merge(as.data.frame(dims.umap), covariates, by="row.names") # umap+covariates
colnames(dims.umap)[1] <- "Sample"
cat("Nb of samples", nrow(dims.umap), "\n")



#############
# PLOT UMAP #
#############

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
  labs(x="UMAP 1", y="UMAP 2", title = " Disease phenotype")+
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




