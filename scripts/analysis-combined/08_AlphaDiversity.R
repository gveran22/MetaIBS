##########################
# Purpose: Alpha-diversity
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
library(reshape2)
library(ggpubr)


# Data
path.phy <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/CLUSTER/PhyloTree/input"
physeq.labus    <- readRDS(file.path(path.phy, "physeq_labus.rds"))
physeq.lopresti <- readRDS(file.path(path.phy, "physeq_lopresti.rds"))
# physeq.ringel   <- readRDS(file.path(path.phy, "physeq_ringel.rds"))
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


# Merge
physeq.all <- merge_phyloseq(physeq.labus,
                             physeq.lopresti,
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
author.order <- c('Labus', 'LoPresti', # 454 pyrosequencing
                  'AGP', 'Liu', 'Pozuelo', # Illumina single end
                  'Fukui', 'Hugerth', 'Mars', 'Zhu', 'Zhuang', # Illumina paired end
                  'Nagel', 'Zeber-Lubecka') # Ion Torrent
sample_data(physeq.all)$author <- factor(sample_data(physeq.all)$author, levels=author.order)




###########
# SHANNON #
###########

# Get Shannon values
plt.shannon <- plot_richness(physeq.all, measures="Shannon")
shannon.df <- plt.shannon$data %>%
  select(c("samples", "value", "host_disease", "host_subtype", "sample_type", "Collection", "author")) %>%
  dplyr::rename(shannon=value)


# Do statistics (HC vs IBS)
for(author in author.order){
  # get the shannon values from each dataset
  print(author)
  df <- shannon.df[shannon.df$author == author & shannon.df$Collection=="1st" & shannon.df$sample_type=="stool",]
  # print number of samples in the subset df (as sanity check)
  cat("# HC/IBS samples:", nrow(df[df$host_disease=="Healthy",]), "/", nrow(df[df$host_disease=="IBS",]), "\n")
  # run wilcox test and print the pval
  if(nrow(df) != 0){ # mars will have 0 samples (because all sigmoid)
    pval <- wilcox.test(df[df$host_disease == "IBS","shannon"],
                        df[df$host_disease == "Healthy","shannon"])$p.value
    cat("P-value:", pval)
  }
  cat("\n\n")
}


# Plot
a <- ggplot(shannon.df %>% filter(Collection=="1st" & sample_type=="stool"),
       aes(x = author, y = shannon, fill=host_disease))+
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), aes(color=host_disease), size=0.2)+
  geom_boxplot(position=position_dodge(width=0.75), outlier.shape = NA, width = 0.4, lwd=0.5, alpha=0.2)+
  # add p-values
  geom_signif(y_position = c(6,6,6,6),
              xmin = c(0.8,4.8,5.8,7.8),
              xmax = c(1.2,5.2,6.2,8.2),
              annotation = c("*","*","***","***"), tip_length = 0)+
  # theme and color
  scale_color_manual(values=c("#3182bd", "#de2d26"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#de2d26"))+
  # ylim(c(-5,17))+ # note: we're cutting one point from Fukui (at y=-12)
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, color="black", hjust=1))+
  labs(x = '', y = "Shannon index", fill="", title="Shannon index")



# ++++++++
# PER HOST_SUBTYPE
# ++++++++

# Subset the dataframe to keep only samples where we know the IBS subtype (and their HC counterpart)
shannon.df.subtype <- shannon.df %>% filter(Collection=="1st",
                                            host_subtype!="IBS-unspecified",
                                            author %in% c("Labus", "LoPresti", "Liu", "Pozuelo", "Mars", "Zhuang", "Nagel", "Zeber-Lubecka"),
                                            !(author=="LoPresti" & sample_type=="sigmoid") # remove the (very) few sigmoid samples of LoPresti
                                            ) %>%
  mutate(author=recode(author, "Mars"="Mars (sigmoid)")) %>%
  mutate(host_subtype=replace(host_subtype, host_subtype=="HC", "Healthy"))

# Do statistics (HC vs IBS)
for(author in c("Labus", "LoPresti", "Liu", "Pozuelo", "Mars (sigmoid)", "Zhuang", "Nagel", "Zeber-Lubecka")){
  # get the shannon values from each dataset
  print(author)
  df <- shannon.df.subtype[shannon.df.subtype$author == author,]
  # print number of samples in the subset df (as sanity check)
  cat("# HC/IBS samples:", nrow(df[df$host_disease=="Healthy",]), "/", nrow(df[df$host_disease=="IBS",]), "\n")
  # Wilcox HC vs IBS-C
  pval <- wilcox.test(df[df$host_subtype == "IBS-D","shannon"],
                      df[df$host_subtype == "Healthy","shannon"])$p.value
  cat("P-value HC vs IBS-D:", pval, "\n")
  # Wilcox HC vs IBS-C
  if(author %in% c("Labus", "LoPresti", "Pozuelo", "Mars (sigmoid)", "Zeber-Lubecka")){
    pval <- wilcox.test(df[df$host_subtype == "IBS-C","shannon"],
                        df[df$host_subtype == "Healthy","shannon"])$p.value
    cat("P-value HC vs IBS-C:", pval, "\n")
  }
  # Wilcox HC vs IBS-M
  if(author %in% c("Labus", "LoPresti", "Pozuelo", "Zeber-Lubecka")){
    pval <- wilcox.test(df[df$host_subtype == "IBS-M","shannon"],
                        df[df$host_subtype == "Healthy","shannon"])$p.value
    cat("P-value HC vs IBS-M:", pval, "\n")
  }
  cat("\n\n")
}

# Plot
annotation.df <- data.frame(author=c("Labus", "Pozuelo"),
                            start=c("Healthy", "Healthy"),
                            end=c("IBS-C", "IBS-D"),
                            label=c("*", "**"),
                            y=c(3,5.5))
c <- ggplot(shannon.df.subtype, aes(x = host_subtype, y = shannon))+
  facet_wrap(~factor(author, levels=c("Labus", "LoPresti", "Liu", "Pozuelo", "Mars (sigmoid)", "Zhuang", "Nagel", "Zeber-Lubecka")),
             scales="free_y", ncol=1, strip.position="right")+
  geom_jitter(aes(color=host_subtype), width=0.05, size=0.5)+
  geom_boxplot(aes(fill=host_subtype), lwd=0.5, width=0.4, outlier.shape=NA, alpha=0.2)+
  geom_signif(data=annotation.df, aes(y_position=y, xmin=start, xmax=end, annotations=label), manual=T, tip_length=0, vjust=0.4)+
  scale_color_manual(values=c("#3182bd", "#a50f15", "#fcae91", "#fb6a4a"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#a50f15", "#fcae91", "#fb6a4a"), guide="none")+
  # ylim(c(0.8,6.4))+
  # THEME
  theme_cowplot()+
  theme(strip.text.y = element_text(angle=0, hjust=0),
        strip.background = element_blank())+
  labs(title="IBS subtypes", x = '', y = "Shannon index")




###########
# SIMPSON #
###########

# Get simpson values
plt.simpson <- plot_richness(physeq.all, measures="Simpson")
simpson.df <- plt.simpson$data %>%
  select(c("samples", "value", "host_disease", "host_subtype", "sample_type", "Collection", "author")) %>%
  rename(simpson=value)


# Do statistics
for(author in author.order){
  # get the shannon values from each dataset
  print(author)
  df <- simpson.df[simpson.df$author == author & simpson.df$Collection=="1st" & simpson.df$sample_type=="stool",]
  # print number of samples in the subset df (as sanity check)
  cat("# HC/IBS samples:", nrow(df[df$host_disease=="Healthy",]), "/", nrow(df[df$host_disease=="IBS",]), "\n")
  # run wilcox test and print the pval
  if(nrow(df) != 0){
    pval <- wilcox.test(df[df$host_disease == "IBS","simpson"],
                        df[df$host_disease == "Healthy","simpson"])$p.value
    cat("P-value:", pval)
  }
  cat("\n\n")
}


# Plot
b <- ggplot(simpson.df %>% filter(Collection=="1st" & sample_type=="stool"),
       aes(x = author, y = simpson, fill=host_disease))+
  geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), aes(color=host_disease), size=0.2)+
  geom_boxplot(position=position_dodge(width=0.75), outlier.shape = NA, width = 0.4, lwd=0.5, alpha=0.2)+
  # add p-values
  geom_signif(y_position = c(1,1,1),
              xmin = c(4.8,5.8,7.8),
              xmax = c(5.2,6.2,8.2),
              annotation = c("*","***","**"), tip_length = 0)+
  # theme and color
  scale_color_manual(values=c("#3182bd", "#de2d26"), guide="none")+
  scale_fill_manual(values=c("#3182bd", "#de2d26"))+
  theme_cowplot()+
  ylim(c(0.5,1))+ # note: we're cutting 17 points (mostly from AGP)
  theme(axis.text.x = element_text(angle = 45, color="black", hjust=1))+
  labs(x = '', y = "Simpson index", fill="", title="Simpson index")




#####################
# PLOT ALL TOGETHER #
#####################

ggdraw() +
  draw_plot(a, x = 0, y = .5, width = .65, height = .5) +
  draw_plot(b, x = 0, y = 0, width = .65, height = .5) +
  draw_plot(c, x = 0.65, y = 0, width = 0.35, height = 1) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0, 0.65), y = c(1, 0.5, 1))
ggsave("~/Projects/IBS_Meta-analysis_16S/data/plots_paper/alpha-diversity.jpg", width=15, height=10)



###############
# AGP SUBTYPE #
###############

# Infer IBS subtype in AGP data
# sample_data(physeq.agp)[sample_data(physeq.agp)$host_disease=="IBS" & sample_data(physeq.agp)$bowel_movement_quality=="Constipated", "host_subtype"] <- "IBS-C"
# sample_data(physeq.agp)[sample_data(physeq.agp)$host_disease=="IBS" & sample_data(physeq.agp)$bowel_movement_quality=="Diarrhea", "host_subtype"] <- "IBS-D"
# sample_data(physeq.agp)[sample_data(physeq.agp)$host_disease=="Healthy" & sample_data(physeq.agp)$bowel_movement_quality!="Normal", "host_subtype"] <- "HC-unknown"

# Plot
# plot_richness(subset_samples(physeq.agp, host_subtype != "HC-unknown" & host_subtype != "IBS-unspecified"),
#               x="host_subtype", measures=c("Shannon", "Simpson")) +
#   geom_boxplot(fill=NA, width=0.3) +
#   theme_bw() +
#   labs(x="", y="", title="AGP")
# ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/08_AlphaDiversity/agp_subtype.jpg", width=4, height=4)