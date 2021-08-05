##########################
# Purpose: Plotting seq depth pre/post DADA2
# Date: August 2021
# Author: Salomé Carcy
##########################


##########
# IMPORT #
##########

# Libraries
library(ggplot2)
library(RColorBrewer)

# Data
path <- "~/Projects/IBS_Meta-analysis_16S/data/analysis-individual/DADA2-FILT"
df <- melt(rbind(readRDS(file.path(path, "seqdepth_agp.rds")),
                  readRDS(file.path(path, "seqdepth_fukui.rds")),
                  readRDS(file.path(path, "seqdepth_hugerth.rds")),
                  readRDS(file.path(path, "seqdepth_labus.rds")),
                  readRDS(file.path(path, "seqdepth_liu.rds")),
                  readRDS(file.path(path, "seqdepth_lopresti.rds")),
                  readRDS(file.path(path, "seqdepth_mars.rds")),
                  readRDS(file.path(path, "seqdepth_nagel.rds")),
                  readRDS(file.path(path, "seqdepth_pozuelo.rds")),
                  #readRDS(file.path(path, "seqdepth_ringel.rds")),
                  readRDS(file.path(path, "seqdepth_zeber.rds")),
                  readRDS(file.path(path, "seqdepth_zhu.rds")),
                  readRDS(file.path(path, "seqdepth_zhuang.rds"))))


########
# PLOT #
########

# Set order
authors_order <- c('Labus', 'LoPresti',
                   'Pozuelo', 'Zhuang', 'Zhu', 'Hugerth', 'Fukui', 'Mars', 'Liu', 'AGP',
                   'Nagel', 'Zeber')

# Plot
ggplot(df, aes(x = variable, y = value, fill=factor(dataset, levels = authors_order)))+
  facet_wrap(~factor(dataset, levels = authors_order), scales="fixed")+
  geom_violin()+
  geom_jitter(width = 0.1, size = 0.1)+
  theme_classic()+
  scale_y_log10()+
  scale_fill_brewer(palette="Set3")+
  # scale_fill_discrete(name = "Datasets")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.line = element_line(arrow = arrow(length = unit(0.1, "inches"))),
        legend.position="None")+
  # scale_x_discrete(breaks=variable_order,
  #                  labels=rep(c("raw", "QC"), times = 9))+
  labs(x = '', y = "Number of reads per sample", title = "Sequencing depth before/after processing raw 16S reads (dada2 pipeline)")

# Save figure
ggsave("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/04_QCplot/seqdepth.jpg", width=7, height=7)

