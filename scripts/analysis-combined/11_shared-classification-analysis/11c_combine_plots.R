# *****************************************************************************
# Purpose: Combine Figures for the paper
# Date: April 2023
# Author: (Minh) Viet Tran
# *****************************************************************************


# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(tidyverse)
library(patchwork)
library(gridExtra)
library(ggtree)

## 1.2. Data ####
path_root <- "~/MetaIBS" # CHANGE THIS ROOT DIRECTORY ON YOUR COMPUTER
path_intermed <- paste0(path_root,
                        "/data/analysis-combined/11_shared_classification/")

common_asvs <- readRDS(paste0(path_intermed, "common_asvs.RData")) +
  theme(plot.margin = unit(c(0,0,0,0) ,"cm"),
        panel.spacing = unit(0,"cm"),
        plot.tag = element_text(face = "bold"))

selected_asvs <- readRDS(paste0(path_intermed, "selected_asvs_table.RData"))
selected_asvs <- tableGrob(selected_asvs, rows = NULL, 
                           theme =
                             ttheme_minimal(base_size = 11,
                                            core =
                                              list(padding = unit(c(4.5, 1.5),
                                                                  "mm"))))
selected_asvs$heights[1] <- selected_asvs$heights[2]
highlight_tree <- readRDS(paste0(path_intermed, "highlight_tree.RData")) + 
  theme(plot.margin = unit(c(0, 0, 0, 0) ,"cm"),
        panel.spacing = unit(0,"cm"), legend.position = "none",
        plot.tag = element_text(face = "bold"))
sparse_log_tree <- readRDS(paste0(path_intermed, "sparse_tree.RData")) +
  theme(legend.position = "bottom", legend.spacing = unit(0.05, "cm"), 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        plot.tag = element_text(face = "bold")) +
  guides(colour = guide_legend(order = 1), 
         fill = guide_legend(order = 2, ncol = 1))

tree_sparse <- wrap_plots(sparse_log_tree, selected_asvs,
                          ncol = 1, heights = c(6, 2.5), guides = "keep")

layout <- "
AACCCC
AACCCC
AACCCC
BBCCCC
BBDDDD
BBDDDD
"

com_plot <- common_asvs + highlight_tree +  sparse_log_tree + selected_asvs +
  plot_layout(design = layout, tag_level = "new") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 12, face = "bold")) 


ggsave(filename = "sparse_log_contrast.jpg", plot = com_plot, device = "jpeg",
       path = path_intermed, width = 13, height = 11, dpi = 700)

ggsave(filename = "sparse_log_contrast.pdf", plot = com_plot, device = "pdf",
       path = path_intermed, width = 13, height = 11, dpi = 700)
