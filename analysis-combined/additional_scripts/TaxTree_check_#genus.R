##########################
# Purpose: Check number of genera for universal tax tree plot
# Date: March 2022
# Author: Salomé Carcy
##########################




##########
# IMPORT #
##########

# Libraries
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(trac)
library(treeio)
library(ggtree)
library(ggnewscale)
library(RColorBrewer)
library(scales)

# Data (agg tables)
path <- "~/Projects/IBS_Meta-analysis_16S/aggregated-tables"
aggTable.ringel   <- read.csv(file.path(path, "Ringel-2015/ringel_genus-agg.csv"))
aggTable.labus    <- read.csv(file.path(path, "Labus-2017/labus_genus-agg.csv"))
aggTable.lopresti <- read.csv(file.path(path, "LoPresti-2019/lopresti_genus-agg.csv"))
aggTable.pozuelo  <- read.csv(file.path(path, "Pozuelo-2015/pozuelo_genus-agg.csv"))
aggTable.zhuang   <- read.csv(file.path(path, "Zhuang-2018/zhuang_genus-agg.csv"))
aggTable.zhu      <- read.csv(file.path(path, "Zhu-2019/zhu_genus-agg.csv"))
aggTable.hugerth  <- read.csv(file.path(path, "Hugerth-2019/hugerth_genus-agg.csv"))
aggTable.fukui    <- read.csv(file.path(path, "Fukui-2020/fukui_genus-agg.csv"))
aggTable.mars     <- read.csv(file.path(path, "Mars-2020/mars_genus-agg.csv"))
aggTable.liu      <- read.csv(file.path(path, "Liu-2020/liu_genus-agg.csv"))
aggTable.agp      <- read.csv(file.path(path, "AGP-2021/agp_genus-agg.csv"))
aggTable.nagel    <- read.csv(file.path(path, "Nagel-2016/nagel_genus-agg.csv"))
aggTable.zeber    <- read.csv(file.path(path, "Zeber-2016/zeber_genus-agg.csv"))


# Aggregated phyloseq object to Genus level
physeq.glom <- readRDS("~/Projects/IBS_Meta-analysis_16S/data/analysis-combined/03_TaxonomicTree/physeq_all_glomGenus.rds")



###########
# COMPARE #
###########

# get a vector with all genera from aggregated tables
agg.genera <- unique(c(aggTable.labus$Genus,
                        aggTable.lopresti$Genus,
                        aggTable.pozuelo$Genus,
                        aggTable.zhuang$Genus,
                        aggTable.zhu$Genus,
                        aggTable.hugerth$Genus,
                        aggTable.fukui$Genus,
                        aggTable.mars$Genus,
                        aggTable.liu$Genus,
                        aggTable.agp$Genus,
                        aggTable.nagel$Genus,
                        aggTable.zeber$Genus))

# get vector with all genera from phyloseq object
phyloglom.genera <- unique(as.vector(tax_table(physeq.glom)[,"Genus"]))

# Compare length
length(agg.genera)
length(phyloglom.genera)

# Find genera that are different
setdiff(agg.genera, phyloglom.genera)
setdiff(phyloglom.genera, agg.genera)

# In the aggregated tables, there isn't the genus "EMP-G18"
# In the phyloseq object there isn't the genus "Dermacoccus", "Methylobacillus", "Caldanaerobacter", "Pelomonas" and "Amaricoccus"
