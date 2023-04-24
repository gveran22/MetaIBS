# *****************************************************************************
# Purpose: Run Sparse Log Contrast on shared ASVS
# Date: April 2023
# Author: (Minh) Viet Tran
# *****************************************************************************

# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(tidyverse)
library(phyloseq)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(reticulate)
library(trac)

## 1.2  Define parameters ####
# define seed for reproducibility
seed <- 123
# define threshold how often ASVs have to appear
threshold <- 0.1

## 1.3. Data ####
path_root <- "~/MetaIBS" # CHANGE THIS ROOT DIRECTORY ON YOUR COMPUTER
path_intermed <- paste0(path_root,
                        "/data/analysis-combined/11_shared_classification/")
path_phylobj <- 
  file.path(path_root,
            "/data/analysis-combined/05_Common-ASVs/physeq_commonASV_Nagel-Pozuelo.rds")
data <- readRDS(path_phylobj)
dim(data@otu_table)

## 1.4. Preprocessing ####

data <- subset_taxa(data , Kingdom == "Bacteria")
data <- subset_samples(data, Collection == "1st")
# data <- subset_samples(data, host_subtype %in% c("HC", "IBS-D"))

# ASV have appear at least in 10% of the observations
freq <- colSums(sign(data@otu_table@.Data))
data <- phyloseq::prune_taxa(freq >= threshold * phyloseq::nsamples(data),
                             data)
dim(data@otu_table)

# Preprocess the data for the sparse log-contrast analysis
preprocess_data <- function(data_phylo) {
  # get tax table from phyloseq file
  tax <- phyloseq::tax_table(data_phylo)
  tax <- as.data.frame(tax)
  # replace NA with blank space --> for further preprocessing
  tax[is.na(tax)] <- ""
  # add prefix to each column
  prefix <- paste0(c("p", "c", "o", "f", "g", "s"), "__")
  for (i in 1:6) tax[, i + 1] <- paste0(prefix[i], tax[, i + 1])
  # add a root node
  # tax <- cbind("Life", tax); colnames(tax)[1] <- "Rank0"
  # add ASVs
  tax <- cbind(tax, rownames(tax))
  colnames(tax)[ncol(tax)] <- "ASV"
  
  # make it so labels are unique
  for (i in seq(2, 7)) {
    # add a number when the type is unknown... e.g. "g__"
    ii <- nchar(tax[, i]) == 3
    if (sum(ii) > 0)
      tax[ii, i] <- paste0(tax[ii, i], 1:sum(ii))
  }
  
  for (i in 2:7) {
    tax[, i] <- paste(tax[, i - 1], tax[, i], sep = "::")
  }
  
  tax <- as.data.frame(tax, stringsAsFactors = TRUE)
  for (i in seq_len(ncol(tax))) {
    tax[, i] <- as.factor(tax[, i])
  }
  
  # form phylo object:
  tree1 <- trac::tax_table_to_phylo(~ Kingdom / Phylum / Class / Order /
                                      Family / Genus / Species / ASV,
                                    data = tax, collapse = TRUE)
  
  # convert this to an A matrix to be used for aggregation:
  a_matrix <- trac::phylo_to_A(tree1)
  
  list(y = sample_data(data_phylo)$host_disease,
       x = otu_table(data_phylo)@.Data,
       sample_data = as_tibble(asS3(sample_data(data_phylo))),
       tax =  as.data.frame(tax_table(data_phylo)@.Data),
       tree = tree1,
       A = a_matrix)
}

dat <- preprocess_data(data_phylo = data)

set.seed(seed)
log_pseudo <- function(x, pseudo_count = 1) log(x + pseudo_count)
z_matrix <- log_pseudo(dat$x)
dim(z_matrix)
y <- dat$y
table(y)
y <- (y == "IBS") * 2 - 1

## 1.5. Run analysis ####


fit_classif <- sparse_log_contrast(Z = z_matrix, y = y, 
                                   method = "classif", min_frac = 1e-03)

cvfit_classif <- cv_sparse_log_contrast(fit_classif, 
                                        Z = z_matrix, y = y, 
                                        nfolds = 5)
authors <- c("Nagel", "Pozuelo")
da_models <- c("LinDA", "ANCOMBC", "sccoda")
suffix <- c("a1", "a2", "combined")

# get non-zero coefficients
show_nonzeros <- function(x) x[x != 0]
nzz <- show_nonzeros(fit_classif$beta[, cvfit_classif$cv$i1se])

# Define named list with all the selected coefficents
nzz_list <- list(combined = nzz)

# plot cv
plot_cv_log_contrast <- function(cvfit) {
  # extract number of non-zero coefficients, mean error and sd
  tibble(non_zero = cvfit$cv$nonzeros,
         error = cvfit$cv$m,
         sd = cvfit$cv$se) %>%
    ggplot(aes(x = non_zero, y = error)) +
    theme_classic() +
    # 1se line
    geom_vline(xintercept = cvfit$cv$nonzeros[cvfit$cv$i1se],
               colour = "gray", linetype = "longdash") +
    # ibest line
    geom_vline(xintercept = cvfit$cv$nonzeros[cvfit$cv$ibest],
               colour = "gray50", linetype = "longdash") +
    # line + point and error bars
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = error - sd, ymax = error + sd), width = .2,
                  position = position_dodge(0.05))
}
plot_cv_log_contrast(cvfit =  cvfit_classif)


# plot heatmap

# get heatmap
get_overview_plot <- function(tax, selected_ASVs) {
  # takes the taxonomic table and a named list of selected ASVs and create a
  # plot
  # 1. transform list to dataframe
  # 2. filter the taxonomic tbl for the selected ASVs
  # 3. combine the name and create heatmap
  
  selected_ASVs <- selected_ASVs %>%
    map_dfr(enframe, .id = "model")
  rownames_to_column(tax) %>%
    filter(rowname %in% selected_ASVs$name) %>%
    rename(name = rowname) %>%
    right_join(selected_ASVs, by = "name") %>%
    select(-name,) %>%
    mutate(across(where(is.character), ~replace_na(.x, "_"))) %>%
    mutate(name = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species,
                        sep = "*")) %>%
    select(name, model, value) %>%
    mutate(value = (value != 0)) %>%
    complete(model, name) %>%
    ggplot(aes(x = model, y = name, fill = value)) +
    geom_tile() +
    scale_fill_brewer(name = "Found?", na.value = "grey") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.text = element_text(size = 6)) +
    labs(x = "Author", y = "ASV")
}

get_da_results <- function(path, authors, da_models, suffix) {
  da_results <- list()
  for (m in da_models) {
    for (s in suffix) {
      da_results[[m]][[s]] <- read.csv(
        paste0(path, "shared_", authors[1], authors[2], "_", m, "_None_",
               s, ".csv")) %>% filter(is_da == "True")
    }
    if (m == "ANCOMBC") {
      da_results[[m]][[s]] <- read.csv(
        paste0(path, "shared_", authors[1], authors[2], "_", m, "_None_",
               s, ".csv")) %>% filter(is_da == "True")
    }
    if (m == "LinDA") {
      da_results[[m]][["combined"]] <- read.csv(
        paste0(path, "shared_", authors[1], authors[2], "_", m, "_None_",
               "combined.csv")) %>% filter(is_da == "True")
    }
    da_results[[m]] <- da_results[[m]] %>%
      map_dfr(.f = I, .id = "author")
  }
  
  da_results %>%
    map_dfr(.f = I, .id = "model") %>%
    mutate(effect_size = case_when(model == "LinDA" ~ stat,
                                   model == "ANCOMBC" ~ beta,
                                   model == "sccoda" ~ Final.Parameter)) %>%
    select(model, author, OTU, Type, effect_size) %>%
    mutate(author = str_replace_all(author, "a1", authors[1])) %>%
    mutate(author = str_replace_all(author, "a2", authors[2])) %>%
    mutate(author = tolower(author))
}


get_overview_plot(dat$tax, nzz_list)

# Compare with DA testing

compare_da <- function(path, authors, da_models, suffix, tax, selected_ASVs) {
  # read in results from DA testing
  da_results <- get_da_results(path = path,
                               authors = authors,
                               da_models = da_models,
                               suffix = suffix)
  
  # process results fro sparse log-contrast model
  selected_ASVs <- selected_ASVs %>%
    map_dfr(enframe, .id = "author")
  
  results_slc <- rownames_to_column(tax) %>%
    filter(rowname %in% selected_ASVs$name) %>%
    rename(name = rowname) %>%
    right_join(selected_ASVs, by = "name") %>%
    mutate(across(where(is.character), ~replace_na(.x, "_"))) %>%
    mutate(Type = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species,
                        sep = "*")) %>%
    rename(OTU = name) %>%
    select(Type, OTU, author) %>%
    mutate(model = "sparse log-contrast") %>%
    bind_rows(da_results) %>%
    mutate(found = "yes")
  results_slc_unique <- distinct(results_slc, OTU, .keep_all = TRUE) %>%
    arrange(desc(Type))
  
  results_slc %>%
    complete(OTU, author, model) %>%
    ggplot(aes(x = model, y = OTU, fill = found)) +
    geom_tile() +
    scale_y_discrete(limit = results_slc_unique$OTU,
                     labels = results_slc_unique$Type) +
    scale_fill_brewer(name = "Found?", na.value = "grey") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.text = element_text(size = 6, angle = 45, hjust = 1)) +
    labs(x = "Method", y = "ASV") +
    facet_grid(~author)
}


# define path to DA analysis results
path_da <- file.path(path_root,
                          "data/analysis-combined/10_DA-analysis/shared_asvs/")
compare_da(path = path_da, authors = authors, da_models = da_models,
           suffix = suffix, tax = dat$tax, selected_ASVs = nzz_list)


generate_tree_sparse_log <- function(tax, tree, path, authors, da_models, suffix,
                                     nzz, author_selected = NULL, path_output) {
  # plot taxonomic tree with points marked by trac and outer ring
  # indicating the results of the other methods
  #
  # tax: output of preprocessing step
  # tree: output of preprocessing step
  # path: path of the DA files
  # authors: name of authors
  # da_models: name of da methods
  # suffix: part of the file name
  # nzz: number of non-zero sparse log-contrast
  # nzz_trac: number of non-zero trac
  # author_selected: for single dataset
  # Define coloring of tree
  phylum_group <- rownames_to_column(tax) %>%
    select(Phylum, rowname) %>%
    mutate(Phylum = case_when(
      Phylum == "Actinobacteriota" ~ "Actinobacteriota",
      Phylum == "Bacteroidota" ~ "Bacteroidota",
      Phylum == "Proteobacteria" ~ "Proteobacteria",
      Phylum == "Firmicutes" ~ "Firmicutes",
      Phylum == "Verrucomicrobiota" ~ "Verrucomicrobiota",
      TRUE ~ "Other"
    ))
  phylum_group <- split(phylum_group$rowname, phylum_group$Phylum)
  tree_1 <- groupOTU(tree, phylum_group)
  
  # generate tree on circle
  
  cols = c("Actinobacteriota" = "#1B9E7799",
           "Bacteroidota" = "#D95F0299",
           "Firmicutes" = "#7570B399",
           "Proteobacteria" = "#E7298A99",
           "Verrucomicrobiota" = "#66A61E99",
           "Other" = "#E6AB0299")
  
  tree_2 <- ggtree(tree_1, aes(color = group),
                   layout = "circular", size = 0.5) +
    scale_color_manual(values = cols, name = "Phylum")
  tree_2 <- open_tree(tree_2, angle = 22)
  tree_2 <- tree_2 %>%
    rotate_tree(-270)
  
  # get DA results
  da_results <- get_da_results(path = path,
                               authors = authors,
                               da_models = da_models,
                               suffix = suffix) %>%
    rename(ID = OTU) %>%
    {if (!is.null(author_selected)) filter(., author == author_selected)
      else filter(., author == "combined")} %>%
    rename(value = effect_size) %>%
    select(ID, model, value) %>%
    mutate(found = "Yes")
  
  # get sparse log-contrasts results and combine with DA results
  sparse_results <- enframe(nzz) %>%
    rename(ID = name) %>%
    mutate(found = "Yes") %>%
    mutate(model = "sparse") %>%
    bind_rows(da_results)
  
  # Add other not selected ASVs --> for outer ring otherwise we can
  # not assign a color
  
  other_asvs <- !(tree$tip.label %in% names(nzz))
  other_asvs_df <- data.frame(ID = tree$tip.label[other_asvs],
                              model = "sparse",
                              found = NA)
  sparse_results <- bind_rows(sparse_results, other_asvs_df)
  sparse_results <- sparse_results %>%
    complete(ID, model) %>%
    mutate(found = replace_na(found, "No")) %>%
    mutate(found = factor(found, levels = c("No", "Yes"))) %>%
    mutate(model = str_replace_all(model, "sccoda", "scCODA")) %>%
    mutate(model = str_replace_all(model, "LinDa", "LinDA")) %>%
    mutate(model = str_replace_all(model, "ANCOMBC", "ANCOM-BC")) %>%
    mutate(model = str_replace_all(model, "sparse", "Sparse")) %>%
    mutate(model = factor(model)) %>%
    mutate(model = fct_relevel(model, "Sparse"))
  # add the points selected by trac and outer ring
  df_subset_n_larger_four <- sparse_results %>%
    group_by(ID) %>%
    filter(found == "Yes") %>%
    summarise(n = n()) %>%
    filter(n == 4)
  
  
  
  tax_subset_n_larger_four <- tax %>%
    rownames_to_column() %>%
    filter(rowname %in% df_subset_n_larger_four$ID) %>%
    mutate(across(.cols = everything(), ~replace_na(., "*"))) %>%
    mutate(name = paste(Phylum, Class, Order, Family,
                        Genus, Species, sep = "_")) %>%
    rename(ID = rowname) %>%
    select(ID, name) %>%
    rename(label2 = name) %>%
    rename(label = ID) %>%
    arrange(desc(label)) %>%
    rownames_to_column() %>%
    rename(label3 = rowname) %>%
    relocate(label)
  
  sparse_results_tmp <- sparse_results %>%
    rename(label = ID) %>%
    filter(!is.na(value)) %>%
    mutate(value = round(value, digits = 2)) %>%
    pivot_wider(id_cols = label, names_from = model, values_from = value)
  
  tax_levels <- c("Phylum", "Class",
                  "Order", "Family", "Genus", "Species")
  
  table_selected <-  tax_subset_n_larger_four %>%
    left_join(sparse_results_tmp, by = "label") %>%
    #      select(-model, -found) %>%
    rename(ID = label3) %>%
    separate(label2,
             into = as.character(tax_levels),
             sep = "_",
             remove = FALSE,
             fill = "right") %>%
    select(-label, -label2) %>%
    mutate(Effect = sign(Sparse))
  all(table_selected$Effect == sign(table_selected$`ANCOM-BC`))
  all(table_selected$Effect == sign(table_selected$scCODA))
  all(table_selected$Effect == sign(table_selected$LinDA))
  table_selected <- table_selected %>%
    mutate(Effect = case_when(Effect == (-1) ~ "-",
                              Effect == 1 ~ "+")) %>%
    select(-"ANCOM-BC", -"LinDA", -"scCODA", -"Sparse")
  
  tree_2 <- tree_2 %<+% tax_subset_n_larger_four
  tree_3 <- tree_2 +
    #    geom_tiplab(aes(label = label2), size = 2, offset = 2.4, color = "black") +
    geom_label2(aes(label = label3), nudge_x = 3, label.size = 0,
                show.legend = FALSE, fill = "grey50", color = "white",
                label.r = unit(0.6, "lines"), label.padding = unit(0.4, "lines")) +
    new_scale_fill() +
    geom_fruit(data = sparse_results, geom = geom_tile,
               mapping = aes(y = ID, x = model, fill = found),
               axis.params = list(axis = "x", text.angle = 0, text.size = 3,
                                  line.size = 0, hjust = 0),
               alpha = 0.9,
               offset = 0.04, size = 0.02, color = "grey50") +
    scale_fill_manual(values = c("ghostwhite", "grey50"), name = "Signal") +
    theme(plot.margin = unit(c(0,0,0,0) ,"cm"),
          panel.spacing = unit(0,"cm"))
  
  saveRDS(table_selected,
          file = paste0(path_output, "selected_asvs_table.RData"))
  
  tree_3
}



sparse_tree <- generate_tree_sparse_log(tax = dat$tax, tree = dat$tree, 
                                        path = path_da, 
                                        authors = authors,
                                        da_models = da_models, suffix = suffix, 
                                        nzz = nzz, path_output = path_intermed)

sparse_tree
saveRDS(sparse_tree, file = paste0(path_intermed, "sparse_tree.RData"))
# open_tree(sparse_tree, angle = 20)
ggsave(filename = paste0(path_intermed, "sparse_tree.pdf"), 
       plot = sparse_tree, device = 'pdf', dpi = 300)

