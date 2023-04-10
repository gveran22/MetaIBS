# data/analysis-individual/

This repository contains `data` generated when you preprocess raw fastq files into ASV and taxonomic tables. Our meta-analysis included 13 datasets, so there are 13 sub-directories (+ `CLUSTER` sub-directory). Even though we used a standardized pipeline to preprocess the data for all datasets, we still have separate directories & scripts as there are subtle differences between datasets (e.g. cleaning up metadata tables, different primer sequences, ...). However, we structured all these sub-directories in a similar way.


# Structure of each dataset sub-directory

Within each sub-directory, you will find:
1. A `01_Dada2-NameDataset/` folder, in which you will save `.rds` files as checkpoints when running the `01_Dada2-NameDataset.Rmd` script (in the corresponding [scripts sub-directory](../../scripts/analysis-individual/));
2. A `filtered1/` and `filtered2/` folders, in which `.fastq` files will be generated during primer removal and quality-filtering of the reads (respectively), again using the `01_Dada2-NameDataset.Rmd` script (in the corresponding [scripts sub-directory](../../scripts/analysis-individual/));
3. A `03_EDA-NameDataset/` folder, in which you will save normalized phyloseq objects as `.rds` files, and also save plots when running the `03_EDA-NameDataset.Rmd` script (in the corresponding [scripts sub-directory](../../scripts/analysis-individual/)).


# CLUSTER sub-directory

This subdirectory is made to store data & scripts to run on your lab's computer clusters (server/high-performance computing). There are two tasks that we recommend running on a cluster:
1. **Assigning taxonomy**: after inferring ASVs, we align their sequence to the Silva database. You can copy the [CLUSTER/taxonomy/](./CLUSTER/taxonomy/) folder to your cluster & submit a job with the `bash_assignTaxonomy.sh` script, it will assign taxonomy for all ASV tables present in the [CLUSTER/taxonomy/input/](./CLUSTER/taxonomy/input/) subfolder. Taxonomic tables will be saved as an `.rds` object in the [CLUSTER/taxonomy/output/](./CLUSTER/taxonomy/output/) subfolder, so you can copy only the `output/` subfolder back to your computer.
2. **Phylogenetic trees**: to infer a phylogenetic tree and add it to your phyloseq object, you can copy the [CLUSTER/phylotree/](./CLUSTER/phylotree/) folder to your cluster & submit a job with the `bash_phylotree_cluster.sh` script. It will compute a phylogenetic tree for all phyloseq objects present in the [CLUSTER/phylotree/input/](./CLUSTER/phylotree/input/) subfolder. Phyloseq objects containing a phylogenetic tree will be saved as an `.rds` object in the [CLUSTER/phylotree/output/](./CLUSTER/phylotree/output/) subfolder, so you can copy only the `output/` subfolder back to your computer.
