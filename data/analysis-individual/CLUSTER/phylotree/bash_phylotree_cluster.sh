#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=64G
#$ -N phylotree_230128
#$ -o phylotree_230128.txt
#$ -e phylotree_230128.txt

module load EBModules
module load R/3.6.2-fosscuda-2019b
Rscript Rscript_phylotree_cluster.R