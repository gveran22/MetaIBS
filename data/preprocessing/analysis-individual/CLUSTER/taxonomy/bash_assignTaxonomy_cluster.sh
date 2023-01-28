#!/bin/bash
#$ -cwd
#$ -pe threads 4
#$ -l m_mem_free=16G
#$ -N assignTaxonomy_230128
#$ -o assignTaxonomy_230128.txt
#$ -e assignTaxonomy_230128.txt

module load EBModules
module load R/3.6.2-fosscuda-2019b
Rscript Rscript_assignTaxonomy_cluster.R