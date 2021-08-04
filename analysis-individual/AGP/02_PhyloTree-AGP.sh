#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=16G
#$ -N AGP_phangorn
#$ -o AGP_phangorn.txt
#$ -e AGP_phangorn.txt

module load EBModules
module load R/3.6.2-fosscuda-2019b
Rscript 02_PhyloTree_AGP.R