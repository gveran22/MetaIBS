#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=1G
#$ -N fukui_phangorn
#$ -o fukui_phangorn.txt
#$ -e fukui_phangorn.txt

module load EBModules
module load R/3.6.2-fosscuda-2019b
Rscript 02_PhyloTree_Fukui.R