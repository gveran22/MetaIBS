#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=8G
#$ -N logratios
#$ -o logratios.txt
#$ -e logratios.txt

module load EBModules
module load R/3.6.2-fosscuda-2019b
Rscript 04_LogRatios-Taxa.R
