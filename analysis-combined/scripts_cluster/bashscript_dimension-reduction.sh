#!/bin/bash
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=8G
#$ -N 2022-01-10_dimension-reduction
#$ -o 2022-01-10_dimension-reduction.txt
#$ -e 2022-01-10_dimension-reduction.txt

module load EBModules
module load R/3.6.2-fosscuda-2019b
Rscript rscript_dimension-reduction.R