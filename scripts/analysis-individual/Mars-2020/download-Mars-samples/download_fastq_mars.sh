#!/bin/bash
#Purpose: download IBS & healthy samples from Mars et al. (2020)
#Author: Salome Carcy June.2021
#--------------------------------------

INPUT=list_files_mars.txt

for FILE in $(cat list_files_mars.txt )
do
	echo "Run: $FILE"
	#~/sra-toolkit/fastq-dump --origfmt -X 2 -Z $FILE
	~/sra-toolkit/fastq-dump --split-files --origfmt --gzip $FILE
done < $INPUT