#!/bin/bash
#Purpose: download Pozuelo dataset
#Author: Salome Carcy March.2023
#--------------------------------------

INPUT=list_files_pozuelo.txt

for FILE in $(cat list_files_pozuelo.txt )
do
	echo "Run: $FILE"
	#~/sra-toolkit/fastq-dump --origfmt -X 2 -Z $FILE
	~/sra-toolkit/fastq-dump --split-files --origfmt --gzip $FILE
done < $INPUT
