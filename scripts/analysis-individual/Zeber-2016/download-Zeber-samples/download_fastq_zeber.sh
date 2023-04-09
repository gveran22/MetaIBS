#!/bin/bash
#Purpose: download Zeber-Lubecka dataset
#Author: Salome Carcy March.2023
#--------------------------------------

INPUT=list_files_zeber.txt

for FILE in $(cat list_files_zeber.txt)
do
	echo "Run: $FILE"
	#~/sra-toolkit/fastq-dump --origfmt -X 2 -Z $FILE
	~/sra-toolkit/fastq-dump --split-files --origfmt --gzip $FILE
done < $INPUT