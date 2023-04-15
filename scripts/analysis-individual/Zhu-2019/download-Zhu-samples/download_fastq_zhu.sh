#!/bin/bash
#Purpose: download Zhu dataset
#Author: Salome Carcy March.2023
#--------------------------------------

INPUT=list_files_zhu.txt

for FILE in $(cat list_files_zhu.txt)
do
	echo "Run: $FILE"
	#~/sratoolkit/bin/fastq-dump --origfmt -X 2 -Z $FILE
	~/sratoolkit/bin/fastq-dump --split-files --origfmt --gzip $FILE
done < $INPUT
