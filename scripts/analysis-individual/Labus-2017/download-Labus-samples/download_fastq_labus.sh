#!/bin/bash
#Purpose: download Labus dataset
#Author: Salome Carcy March.2021
#--------------------------------------

INPUT=list_files_labus.txt

for FILE in $(cat list_files.txt )
do
	echo "Run: $FILE"
	#~/sra-toolkit/fastq-dump --origfmt -X 2 -Z $FILE
	~/sra-toolkit/fastq-dump --split-files --origfmt --gzip $FILE
done < $INPUT
