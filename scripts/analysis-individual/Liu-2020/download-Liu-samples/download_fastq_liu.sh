#!/bin/bash
#Purpose: download Liu dataset
#Author: Salome Carcy March.2021
#--------------------------------------

INPUT=list_files_liu.txt

for FILE in $(cat list_files_liu.txt )
do
	echo "Run: $FILE"
	#~/sra-toolkit/fastq-dump --origfmt -X 2 -Z $FILE
	~/sra-toolkit/fastq-dump --origfmt --gzip $FILE
done < $INPUT