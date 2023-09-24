#!/bin/bash
#Purpose: download Fukui dataset
#Author: Salome Carcy March.2021
#--------------------------------------

INPUT=list_files_fukui.txt

for FILE in $(cat list_files_fukui.txt)
do
	echo "Run: $FILE"
	#~/sratoolkit/bin/fastq-dump --origfmt -X 2 -Z $FILE
	#~/sratoolkit/bin/fastq-dump --split-files --origfmt --gzip $FILE
	wget $FILE
done < $INPUT
