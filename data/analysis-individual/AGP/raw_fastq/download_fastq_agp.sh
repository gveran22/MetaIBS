#!/bin/bash
#Purpose: download IBS & healthy samples from AGP
#Author: Salome Carcy June.2021
#--------------------------------------

INPUT=list_files_agp.txt

for FILE in $(cat list_files_agp.txt)
do
	echo "Run: $FILE"
	#~/sratoolkit/bin/fastq-dump --origfmt -X 2 -Z $FILE
	#~/sratoolkit/bin/fastq-dump --origfmt --gzip $FILE
	wget $FILE
done < $INPUT