#!/bin/bash
#Purpose: download LoPresti dataset
#Author: Salome Carcy March.2021
#--------------------------------------

INPUT=list_files_lopresti.txt

for FILE in $(cat list_files_lopresti.txt )
do
	echo "Run: $FILE"
	#~/sratoolkit/bin/fastq-dump --origfmt -X 2 -Z $FILE
	# ~/sratoolkit/bin/fastq-dump --origfmt --gzip $FILE
	wget $FILE
done < $INPUT
