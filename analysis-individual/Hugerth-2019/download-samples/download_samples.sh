#!/bin/bash
#Purpose: download Hugerth dataset
#Author: Salome Carcy Sept.2020
#--------------------------------------

INPUT=list-samples.txt

for FILE in $(cat list-samples.txt )
do
	echo "Run: $FILE"
	#~/sra-toolkit/fastq-dump --origfmt -X 2 -Z $FILE # print first 2 reads
	~/sra-toolkit/fastq-dump --split-files --origfmt --gzip $FILE
done < $INPUT
