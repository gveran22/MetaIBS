#!/bin/bash
#Purpose: demultiplex Nagel data
#Author: Salome Carcy September.2023
#--------------------------------------

# Install sabre
echo "installing sabre"
unzip master.zip
cd sabre-master
make
cd ..

# Demultiplex with sabre
echo "demultiplex"
./sabre-master/sabre se -f ./multiplexed_data/TOR-143-Pool_1.fastq -b sabre_formatted_barcode_Pool1.txt -u unknown_barcode_Pool1.fastq
./sabre-master/sabre se -f ./multiplexed_data/ION-184-Pool_2.fastq -b sabre_formatted_barcode_Pool2.txt -u unknown_barcode_Pool2.fastq

# Move all fastq files to "original" directory
echo "move fastq files"
mv *.fastq ../raw_fastq/
cd ../raw_fastq/

# Keep the unknown barcodes in the "original_demultiplex" directory
mv unknown* ../raw_fastq_multiplexed/
# Remove the blastocystis positive samples
rm B*.fastq
rm HB*.fastq
rm IBD*.fastq