#!/bin/bash

echo "Unzipping ena_files.zip"
unzip ena_files.zip
echo "Move fastq files to current directory"
find . -type f -name "*.fastq.gz" -exec mv {} . \;
echo "Remove empty directories"
find . -type d -empty -delete
echo "Remove ena_files.zip"
rm ena_files.zip