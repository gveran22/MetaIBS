This directory contains metadata for the LoPresti dataset.
- `LoPrestiSraRunTable.csv`: metadata table downloaded from the [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA391149&o=acc_s%3Aa) (accession number PRJNA391149);
- `Metadata-LoPresti-final.csv`: final metadata table (without samples deleted through DADA2 pipeline).
- `Metadata-LoPresti.csv`: metadata table with essential information, cleaned up from `LoPrestiSraRunTable.csv` with the [00_Metadata-LoPresti.Rmd](../../../../scripts/analysis-individual/LoPresti-2019/00_Metadata-LoPresti.Rmd) script.