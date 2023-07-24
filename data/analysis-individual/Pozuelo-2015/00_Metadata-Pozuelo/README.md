This directory contains metadata for the Pozuelo dataset.
- `Metadata-Pozuelo.csv`: final metadata table with essential information, cleaned up from `PozueloSraRunTable.txt` and merging with `host_subtype.csv`, using the [00_Metadata-Pozuelo.R](../../../../scripts/analysis-individual/Pozuelo-2015/00_Metadata-Pozuelo.R) script;
- `Pozuelo_hostsubtype.csv`: metadata table with IBS subtype, obtained from the [00_Metadata-Pozuelo.ipynb](../../../../scripts/analysis-individual/Pozuelo-2015/00_Metadata-Pozuelo.ipynb) notebook (need to get IBS subtype from BioSample "Description" slot on SRA website);
- `PozueloSraRunTable.csv`: metadata table downloaded from the [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA268708&o=acc_s%3Aa) (accession number PRJNA268708).
