This directory contains metadata for the Pozuelo dataset.
- `host_subtype.csv`: metadata table with IBS subtype, obtained from the [00_Metadata-Pozuelo.ipynb](../00_Metadata-Pozuelo.ipynb) notebook (need to get IBS subtype from BioSample "Description" slot on SRA website);
- `Metadata-Pozuelo.csv`: final metadata table with essential information, cleaned up from `PozueloSraRunTable.txt` and merging with `host_subtype.csv`, using the [00_Metadata-Pozuelo.R](../00_Metadata-Pozuelo.R) script;
- `PozueloSraRunTable.csv`: metadata table downloaded from the [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA268708&o=acc_s%3Aa) (accession number PRJNA268708).
