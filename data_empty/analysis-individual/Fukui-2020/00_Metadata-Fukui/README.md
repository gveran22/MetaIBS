This directory will contain metadata for the Fukui dataset.
1. Download the SraRunTable from the [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA637763&o=acc_s%3Aa) (accession number PRJNA637763);
2. Rename the SraRunTable to `FukuiSraRunTable.txt`;
3. Use the [00_Metadata-Fukui.R](../../../../scripts/analysis-individual/Fukui-2020/00_Metadata-Pozuelo.R) script to obtain the final `Metadata-Fukui.csv` table.