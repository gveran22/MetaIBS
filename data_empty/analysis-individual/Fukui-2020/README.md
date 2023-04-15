# Fukui preprocessing

This folder will auto-populate when you run the [01_Dada2-Fukui.Rmd](../../../scripts/analysis-individual/Fukui-2020/01_Dada2-Fukui.Rmd) and [03_EDA-Fukui.Rmd](../../../scripts/analysis-individual/Fukui-2020/03_EDA-Fukui.Rmd) scripts.
- `00_Metadata-Fukui/` will contain the SRA metadata dataframe and the final metadata table saved as .csv file from the 00_Metadata-Fukui script;
- `01_Dada2-Fukui/` will contain .rds files saved as checkpoints & also for a quicker html output of the 01_Dada2-Fukui.Rmd notebook;
- `03_EDA-Fukui/` will contain .rds files saved as checkpoints & also for a quicker html output of the 03_EDA-Fukui.Rmd notebook;
- `filtered_fastq/` will contain the .fastq files generated from quality filtering;