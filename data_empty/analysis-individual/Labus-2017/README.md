# Labus preprocessing

This folder will auto-populate when you run the [01_Dada2-Labus.Rmd](../../../scripts/analysis-individual/Labus-2017/01_Dada2-Labus.Rmd) and [03_EDA-Labus.Rmd](../../../scripts/analysis-individual/Labus-2017/03_EDA-Labus.Rmd) scripts.
- `01_Dada2-Labus/` will contain .rds files saved as checkpoints & also for a quicker html output of the 01_Dada2-Labus.Rmd notebook;
- `03_EDA-Labus/` will contain .rds files saved as checkpoints & also for a quicker html output of the 03_EDA-Labus.Rmd notebook, and it will also contain plots generated from that script;
- `filtered1/` will contain the .fastq files generated from primer removal;
- `filtered2/` will contain the .fastq files generated from quality filtering;
- `raw_fastq/` will contain the raw .fastq files downloaded from the ENA.