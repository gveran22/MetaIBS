# Nagel preprocessing

This folder will auto-populate when you run the [01_Dada2-Nagel.Rmd](../../../../scripts/analysis-individual/Nagel-2016/01_Dada2-Nagel.Rmd) and [03_EDA-Nagel.Rmd](../../../../scripts/analysis-individual/Nagel-2016/03_EDA-Nagel.Rmd) scripts.
- `01_Dada2-Nagel/` will contain .rds files saved as checkpoints & also for a quicker html output of the 01_Dada2-Nagel.Rmd notebook;
- `03_EDA-Nagel/` will contain .rds files saved as checkpoints & also for a quicker html output of the 03_EDA-Nagel.Rmd notebook, and it will also contain plots generated from that script;
- `filtered1/` will contain the .fastq files generated from removing primers;
- `filtered2/` will contain the .fastq files generated from quality filtering;
- `original_trim/` will contain the .fastq files generated from removing reads with less than 50bp (raw fastq files, with only reads >50bp).