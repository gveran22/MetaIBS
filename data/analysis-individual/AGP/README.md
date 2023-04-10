# AGP preprocessing

This folder will auto-populate when you run the [01_Dada2-AGP.Rmd](../../../scripts/analysis-individual/AGP/01_Dada2-AGP.Rmd) and [03_EDA-AGP.Rmd](../../../scripts/analysis-individual/AGP/03_EDA-AGP.Rmd) scripts.
- `01_Dada2-AGP/` will contain .rds files saved as checkpoints & also for a quicker html output of the 01_Dada2-AGP.Rmd notebook;
- `03_EDA-AGP/` will contain .rds files saved as checkpoints & also for a quicker html output of the 03_EDA-AGP.Rmd notebook;
- `filtered_fastq/` will contain the .fastq files generated from quality filtering;
- `newbloom.all.fna` contains a list of OTUs known to grow at room temperature, that [bias the microbiome composition of samples from the AGP dataset](https://journals.asm.org/doi/10.1128/mSystems.00199-16#B5). This file was downloaded from the [github repository](https://github.com/knightlab-analyses/bloom-analyses) of the publication.