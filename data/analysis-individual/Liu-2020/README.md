# Liu preprocessing

This folder will auto-populate when you run the [01_Dada2-Liu.Rmd](../../../../scripts/analysis-individual/Liu-2020/01_Dada2-Liu.Rmd) and [03_EDA-Liu.Rmd](../../../../scripts/analysis-individual/Liu-2020/03_EDA-Liu.Rmd) scripts.
- `01_Dada2-Liu/` will contain .rds files saved as checkpoints & also for a quicker html output of the 01_Dada2-Liu.Rmd notebook;
- `03_EDA-Liu/` will contain .rds files saved as checkpoints & also for a quicker html output of the 03_EDA-Liu.Rmd notebook, and it will also contain plots generated from that script;
- `filtered1(revcomp)/` will contain the .fastq files generated from primer removal (in reverse complement orientation);
- `filtered1/` will contain the .fastq files generated from primer removal (in original orientation);
- `filtered2/` will contain the .fastq files generated from quality filtering (in original orientation);
- `original(revcomp)` will contain the raw .fastq files (in reverse complement orientation).

For clarification: in the [01_Dada2-Liu.Rmd](../../../../scripts/analysis-individual/Liu-2020/01_Dada2-Liu.Rmd) script, this is how the raw .fastq files are processed:
1. Reverse complement the original .fastq files => saved in [original(revcomp)/](./original(revcomp)/)
2. Trim off reverse primer => saved in [filtered1(revcomp)/](./filtered1(revcomp)/)
3. Reverse complement the reads to get back to original orientation => saved in [filtered1/](./filtered1/)
4. Quality filter the reads => saved in [filtered2/](./filtered2/)