# original_revcomp (Liu)

This directory will contain .fastq files generated from reverse complementing the raw .fastq files with the [01_Dada2-Liu.Rmd](../../../../scripts/analysis-individual/Liu-2020/01_Dada2-Liu.Rmd) script.

The raw .fastq files only contain the reverse primer (no forward primer). We need to provide a `primer.fwd` argument to the DADA2 function to remove primers. We are thus reverse complementing the raw reads, to be able to use the reverse primers in the `primer.fwd` parameter. The reverse primers will be trimmed off (=> [filtered1(revcomp)](../filtered1(revcomp)/)), then we will reverse-complement again the reads to go back to the original direction (=> [filtered1](../filtered1/))