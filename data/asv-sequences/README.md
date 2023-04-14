# data/asv-sequences

This directory contains the ASVs sequences that were imputed from the DADA2 algorithm for each dataset. To note:
- chimeras were already removed (see the section "4.4. Remove chimeras" in "4.Construct ASV table" in any of the Dada2 R notebooks, e.g. [01_Dada2_Labus.Rmd](../../scripts/analysis-individual/Labus-2017/01_Dada2-Labus.Rmd));
- these ASV sequences are the ones that we align to the SILVA taxonomy (v.138) with the function `assignTaxonomy` (see the section "6.Taxonomic table" in any of the Dada2 R notebooks).
Here, I have saved these ASV sequences as FASTA files, so that taxonomy can be assigned with any other taxonomic reference (e.g. trying the AnnotIEM tool from other lab). You may import these FASTA files with the `readDNAStringSet("file-path-to-fasta-file")` function from the _Biostrings_ package.

To compare the output from using solely the SILVA reference database, to the output given by another reference database (such as AnnotIEM), I have saved the taxonomic tables of these ASVs in the [silva-tax](./silva-taxtable) directory, as _.rds_ files.

Very important: these ASV sequences or taxonomic tables **do not correspond to the final list of ASVs used in downstream statistical analyses**. After assigning taxonomy, a few more steps of data preprocessing were performed:
- I removed ASVs that belonged to Eukaryota and/or that had an unknown Phylum;
- I deleted samples with less than 500 total counts, and as some ASVs were present only in these low-count samples, they got removed as well with the “bad-quality” samples;
- for the `Pozuelo` dataset in particular, I removed 17 samples that were neither healthy nor IBS, and as some ASVs were present only in these non-healthy non-IBS samples, they got removed as well;
- for the `AGP` dataset in particular, after taxonomy assignment I had to remove bloom sequences as advised (see the section "7.3. Remove bloom sequences" in "7.Last preprocessing" in the [AGP R notebook](../../scripts/analysis-individual/AGP/01_Dada2-AGP.Rmd)).

Thus, there are more ASVs in this directory (both in the .fasta files and .rds taxtable files) than in the "final" [phyloseq objects](../phyloseq-objects/) shared on this github repository.  Sometimes only 1 ASV got removed (e.g. in the Nagel dataset), sometimes up to 1,928 ASVs got removed (e.g. in the AGP dataset).


## Additional note(s)
For the `Zhu` dataset, there was some sequencing bias (forward primer was found in the middle of the reads, see more explanation in the [01_Dada2-Zhu_Notes.md](../../scripts/analysis-individual/Zhu-2019/01_Dada2-Zhu_Notes.md)). After inferring the ASVs, I have attempted to "cut" the ASVs after the forward primer sequence, which showed promising result as the length of the "new" ASVs post-cut was 250bp (the expected length of the V4 region sequenced in this study). I have saved here these "new" ASVs in the `asv_zhu_cutFWDprimer.fasta` file. It would be interesting to compare the taxonomic assignment (maybe improved?) with the `asv_zhu.fasta` file that contains the "original" ASVs (with the sequencing bias).