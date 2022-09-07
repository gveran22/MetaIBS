# ASV-sequences

This directory contains the ASVs sequences that were imputed from the DADA2 algorithm. Chimeras were already removed (see the section "3. Remove chimeras" in "Construct sequence table" in any of the Dada2 R notebooks, i.e. [01_Dada2_Labus.Rmd](../analysis-individual/Labus-2017/01_Dada2-Labus.Rmd)).
In the simple pipeline, these ASV sequences were aligned to the SILVA taxonomy (v.138) with the function `assignTaxonomy`. Here, I have saved these ASV sequences as FASTA files, so that taxonomy can be assigned with the new AnnotIEM tool. You may import these FASTA files with the `readDNAStringSet("file-path-to-fasta-file")` function from the _Biostrings_ package.

To compare the output from using solely the SILVA reference database, to the output given by the AnnotIEM, I have saved the taxonomic tables of these ASVs in the [silva-tax](./silva-taxtable) directory, as _.rds_ files.

Be careful, these ASV sequences or taxonomic tables **do not correspond to the final list of ASVs used in downstream statistical analyses**. After assigning taxonomy, a few more steps of data preprocessing were performed:
- I removed ASVs that belonged to Eukaryota and/or that had an unknown Phylum;
- I deleted samples with less than 500 total counts, and as some ASVs were present only in these low-count samples, they got removed as well with the “bad-quality” samples;
- In the Pozuelo dataset in particular, I removed 17 samples that were neither healthy nor IBS, and as some ASVs were present only in these non-healthy non-IBS samples, they got removed as well.

Thus, there are more ASVs in this directory (both in the .fasta files and .rds taxtable files) than in the "final" [phyloseq objects](../phyloseq-objects/) shared on this github repository.  Sometimes only 1 ASV got removed (e.g. in the Nagel dataset), sometimes up to 1,928 ASVs got removed (e.g. in the AGP dataset).


## To note!!!
- for the **AGP** dataset, after taxonomy assignment, I would remove the bloom sequences as advised (see the section "4. Remove bloom sequences" in "Last preprocessing" in the [AGP R notebook](../analysis-individual/AGP/01_Dada2-AGP.Rmd)). Then I would check that the identified "bloom ASVs" belong mostly to Proteobacteria, as expected.
- for the **Zhu** dataset, there was some sequencing bias (forward primer was found in the middle of the reads, see more explanation in the [01_Dada2-Zhu_Notes.md](../analysis-individual/Zhu-2019/01_Dada2-Zhu_Notes.md), with the image). After inferring the ASVs, I have attempted to "cut" the ASVs after the forward primer sequence, which showed promising result as the length of the "new" ASVs post-cut was 250bp (the expected length of the V4 region sequenced in this study). I have saved here these "new" ASVs in the _asv_zhu_cutFWDprimer.fasta_ file. It would be interesting to compare the taxonomic assignment (maybe improved?) with the _asv_zhu.fasta_ file that contains the "original" ASVs (with the sequencing bias).