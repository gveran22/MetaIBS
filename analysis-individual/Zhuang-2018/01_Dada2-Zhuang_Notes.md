# DADA2 - ZHUANG DATASET

Zhuang et al. (_Frontiers in Microbiology_, 2018) - [Fecal microbiota alterations associated with diarrhea-predominant irritable bowel syndrome][1]

[1]: https://www.frontiersin.org/articles/10.3389/fmicb.2018.01600/full


## Samples
- **30 total**
- 15 IBS & 15 HC


## Data Quality
- **Technology** - Illumina MiSeq, paired end
- **Nb of reads per sample** - mean of 64,849 reads per sample (6,793 - 84,048)
- **Read length** - ~300 bp
- **Quality** - excellent


## Primers
- V3-V4 variable region (about 430bp)
- FWD - F338 - 5’ - ACTCCTACGGGAGGCAGCA - 3’
- REV - R806 - 5’ - GGACTACHVGGGTWTCTAAT - 3’
- forward primer found at the beginning of forward reads (same for reverse primer & reverse reads)


## Filtering
- **primer removal**: fastq files deposited on the SRA/ENA databases don't contain Ilumina headers, which are needed to merge paired reads later on (if there is a different number of forward and reverse reads for each sample). Thus, reads not containing primer sequences could not be filtered out. However, >95% of the reads contain both primers, so we simply trimmed the first 25bp of the reads to remove the primers.
- **quality filtering**: \~70% of reads are kept


## Learn error rates
- parametric error model fits data


## Construct ASV table
### a) Infer sequence variants
- 5,613 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 942 seq variants (but still >88% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - ...
- Archaea - ...
- Eukaryota - ...

All Eukaryota or unassigned phyla were removed (n=...). The final ASV table contains **... sequence variants**.

## Metadata
- IBS subtype (all IBS-D)