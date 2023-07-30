# DADA2 - ZHUANG DATASET

Zhuang et al. (_Frontiers in Microbiology_, 2018) - [Fecal microbiota alterations associated with diarrhea-predominant irritable bowel syndrome][1]

[1]: https://www.frontiersin.org/articles/10.3389/fmicb.2018.01600/full


## Samples
- **30 total**
- 20 IBS-D & 10 HC


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
- **primer removal**: 100% of forward reads kept, \~99% of reverse reads kept.
- **quality filtering**: \~69% of reads are kept


## Learn error rates
- parametric error model fits data


## Construct ASV table
### a) Infer sequence variants
- 6,099 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 948 seq variants (but still >88% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - 948
- Archaea - 0
- Eukaryota - 0

All Eukaryota or unassigned phyla were removed (n=0). The final ASV table contains **948 sequence variants**.

## Metadata
- IBS subtype (all IBS-D)