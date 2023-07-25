# DADA2 - POZUELO DATASET

Pozuelo et al. (_Scientific Reports_, 2015) - [Reduction of butyrate- and methane-producing microorganisms in patients with Irritable Bowel Syndrome][1]

[1]: https://www.nature.com/articles/srep12693#Abs1


## Samples
- **290 total**
|                | IBS | Healthy |
| -------------- | ---:| ------: |
| Time point 1   | 113 |   66    |
| Time point 2   |  72 |   22    |
|   **Total**    | **185** | **88** |

## Data Quality
- **Technology** - Illumina MiSeq, single end
- **Nb of reads per sample** - mean of 59,319 reads per sample (23,370 - 96,888)
- **Read length** - ~300 bp
- **Quality** - excellent

## Primers
- V4 variable regions (about 250bp)
- FWD - 515F - 5’ - GTGCCAGCMGCCGCGGTAA - 3’
- REV -  806R - 5’ - GGACTACHVGGGTWTCTAAT - 3’
- reverse primer found at the end of most of the reads (in its revcomp orientation)

## Filtering
- **primers removal** - \~94% reads kept per sample.
- **quality filter** - \~99.97% reads kept per sample.

## Learn error rates
- parametric error model fits data

## Construct ASV table
### a) Infer sequence variants
- 57,122 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 23,628 seq variants (but still >77% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - 23,581
- Archaea - 21
- Eukaryota - 0
There was 527 unassigned phyla. The final ASV table contains **23,101 sequence variants** in the phyloseq object.

## Metadata
- 1st or 2nd collection time (1 month apart)
