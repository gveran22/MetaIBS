# DADA2 - POZUELO DATASET

Pozuelo et al. (_Scientific Reports_, 2015) - [Reduction of butyrate- and methane-producing microorganisms in patients with Irritable Bowel Syndrome][1]

[1]: https://www.nature.com/articles/srep12693#Abs1


## Samples
- **290 total**
|                | IBS | Healthy |
| -------------- | ---:| ------: |
| Time point 1   | 125 |   66    |
| Time point 2   |  77 |   22    |
|   **Total**    | **202** | **88** |

## Data Quality
- **Technology** - Illumina MiSeq, single end
- **Nb of reads per sample** - mean of 59,556 reads per sample (23,370 - 96,888)
- **Read length** - ~300 bp
- **Quality** - excellent

## Primers
- V4 variable regions (about 250bp)
- FWD - 515F - 5’ - GTGCCAGCMGCCGCGGTAA - 3’
- REV -  806R - 5’ - GGACTACHVGGGTWTCTAAT - 3’
- reverse primer found at the end of most of the reads (in its revcomp orientation)

## Filtering
- **primers removal** - \~94% reads kept per sample.
- **quality filter** - \~100% reads kept per sample.

## Learn error rates
- parametric error model fits data

## Construct ASV table
### a) Infer sequence variants
- 57,879 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 23,849 seq variants (but still >77% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - ...
- Archaea - ...
- Eukaryota - ...
There was ... unassigned phyla. The final ASV table contains **... sequence variants**.

## Metadata
- 1st or 2nd collection time (1 month apart)
