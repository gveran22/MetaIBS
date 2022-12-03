# DADA2 - LO PRESTI DATASET

Lo Presti et al. (_Frontiers in Microbiology_, 2019) - [Fecal and mucosal microbiota profiling in irritable bowel syndrome and inflammatory bowel disease][1]

[1]: https://www.frontiersin.org/articles/10.3389/fmicb.2019.01655/full


## Samples
- **163 total**
|                | IBS | Healthy |
| -------------- | ---:| ------: |
|      Stool     |  36 |   40    |
| Sigmoïd mucosa |  44 |   46    |
|    **Total**   | **77** | **86** |

## Data Quality
- **Technology** - 454 pyrosequencing, single end
- **Nb of reads per sample** - mean of 1,468 reads per sample (15 - 8,792)
- **Read length** - ~400 bp
- **Quality** - average (Q around 25)

## Primers
- V1-V3 variable regions (about 450bp)
- FWD - 28F - 5’ - TTTGATCNTGGCTCAG - 3’
- REV -  519R - 5’ - GTNTTACNGCGGCKGCTG - 3’
- forward primer found at the beginning of most of the reads (in its forward orientation)

## Filtering
- **primers removal** - \~84% reads kept per sample.
- **quality filter** - \~54% reads kept per sample.

## Learn error rates
- parametric error model fits data

## Construct ASV table
### a) Infer sequence variants
- 1,094 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 1,085 seq variants (but still >99% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - 1,085
- Archaea - 0
- Eukaryota - 0
There was 0 unassigned phyla. There were 106 samples with total read count below 500. The final ASV table contains **1,085 sequence variants**.

## Metadata
- age
- gender
- IBS subtype
- sample type (sigmoid biopsy/stool)