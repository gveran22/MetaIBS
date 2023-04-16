# DADA2 - LABUS DATASET

Labus et al. (_Microbiome_, 2017) - [Differences in gut microbial composition correlate with regional brain volumes in irritable bowel syndrome][1]

[1]: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0260-z


## Samples
- **52 total**
- 23 Healthy & 29 IBS

## Data Quality
- **Technology** - 454 pyrosequencing, single end
- **Nb of reads per sample** - mean of 3,526 reads per sample (1,611 - 9,129)
- **Read length** - ~350 bp
- **Quality** - average (Q around 30)

## Primers
- V3-V5 variable regions (about 550bp)
- FWD - 341F - 5’ - CCTACGGGAGGCAGCAG - 3’
- REV -  926R - 5’ - CCGTCAATTCMTTTRAGT - 3’
- reverse primer found at the beginning of all of the reads (in its forward orientation)

## Filtering
- **primers removal** - 100% reads kept per sample.
- **quality filter** - \~73% reads kept per sample.

## Learn error rates
- parametric error model fits data

## Construct ASV table
### a) Infer sequence variants
- 726 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 710 seq variants (but still >98% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - 710
- Archaea - 0
- Eukaryota - 0
There were 3 unassigned phyla. The final ASV table contains **707 sequence variants**.

## Metadata
- age
- gender
- BMI
- IBS subtype