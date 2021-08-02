# DADA2 - RINGEL-KULKA DATASET

Ringel-Kulka et al. (_Am J Physiol Gastrointest Liver Physiol_, 2016) - [Molecular characterization of the intestinal microbiota in patients with and without abdominal bloating][1]

[1]: https://journals.physiology.org/doi/full/10.1152/ajpgi.00044.2015


## Samples
- **76 total**
- 20 Healthy & 56 IBS

## Data Quality
- **Technology** - 454 pyrosequencing, single end
- **Nb of reads per sample** - mean of 13,530 reads per sample (118 - 23,762)
- **Read length** - ~350 bp
- **Quality** - average (Q around 30)

## Primers
- V1-V2 variable regions (about 350bp)
- FWD - 8F - 5’ - TTTGATCMTGGCTCAG - 3’
- REV -  357R - 5’ - CTGCTGCCTYCCGTA - 3’
- forward primer found at the beginning of all of the reads & reverse primer found at the end of 50-60% of the reads

## Filtering
- **primers removal** - \~58% reads kept per sample.
- **quality filter** - \~98% reads kept per sample.

## Learn error rates
- parametric error model fits data

## Construct ASV table
### a) Infer sequence variants
- 3,438 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 3,264 seq variants (but still >98% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - 3,264
- Archaea - 0
- Eukaryota - 0
There was 2 unassigned phyla. The final ASV table contains **3,262 sequence variants**.

## Metadata
- none (no IBS/HC label)