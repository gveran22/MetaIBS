# DADA2 - ZHU DATASET

Zhu (_Frontiers in Cellular and Infection Microbiology_, 2019) - [Identification of Gut Microbiota and Metabolites Signature in Patients With Irritable Bowel Syndrome][1]

[1]: https://www.frontiersin.org/articles/10.3389/fcimb.2019.00346/full 


## Samples
- **29 total**
- 15 IBS & 14 HC

## Data Quality
- **Technology** - Illumina MiSeq, paired end
- **Nb of reads per sample** - mean of 44,384 reads per sample (21,956 - 55,164)
- **Read length** - ~220 bp
- **Quality** - excellent

## Primers
- V4 variable region (about 250bp)
- FWD - F515 - 5’ - GTGCCAGCMGCCGCGGTAA - 3’
- REV - R806 - 5’ - GGACTACHVGGGTWTCTAAT - 3’
- **Anomaly**: forward primer found in the _middle_ of >95% of forward reads (for all samples); but also found at the _end_ of reverse reads (only in half of the samples). Reverse primer not found (even with more mismatch allowed).

<p align="center">
<img src="./primer_anomaly.png" width="500" title="Primer-Anomaly-Schematic">
</p>


## Filtering
- **primer removal**: as about half the samples don't contain any primers in the reverse reads, primers were not removed. However, we will try to cut the FWD primer from the infered ASVs later on.
- **quality filtering**: \~80% of reads are kept

## Learn error rates
- parametric error model fits data

## Construct ASV table
### a) Infer sequence variants
- 5,307 amplicon sequence variants (ASVs)
The reads are much longer than the V4 region (~400 bp instead of expected 250bp). We didn't remove primers previously, but we know that the FWD primer was found in the middle of the FWD reads, so we cut off any sequence preceding the FWD primer in the inferred ASVs. The new (shortened) ASVs had duplicates, and we removed any ASV longer than 300bp, so we ended up with **1,631 ASVs** (all of them around 250bp).

### b) Remove chimeras
- 672 seq variants (but still >94% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - 672
- Archaea - 0
- Eukaryota - 0

All Eukaryota or unassigned phyla were removed (n=0). The final ASV table contains **672 sequence variants**.

## Metadata
- age
- gender



