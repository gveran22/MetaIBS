# DADA2 - FUKUI DATASET

Fukui (_J. Clin. Med_) - [Usefulness of Machine Learning-Based Gut Microbiome Analysis for Identifying Patients with Irritable Bowels Syndrome][1]

[1]: https://www.mdpi.com/2077-0383/9/8/2403


## Samples
- **111 total**
- 85 IBS & 26 HC

## Data Quality
- **Technology** - Illumina MiSeq, paired end
- **Nb of reads per sample** - mean of 41,870 reads per sample (22,487 - 72,649)
- **Read length** - ~250 bp
- **Quality** - excellent

## Primers
- V1-V2 variable region (about 300bp)
- FWD - F27 - 5' - TTTGATCCTGGCTCAG - 3'
- REV - R338 - 5' - TGCCTCCCGTAGGAGT - 3'
- forward primer found at the beginning of forward reads (same for reverse primer in reverse reads). Primers found in >95% of reads.

## Filtering
- **primer removal**: \~98% of reads are kept
- **quality filtering**: \~85% of reads are kept

## Learn error rates
- parametric error model fits data

## Construct ASV table
### a) Infer sequence variants
- 10,363 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 6,801 seq variants (but still >96% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - 6,801
- Archaea - 0
- Eukaryota - 0
There were 4 unassigned phyla removed There was 1 sample below 500 total reads.
The final ASV table contains **6,797 sequence variants**.

## Metadata
- none (IBS vs HC)



