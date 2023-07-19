# DADA2 - MARS DATASET

Mars et al. (_Cell_, 2020) - [Longitudinal Multi-omics Reveals Subset-Specific Mechanisms Underlying Irritable Bowel Syndrome][1]

[1]: https://www.sciencedirect.com/science/article/pii/S0092867420309983


## Samples
- **72 total**
- 1st time point: 12 healthy & 27 IBS
- 2nd time point: 12 healthy & 21 IBS
=> Total 24 IBS & 48 healthy samples

Paper describes that there are **42 individuals** that underwent flexible sigmoidoscopy, some of them a second time 6 months later, for a total of **72 samples** (Fig1A). Details from supplementary table 1 (Table S1):
  - T1: 42 samples (14 HC - 28 IBS)
  - T2: 29 samples (9 HC - 20 IBS)
  => total of 71 samples (not 72...)

However, from the metadata table downloaded on the SRA, there are **72 total samples** (24 HC - 48 IBS). Three of the patient ID from the Table S1 could not be found in the SRA metadata table. More details on this in the [R markdown file](00_Metadata-Mars.Rmd) that explores this issue and builds a metadata table usable for downstream analyses.


## Data Quality
- **Technology** - Illumina MiSeq, paired end
- **Nb of reads per sample** - mean of 68,693 reads per sample (1,633 - 129,697)
- **Read length** - ~300 bp
- **Quality** - excellent

## Primers
- V4 variable region (about 250bp)
- FWD - F515 - 5’ - GTGCCAGCMGCCGCGGTAA - 3’
- REV - R806 - 5’ - GGACTACHVGGGTWTCTAAT - 3’
- forward primer found at the beginning of forward reads, and reverse complement of the reverse primer found at the end-ish of the forward reads (the opposite is true for reverse reads). Primers found in >90% of reads.

## Filtering
- **primer removal**: \~89% of forward reads kept, \~79% of reverse reads kept.
- **quality filtering**: \~70% of reads are kept.

## Learn error rates
- parametric error model fits data

## Construct ASV table
### a) Infer sequence variants
- 1,548 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 1,529 seq variants (but still >99% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - 1,526
- Archaea - 3
- Eukaryota - 0
There was no unassigned phyla, and 3 samples below 500 reads. The final ASV table contains **1,524 sequence variants**.

## Metadata
- IBS subtype
- gender
- age
- BMI
- collection time (0 - 6 months)