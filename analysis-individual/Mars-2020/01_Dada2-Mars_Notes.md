# DADA2 - MARS DATASET

Mars et al. (_Cell_, 2020) - [Longitudinal Multi-omics Reveals Subset-Specific Mechanisms Underlying Irritable Bowel Syndrome][1]

[1]: https://www.sciencedirect.com/science/article/pii/S0092867420309983


## Samples
- **72 total**
- 1st time point: 12 healthy & 27 IBS
- 2nd time point: 12 healthy & 21 IBS
=> Total 24 IBS & 48 healthy samples

## Data Quality
- **Technology** - Illumina MiSeq, paired end
- **Nb of reads per sample** - mean of 68,693 reads per sample (1,633 - 129,697)
- **Read length** - ~300 bp
- **Quality** - excellent
- **Per base sequence content** - first few bases are always the same (corresponding to primer), then base content is fluctuating, as expected.

## Primers
- V4 variable region (about 250bp)
- FWD - F515 - 5’ - GTGCCAGCMGCCGCGGTAA - 3’
- REV - R806 - 5’ - GGACTACHVGGGTWTCTAAT - 3’
- forward primer found at the beginning of forward reads, and reverse complement of the reverse primer found at the end of the forward reads (the opposite is true for reverse reads). Primers found in >90% of reads.

## Filtering
- **primer removal**: fastq files deposited on the SRA/ENA databases don't contain Ilumina headers, which are needed to merge paired reads later on. Thus, primer sequences could not be removed (the first 20bp and last 30bp were trimmed instead).
- **quality filtering**: \~64% of reads are kept.