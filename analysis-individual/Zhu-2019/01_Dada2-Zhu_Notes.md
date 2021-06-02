# DADA2 - ZHU DATASET

Zhu (_Frontiers in Cellular and Infection Microbiology_, 2019) - [Identification of Gut Microbiota and Metabolites Signature in Patients With Irritable Bowel Syndrome][1]

[1]: https://www.frontiersin.org/articles/10.3389/fcimb.2019.00346/full 


## Samples
- **29 total**
- 15 IBS & 14 HC

## Data Quality
- **Technology** - Illumina MiSeq, paired end
- **Nb of reads per sample** - mean of 44,384 reads per sample (21,956 - 55,164)
- **Read length** - ~250 bp
- **Quality** - excellent
- **Per base sequence content** - reverse reads seem to mostly have same sequence?...


## Primers
- V4 variable region (about 250bp)
- FWD - F515 - 5’ - GTGCCAGCMGCCGCGGTAA - 3’
- REV - R806 - 5’ - GGACTACHVGGGTWTCTAAT - 3’
- **Anomaly**: forward primer found in the _middle_ of >95% of forward reads (for all samples); but also found at the _end_ of reverse reads only in half of the samples. Reverse primer not found (even with more mismatch allowed).

<p align="center">
<img src="plots-zhu/primer_anomaly.png" width="500" title="Primer-Anomaly-Schematic">
</p>


## Filtering
- **primer removal**: fastq files deposited on the SRA/ENA databases don't contain Ilumina headers, which are needed to merge paired reads later on (if there is a different number of forward and reverse reads for each sample). Thus, primer sequences could not be removed.
- **quality filtering**: \~62% of reads are kept

## Learn error rates
- parametric error model fits data

## Construct ASV table
### a) Infer sequence variants
- 4,560 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 959 seq variants (but still >90% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - 944
- Archaea - 0
- Eukaryota - 15

All Eukaryota or unassigned phyla were removed (n=16). The final ASV table contains **943 sequence variants**.

## Metadata
- age
- gender



