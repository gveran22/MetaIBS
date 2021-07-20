# DADA2 - NAGEL DATASET

Nagel et al. (_Microbiome_, 2016) - [Comparison of faecal microbiota in Blastocystis-positive and Blastocystisnegative irritable bowel syndrome patients][1]

[1]: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0191-0


## Samples
- **30 total**
- 15 Healthy & 15 IBS

## Data Quality
- **Technology** - Ion Torrent, single end
- **Nb of reads per sample** - mean of 40,420 reads per sample (18,838 - 52,361)
- **Read length** - ~300 bp
- **Quality** - average (Q around 30)

## Primers
- V4 variable regions (about 450bp)
- FWD - 515F - 5’ - GTGCCAGCMGCCGCGGTAA - 3’
- REV -  806R - 5’ - GGACTACHVGGGTWTCTAAT - 3’
- forward primer found at the beginning of the reads & reverse primer found at the end of \~50% of the reads

## Filtering
- **primers removal** - \~53% reads kept per sample.
- **quality filter** - \~83% reads kept per sample.

## Learn error rates
- parametric error model was manually modified to fit better the data

## Construct ASV table
### a) Infer sequence variants
- 1,322 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 1,092 seq variants (>95% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - ...
- Archaea - ...
- Eukaryota - ...

All unassigned phyla were removed (n=...), samples below 500 total reads (n=...). The final ASV table contains **... sequence variants**.

## Metadata
- age
- gender
- IBS subtype (all IBS-D)