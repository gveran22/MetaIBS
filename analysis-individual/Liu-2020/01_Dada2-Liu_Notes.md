# DADA2 - LIU DATASET

Liu et al. (_BMC Microbiology_, 2020) - [Microbial and metabolomic profiles in correlation with depression and anxiety co-morbidities in diarrhoea-predominant IBS patients][1]

[1]: https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-020-01841-4


## Samples
- **128 total**
- 44 Healthy & 84 IBS-Diarrhea

## Data Quality
- **Technology** - Illumina MiSeq, single end
- **Nb of reads per sample** - mean of 49,952 reads per sample (30,739 - 73,625)
- **Read length** - ~420 bp
- **Quality** - excellent

## Primers
- V3-V4 variable regions (about 450bp)
- FWD - 338F - 5’ - ACTCCTACGGGAGGCAGCA - 3’
- REV -  806R - 5’ - GGACTACHVGGGTWTCTAAT - 3’
- reverse complement of the reverse primer is found at the very end of >95% of the reads

## Filtering
- **primers removal** - not applied (requires a forward primer)
- **quality filter** - \~81% reads kept per sample.

## Learn error rates
- parametric error model fits data

## Construct ASV table
### a) Infer sequence variants
- 14,605 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 5,551 seq variants (but still >90% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - ...
- Archaea - ...
- Eukaryota - ...

All unassigned phyla were removed (n=...), samples below 500 total reads (n=...). The final ASV table contains **... sequence variants**.

## Metadata
- age
- BMI
- gender
- Bristol stool scale
- IBS symptom severity scale