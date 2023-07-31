# Dataset
Hugerth (_Gut Microbiota_, 2019) - [No distinct microbiome signature of irritable bowel syndrome found in a Swedish random population][1]

[1]: https://gut.bmj.com/content/69/6/1076


## Samples
Paper describes that there is:
  - **561 samples**
  - 376 sigmoid mucosa samples (63 IBS - 313 HC)
  - 185 fecal samples (32 IBS - 153 HC)

However, from the metadata table downloaded on the SRA, there are **607 total samples** (140 IBS - 467 HC). More details on this in the [R markdown file](00_Metadata-Hugerth.Rmd) that explores this issue and builds a metadata table usable for downstream analyses.


## Data Quality
- **Technology** - Illumina MiSeq
- **Nb of reads per sample** - mean of 31,163 reads per sample (9 - 306,468)
- **Read length** - ~300 bp
- **Quality** - very good


## Primers
- V3-V4 variable regions (about 420bp)
- FWD - 341F - 5’ - CCTACGGGNGGCWGCAG - 3’
- REV -  805R - 5’ - GACTACHVGGGTATCTAATCC - 3’
- forward primer found at the beginning of forward reads (same for reverse reads & reverse primer).

A few samples have some reads with reverse complement of forward/reverse primer (\~500 reads). They were not removed.


## Filtering
- **primers removal** - \~98% reads kept
- **quality filter** - \~72% reads kept => 3 samples had no identifier to allow to match forward & reverse reads (ERR3586005, ERR3586042 and ERR3586312 were removed).


## Learn error rates
- parametric error model fits data

## Construct ASV table
### a) Infer sequence variants
- 16,478 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 5,902 seq variants (but still >97% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - 5,896
- Archaea - 6
- Eukaryota - 0

All unassigned phyla were removed (n=12), samples below 500 total reads (n=79). The final ASV table contains **5,885 sequence variants**.

## Metadata
- age
- gender
- BMI
- psychological disorder
- sample type (sigmoid mucosa vs fecal)
- collection date