# Dataset
Hugerth (_Gut Microbiota_, 2019) - [No distinct microbiome signature of irritable bowel syndrome found in a Swedish random population][1]

[1]: https://gut.bmj.com/content/69/6/1076


-------

## Samples

XXXX PART TO MODIFY XXX (so it's more clear - add the drawing from keynote)


Description in paper:
  - **561 samples**
  - 376 sigmoid mucosa samples (63 IBS - 313 HC)
  - 185 fecal samples (32 IBS - 153 HC)
  - individuals who provided biopsy+feces : 29 IBS + 149 HC

In the supplementary data from the paper, they provide the subject_id of included samples. Cross-referencing these subject_ids with the metadata table from the SRA database:
  - **561 samples**
  - 376 sigmoid mucosa samples (62 IBS - 217 HC - 6 unlabelled - 11 _C.difficile_ infection - 2 _C.difficile_ infection+polyps - 3 _C.difficile_ infection+IBS - 59 polyps - 15 polyps+IBS - 1 Salmonella+Campylobacter+IBS)
  - 185 fecal samples (35 IBS - 107 HC - 4 unlabelled - 33 polyps - 6 polyps+IBS)

------

From ENA or SRA database (without filtering the annotations & before quality control):
  - **607 samples**
  - 389 sigmoid mucosa samples (83 IBS - 306 HC)
  - 218 fecal samples (57 IBS - 161 HC)
  - individuals who provided biopsy+feces : 35 IBS + 112 HC


From ENA database (after modifying the annotations and removing inconsistencies as described below):
   - **585 samples**
   - 389 sigmoid mucosa samples (83 IBS - 306 HC)
   - 196 fecal samples (50 IBS - 146 HC)
  - individuals who provided biopsy+feces : 41 IBS + 127 HC


### Modifications done on metadata
Firstly, samples carrying the **same subject_id** were checked. If there were 2 or more inconsistencies, samples would be given different subject_id (different individuals). If there were just 1 inconsistency in the psychological disorder (usually a missing value _NA_), they would keep the same subject_id (see diagram below).

Then, samples were re-ordered by ascending BMI and age, to identify samples with 100% same demographic characteristics but carrying **different subject_id**. There was only 1 case of 2 samples being almost identical, except for one missing value _NA_ in the psychological disorder, but they were considered as belonging to the same patient anyway (see diagram below).

<p align="center">
    <img src="../Metadata_Strategy.jpg" width="1000" title="Strategy to check for inconsistencies in metadata">
</p>



XXXX END PART TO MODIFY XXXXX


## Data Quality
- **Technology** - Illumina MiSeq
- **Nb of reads per sample** - mean of 32,663 reads per sample (9 min - 306,468 max)
- **Read length** - ± 300 bp
- **Quality** - very good (forward reads have a drop in quality in the first 20 bp, then high quality until the last 20 bp // reverse reads very high quality until the last ~50bp).
- **Per base sequence content** - the first bases are always the same and correspond to the primer sequence (FWD primer in forward reads, REV primer sequence in reverse reads)


## Primers
- **V3-V4** variable regions
- FWD - 341F - 5’ - CCTACGGGNGGCWGCAG - 3’
- REV -  805R - 5’ - GACTACHVGGGTATCTAATCC - 3’
- **FWD primer is found at the beginning of forward reads**. **REV primer is found at the beginning of reverse reads**.
- After removing primers, a few samples have some reads with reverse complement of forward/reverse primer (~500 reads). They were not removed.

## Filtering
- primers removal - ~99% reads kept
- quality filter - ~48% reads kept => 3 samples had no identifier to allow to match forward & reverse reads (ERR3586004, ERR3586042 and ERR3586312 were removed)

## Learn error rates
- parametric error model fits data

## Construct ASV table
### a) Infer sequence variants
- 16,443 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 6,033 seq variants (but still >97% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - 6,033
- Archaea - 0
- Eukaryota - 0

All unassigned phyla were removed (n=14), samples below 500 total reads (n=85). The final ASV table contains **5,294 sequence variants**.

## Metadata
- age
- gender
- BMI
- psychological disorder
- sample type (sigmoid mucosa vs fecal)
- collection date