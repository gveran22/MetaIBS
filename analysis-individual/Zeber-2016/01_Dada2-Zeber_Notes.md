# DADA2 - ZEBER DATASET

Zeber-Lubecka et al. (_Gut Microbes_, 2016) - [Limited prolonged effects of rifaximin treatment on irritable bowel syndrome-related differences in the fecal microbiome and metabolome][1]

[1]: https://www.tandfonline.com/doi/full/10.1080/19490976.2016.1215805


## Samples
- **90 total**
- 17 Healthy & 73 IBS

## Data Quality
- **Technology** - Ion Torrent, single end
- **Nb of reads per sample** - mean of 108,088 reads per sample (50,670 - 294,980)
- **Read length** - ~200 bp
- **Quality** - average (Q around 30)

## Primers
- Ion 16S Metagenomics Kit (V2-4-8, V3-6, V7-9)
- Primers are not known and reads from all variable regions are pooled together

## Filtering
- **primers removal** - because several variable regions with unknown primers were sequenced, reads were filtered following recommendations from QIIME forum: https://forum.qiime2.org/t/possible-analysis-pipeline-for-ion-torrent-16s-metagenomics-kit-data-in-qiime2/13476
- **quality filter** - \~50% reads kept per sample.

## Learn error rates
- parametric error model was manually modified to fit better the data

## Construct ASV table
### a) Infer sequence variants
- 15,952 amplicon sequence variants (ASVs)

### b) Remove chimeras
- 14,842 seq variants (>96% reads kept)

### c) Assign taxonomy
Taxonomy assigned with Silva v138.
- Bacteria - 14,834
- Archaea - 0
- Eukaryota - 5

All unassigned phyla were removed (n=30), samples below 500 total reads (n=0). The final ASV table contains **14,812 sequence variants**.

## Metadata
- IBS subtype
- gender