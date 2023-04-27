# Differential abundance analysis (Cluster part)

This directory is to compute the differential abundance testing for scCODA on the AGP dataset and the Nagel/Pozuelo shared ASV data.

We recommend you copy this whole directory to a computer cluster to run the scripts. 

**You must also have a python environment set up on the cluster in the same way as for the rest of the analysis!**


**Necessary elements**:
1. **Input tables** in the `input/` subdirectory containing all aggregated data for AGP and the shared ASV data. 
In total, you should have six files here:
`agp_class-agg.csv`, `agp_family-agg.csv`, `agp_genus-agg.csv`, `agp_order-agg.csv`, `agp_phylum-agg.csv`, `commonASV_Nagel-Pozuelo.csv`
3. **Python scripts** to compute differential abundance => `run_AGP_scCODA.py` and `run_common_NagPoz.py` are already provided in this directory, but _make sure to change the file paths at the top of the scripts!!_
4. **Bash scripts** to execute the python scripts on your computer cluster (recommended memory/CPU settings are given in the files!). Make sure to modify the bash scripts for [the analysis of the AGP data](./scCODA_AGP_job.sh) and [the analysis of the common ASV data](./scCODA_AGP_job.sh) provided in this directory to fit your cluster.
