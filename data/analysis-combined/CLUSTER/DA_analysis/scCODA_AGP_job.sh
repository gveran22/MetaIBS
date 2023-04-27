#!/bin/bash
#SBATCH -o /home/CLUSTER/DA_analysis/slurm_out.o
#SBATCH -e /home/CLUSTER/DA_analysis/slurm_error.e
#SBATCH -p cpu_p
#SBATCH -c 8
#SBATCH --mem=500000
#SBATCH --nice=10000
#SBATCH -t 3-00:00:00

# This script was written for a cluster running slurm. Please replace all paths (input, output, data and python installation), partition names, etc. in this script with the ones fitting your cluster configuration
# Also change the data/output paths in the run_AGP_scCODA script!

/home/anaconda3/envs/metaIBS/bin/python /home/CLUSTER/DA_analysis/run_AGP_scCODA.py