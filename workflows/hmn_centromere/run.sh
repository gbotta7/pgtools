#!/bin/bash
#SBATCH --job-name=run_pgtools_cm
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=5G 
#SBATCH --time=200:00:00      
#SBATCH --output run_pgtools_cm.log
#SBATCH --mail-type=END
#SBATCH --mail-user=gianfranco@ds.dfci.harvard.edu

source /homes2/gianfranco/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

# module load htslib

snakemake --unlock
snakemake --use-conda --cores 64 --rerun-incomplete