source /homes2/gianfranco/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --use-conda --cores 4 --rerun-incomplete 2>&1 | tee snakemake_$(date +%Y%m%d_%H%M%S).log