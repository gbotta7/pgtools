source /homes2/gianfranco/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --use-conda --cores 10 --rerun-incomplete