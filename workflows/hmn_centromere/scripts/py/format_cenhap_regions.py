import pandas as pd
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

# Find the real header line, skipping any title/blank lines above it
with open(snakemake.input["tds"]) as f:
    lines = f.readlines()

header_idx = next(i for i, line in enumerate(lines) if line.strip().startswith("chromosome"))

df = pd.read_csv(
    snakemake.input["tds"],
    sep=r"\s+",
    na_values="NA",
    skiprows=header_idx
)

bed_rows = []
for _, row in df.iterrows():
    chrom = row["chromosome"]

    if pd.notna(row["p_c"]) and pd.notna(row["p_end"]):
        bed_rows.append((chrom, int(row["p_c"]), int(row["p_end"])))

    if pd.notna(row["q_begin"]) and pd.notna(row["q_c"]):
        bed_rows.append((chrom, int(row["q_begin"]), int(row["q_c"])))

bed_df = pd.DataFrame(bed_rows)
bed_df.to_csv(snakemake.output["bed"], sep="\t", header=False, index=False)