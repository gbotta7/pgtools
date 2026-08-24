import numpy as np
import gzip
import glob
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

def chrom_sort_key(c):
    name = c.replace("chr", "")
    if name.isdigit():
        return (0, int(name))
    order = {"X": 23, "Y": 24, "M": 25}
    return (1, order.get(name, 99))

def _open_vcf(vcf_path):
    """Transparently open a plain-text or bgzip/gzip-compressed VCF."""
    if vcf_path.endswith(".gz") or vcf_path.endswith(".bgz"):
        return gzip.open(vcf_path, "rt")
    return open(vcf_path, "rt")
    
def load_kc_matrix(vcf_path):
    """
    Parse a VCF with a GT:KC FORMAT field into a 3D numpy array.
    Returns kc[i, j, 0] = first allele count, kc[i, j, 1] = second allele
    count, for variant i, sample j.
    """
    rows = []
    variant_info = []
    sample_ids = None

    with _open_vcf(vcf_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_ids = line.lstrip("#").split("\t")[9:]
                continue

            fields = line.split("\t")
            chrom, pos, vid, ref, alt = fields[0:5]
            fmt = fields[8].split(":")
            kc_idx = fmt.index("KC")

            sample_fields = fields[9:]
            counts = np.zeros((len(sample_fields), 2), dtype=np.uint16)
            for j, sf in enumerate(sample_fields):
                kc_str = sf.split(":")[kc_idx]
                a1, a2 = kc_str.split(",")
                counts[j, 0] = int(a1)
                counts[j, 1] = int(a2)

            rows.append(counts)
            variant_info.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt})

    kc = np.stack(rows, axis=0)   # shape: (n_variants, n_samples, 2)
    return kc, variant_info, sample_ids

def sum_bed_lengths(bed_dir: str):
    bed_files = sorted(glob.glob(os.path.join(bed_dir, "*.bed")))

    per_file_totals = {}
    grand_total = 0

    for bed_path in bed_files:
        file_total = 0
        n_rows = 0
        with open(bed_path) as f:
            for line_num, line in enumerate(f, 1):
                line = line.rstrip("\n")
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.split("\t")
                start = int(fields[1])
                end = int(fields[2])
                
                file_total += (end - start)
                n_rows += 1

        per_file_totals[bed_path] = (file_total, n_rows)
        grand_total += file_total

    return grand_total, per_file_totals

n_snps_list = {}
as_avg_lengths_list = {}

for vcf, bed_dir in zip(snakemake.input["vcf"], snakemake.input["bed_dir"]):
    chrom = os.path.basename(os.path.dirname(vcf))

    # Count number of SNPs
    kc, variant_info, sample_ids = load_kc_matrix(vcf)
    n_snps_list[chrom] = kc.shape[0]

    # Count the lenght of each centromere
    grand_total, per_file_totals = sum_bed_lengths(bed_dir)
    as_avg_lengths_list[chrom] = grand_total/len(per_file_totals)


# Plot stats
df = pd.DataFrame({
    "chrom": list(n_snps_list.keys()),
    "n_snps": list(n_snps_list.values()),
    "as_avg_length": [as_avg_lengths_list[c] for c in n_snps_list.keys()],
})
# SNP density = average centromere length per SNP (bp/SNP)
df["snp_density"] = df["as_avg_length"] / df["n_snps"]
df = df.sort_values("chrom", key=lambda col: col.map(chrom_sort_key)).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(df["chrom"], df["snp_density"], color="#4C72B0")
ax.set_xlabel("Chromosome")
ax.set_ylabel("SNP density (bp / SNP)")
ax.set_title("Centromeric a-satellite SNP density by chromosome")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(snakemake.output["plot"], dpi=300)
plt.show()