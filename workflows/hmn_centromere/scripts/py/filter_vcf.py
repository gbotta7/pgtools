import argparse
import gzip
import re
import sys

PASSED_COL = 4
CHROM1_COL = 6
IDX_COL = 0


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def contig_sort_key(name):
    m = re.fullmatch(r"(?:chr)?(\d+)", name)
    if m:
        return (0, int(m.group(1)), "")
    m = re.fullmatch(r"(?:chr)?([XYM]|MT)", name, re.IGNORECASE)
    if m:
        order = {"X": 0, "Y": 1, "M": 2, "MT": 2}
        return (1, order[m.group(1).upper()], "")
    return (2, 0, name)

def load_index(path):
    keep = {}
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n").rstrip("\r")
            if not line or line.startswith("#") or lineno == 1:
                continue
            f = line.split("\t")
            if f[PASSED_COL].strip().upper() != "TRUE":
                continue
            idx = int(f[IDX_COL])
            chrom1 = f[CHROM1_COL].strip()
            keep[idx] = chrom1
    return keep

# Parse args
ap = argparse.ArgumentParser()
ap.add_argument("vcf")
ap.add_argument("info")
ap.add_argument("-o", "--output", default="-", help="output path (default: stdout)")

args = ap.parse_args()

keep = load_index(args.info)

header = []
records = []
n = 0

with opener(args.vcf) as fh:
    for line in fh:
        if line.startswith("##"):
            if not line.startswith("##contig="):
                header.append(line.rstrip("\n"))
            continue
        if line.startswith("#CHROM"):
            chrom_line = line.rstrip("\n")
            continue
        n += 1
        if n not in keep:
            continue
        cols = line.rstrip("\n").split("\t")
        cols[0] = keep[n]
        records.append((int(cols[1]), cols))

contigs = sorted({c[1][0] for c in records}, key=contig_sort_key)
rank = {c: i for i, c in enumerate(contigs)}
records.sort(key=lambda r: (rank[r[1][0]], r[0]))

out = sys.stdout if args.output == "-" else open(args.output, "w")
for h in header:
    out.write(h + "\n")
for c in contigs:
    out.write(f"##contig=<ID={c}>\n")
out.write(chrom_line + "\n")
for _, cols in records:
    out.write("\t".join(cols) + "\n")

if out is not sys.stdout:
    out.close()