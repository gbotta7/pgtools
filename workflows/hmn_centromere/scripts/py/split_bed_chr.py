#!/usr/bin/env python3
import os
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

# Suffixes tried, in order, for each haplotype
HAP_SUFFIXES = {
    "1": ["1", "pat"],
    "2": ["2", "mat"],
}


def chrom_sort_key(name):
    core = name[3:] if name.startswith("chr") else name
    return (0, int(core), "") if core.isdigit() else (1, 0, name)


def process_bed(bed_path, chromalias, out_path, chrom_filt):
    missing_contigs = set()
    header = []
    data_lines = []

    with open(bed_path) as fin:
        for line in fin:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                header.append(line)
                continue

            fields = line.split("\t")
            contig = fields[0]
            if contig.startswith("chr"):
                ucsc = contig
            else:
                ucsc = chromalias.get(contig)

            if ucsc is None:
                missing_contigs.add(contig)
                continue

            if ucsc == chrom_filt:
                data_lines.append("\t".join(fields + [ucsc]))

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as out_f:
        for h in header:
            out_f.write(h + "\n")
        for d in data_lines:
            out_f.write(d + "\n")

    if missing_contigs:
        print(f"Missing contigs in file {bed_path}:\n" + "\n".join(sorted(missing_contigs)))


bed_path = snakemake.input["bed"]
chromalias_path = snakemake.input["chrom_table"]
out_path = snakemake.output["bed"]

# Load aliases
chromalias = {}
with open(chromalias_path) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        chromalias[fields[0]] = fields[1].split("_")[0]

process_bed(bed_path, chromalias, out_path, snakemake.params["chrom"])