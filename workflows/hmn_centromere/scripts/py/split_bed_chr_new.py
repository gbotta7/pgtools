#!/usr/bin/env python3
import os
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

# Suffixes tried, in order, for each haplotype
HAP_SUFFIXES = {
    "1": ["1", "pat"],
    "2": ["2", "mat"],
}

NEAR_EDGE_THRESHOLD = 5_000_000  # 5 Mb


def load_fai_lengths(fai_path):
    lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def load_mashmap_order(mashmap_path, chrom_filt):
    order = {}
    with open(mashmap_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            query_name = fields[0]
            target_name = fields[5]
            target_start = int(fields[7])

            if target_name != chrom_filt:
                continue

            if query_name not in order or target_start < order[query_name]:
                order[query_name] = target_start
    return order


def build_full_order(mashmap_order):
    full_order = sorted(mashmap_order.keys(), key=lambda c: mashmap_order[c])
    contig_to_index = {c: idx for idx, c in enumerate(full_order)}
    return full_order, contig_to_index


def find_neighbor(full_order, i, direction, contig, fai_lengths, min_len=NEAR_EDGE_THRESHOLD):
    """
    Walk from index i in full_order (the complete mashmap-derived contig
    order for this chromosome) in the given direction (+1 or -1), skipping:
      - the same contig itself (defensive - shouldn't recur in full_order)
      - contigs shorter than min_len (too small to be a useful anchor)
    full_order is already restricted to a single chromosome by construction,
    so there's no cross-chromosome boundary to check here.
    Returns the first qualifying contig name, or "" if none found.
    """
    j = i + direction
    while 0 <= j < len(full_order):
        cand_contig = full_order[j]

        if cand_contig != contig:
            cand_len = fai_lengths.get(cand_contig)
            if cand_len is not None and cand_len >= min_len:
                return cand_contig
            # else: unknown length or too small - keep looking past it

        j += direction

    return ""


def process_bed(bed_path, chromalias, out_path, chrom_filt, mashmap_order, fai_lengths):
    missing_contigs = set()
    not_in_mashmap = set()
    missing_length = set()
    header = []
    data_lines = []  # list of (contig, start, end, fields, ucsc)

    full_order, contig_to_index = build_full_order(mashmap_order)

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

            if ucsc != chrom_filt:
                continue

            if contig not in contig_to_index:
                # Contig has no alignment to this chromosome in the mashmap file
                not_in_mashmap.add(contig)
                continue

            start = int(fields[1])
            end = int(fields[2])

            data_lines.append((contig, start, end, fields, ucsc))

    # Preserve output ordering by position in the full chromosome order,
    # rather than by whatever order bed_path happened to list them in.
    data_lines.sort(key=lambda x: contig_to_index[x[0]])

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as out_f:
        for h in header:
            out_f.write(h + "\n")

        for contig, start, end, fields, ucsc in data_lines:
            contig_len = fai_lengths.get(contig)
            preceding_contig = ""
            following_contig = ""

            if contig_len is None:
                missing_length.add(contig)
            else:
                near_start = start < NEAR_EDGE_THRESHOLD
                near_end = (contig_len - end) < NEAR_EDGE_THRESHOLD

                i = contig_to_index[contig]

                if near_start:
                    preceding_contig = find_neighbor(
                        full_order, i, -1, contig, fai_lengths
                    )

                if near_end:
                    following_contig = find_neighbor(
                        full_order, i, +1, contig, fai_lengths
                    )

            out_f.write("\t".join(fields + [ucsc, preceding_contig, following_contig]) + "\n")

    if missing_contigs:
        print(f"Missing contigs (no chromalias entry) in file {bed_path}:\n" +
              "\n".join(sorted(missing_contigs)))
    if not_in_mashmap:
        print(f"Contigs discarded (not found in mashmap for {chrom_filt}) in file {bed_path}:\n" +
              "\n".join(sorted(not_in_mashmap)))
    if missing_length:
        print(f"Contigs missing from fai index in file {bed_path}:\n" +
              "\n".join(sorted(missing_length)))


bed_path = snakemake.input["bed"]
chromalias_path = snakemake.input["chrom_table"]
mashmap_path = snakemake.input["mashmap"]
fai_idx = snakemake.input["fai_idx"]
out_path = snakemake.output["bed"]
chrom_filt = snakemake.params["chrom"]

# Load aliases
chromalias = {}
with open(chromalias_path) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        chromalias[fields[0]] = fields[1].split("_")[0]

# Load contig -> CHM13 target-start ordering from the mashmap file
mashmap_order = load_mashmap_order(mashmap_path, chrom_filt)

# Load contig lengths from fai
fai_lengths = load_fai_lengths(fai_idx)

process_bed(bed_path, chromalias, out_path, chrom_filt, mashmap_order, fai_lengths)