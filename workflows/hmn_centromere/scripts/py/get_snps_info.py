#!/usr/bin/env python3
import csv
import gzip
import pandas as pd
import pysam
import sys

from collections import Counter
from functools import lru_cache

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

chromalias_path = snakemake.input["chrom_table"]
bed_path = snakemake.input["anno_table"]  # bgzipped + tabix-indexed
hits_path = snakemake.input["hits"]

KMER_LENGTH = snakemake.params["k"]
THRESH = snakemake.params["threshold"]
DISCARD = {"chrUn", "chrM", "chrX", "chrY"}
BED_COLS = ["chrom", "start", "end", "array_type", "score", "strand", "start1", "end1", "rgb", "sf"]


def parse(line):
    p = line.rstrip("\n").split("\t")
    return p[0], int(p[3]), int(p[4]), p[5:]


SF_COL = BED_COLS.index("sf")


@lru_cache(maxsize=1_000_000)
def lookup_sf(contig, pos, k=KMER_LENGTH):
    """SF of the bed interval with the largest overlap of [pos, pos+k), or None."""
    q_start, q_end = pos, pos + k
    try:
        rows = BED.fetch(contig, q_start, q_end, parser=pysam.asTuple())
    except ValueError:
        # contig absent from the tabix index
        MISSING_CONTIGS.add(contig)
        return None

    best_sf, best_ov = None, 0
    for row in rows:
        ov = min(int(row[2]), q_end) - max(int(row[1]), q_start)
        if ov > best_ov:
            best_sf, best_ov = row[SF_COL], ov
    return best_sf


def get_info(hits, mapping):
    c_out, sf_out, unmapped = [], [], 0
    for h in hits:
        # strip trailing :strand:pos
        contig, _strand, pos = h.rsplit(":", 2)
        c = mapping.get(contig)
        if c is None:
            unmapped += 1
            continue
        c_out.append(c)
        # alpha satellite superfamily
        sf_out.append(lookup_sf(contig, int(pos)))

    return c_out, sf_out, unmapped


def get_top_item(vec):
    vec = [v for v in vec if v is not None]
    if not vec:
        return None, 0.0
    c, n = Counter(vec).most_common(1)[0]
    return c, n / len(vec)


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else -1


def pairs(fh):
    """Yield consecutive non-blank, non-comment line pairs."""
    it = (ln for ln in fh if ln.strip() and not ln.startswith("#"))
    return zip(it, it)


with open(chromalias_path, newline="") as fh:
    reader = csv.reader(fh, delimiter="\t")
    next(reader)  # header
    assembly_to_ucsc = {row[0]: row[1].split("_")[0] for row in reader}

MISSING_CONTIGS = set()
BED = pysam.TabixFile(bed_path)
print(f"bed index: {len(BED.contigs)} contigs")

kept, records = [], []

with gzip.open(hits_path, "rt") as fh:
    for idx, (l1, l2) in enumerate(pairs(fh), start=1):
        k1, tot1, lim1, h1 = parse(l1)
        k2, tot2, lim2, h2 = parse(l2)

        hd = hamming(k1, k2)
        if hd != 1:
            raise ValueError(
                f"pair {idx}: expected 1 substitution, got {hd}\n  {k1}\n  {k2}"
            )

        c1, sf1, u1 = get_info(h1, assembly_to_ucsc)
        c2, sf2, u2 = get_info(h2, assembly_to_ucsc)
        top_c1, f_c1 = get_top_item(c1)
        top_c2, f_c2 = get_top_item(c2)
        top_sf1, f_sf1 = get_top_item(sf1)
        top_sf2, f_sf2 = get_top_item(sf2)

        if top_c1 is None or top_c2 is None:
            reason = "no_hits"
        elif top_c1 in DISCARD or top_c2 in DISCARD:
            reason = f"discarded_{top_c1}"
        elif f_c1 < THRESH or f_c2 < THRESH:
            reason = "no_full_agreement"
        elif top_c1 != top_c2:
            reason = "chrom_mismatch"
        else:
            reason = None

        ok = reason is None
        if ok:
            kept.append(idx)

        records.append(dict(
            idx=idx, kmer1=k1, kmer2=k2,
            chrom=top_c1 if ok else None, passed=ok, reason=reason,
            chrom1=top_c1, c_frac1=round(f_c1, 4), n1=len(c1), tot1=tot1,
            sf1=top_sf1, sf_frac1=round(f_sf1, 4),
            chrom2=top_c2, c_frac2=round(f_c2, 4), n2=len(c2), tot2=tot2,
            sf2=top_sf2, sf_frac2=round(f_sf2, 4),
            unmapped=u1 + u2,
        ))

        if not idx % 1000:
            print(f"Extracted info for {idx} SNPmers")

BED.close()

qc = pd.DataFrame(records)
qc.to_csv(snakemake.output["info"], sep="\t", index=False)

print(f"\nsf lookup cache: {lookup_sf.cache_info()}")
if MISSING_CONTIGS:
    print(f"contigs absent from bed index: {len(MISSING_CONTIGS)}")
    for c in sorted(MISSING_CONTIGS)[:20]:
        print(f"{c}")

print(f"pairs analysed : {len(qc)}")
print(f"pairs kept     : {len(kept)}")
print(f"unmapped hits  : {qc.unmapped.sum()}")
print("\nfailure reasons:")
print(qc.loc[~qc.passed, "reason"].value_counts().to_string())
print("\nkept per chromosome:")
print(qc.loc[qc.passed, "chrom"].value_counts().to_string())