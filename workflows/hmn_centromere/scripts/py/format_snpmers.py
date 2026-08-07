#!/usr/bin/env python3
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.input["snpmers"]) as fin, open(snakemake.output["snpmers"], "w") as fout:
    for line in fin:
        kmer, ref, alt = line.rstrip("\n").split("\t")[:3]
        mid = len(kmer) // 2
        kmer_ref = kmer[:mid] + ref + kmer[mid + 1:]
        kmer_alt = kmer[:mid] + alt + kmer[mid + 1:]
        fout.write(">%s\n%s\n" % (kmer_ref, kmer_ref))
        fout.write(">%s\n%s\n" % (kmer_alt, kmer_alt))