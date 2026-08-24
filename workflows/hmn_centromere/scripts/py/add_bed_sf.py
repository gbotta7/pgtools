#!/usr/bin/env python3
import os
import re
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

def subtype_to_superfamily(name):
    m = re.search(r"\((.+)\)$", name)
    inner = m.group(1) if m else name
    first = inner.split(",")[0]
    base = first.split(".")[0]

    if inner in sf_map:
        return sf_map[inner]
    if first in sf_map:
        return sf_map[first]
    if base in sf_map:
        return sf_map[base]

    match = re.match(r"S(\d+)", base)
    if match:
        return "SF" + match.group(1)
    return name

inbed = snakemake.input["bed"]
sf_file = snakemake.input["sf"]
outbed = snakemake.output["bed"]

sf_map = {}
with open(sf_file) as fh:
    next(fh)
    for line in fh:
        fields = line.rstrip("\n").split("\t")
        sf_map[fields[0]] = fields[1]

with open(inbed) as fin, \
        open(outbed, "w") as fout:
    for line in fin:
        fields = line.rstrip("\n").split("\t")
        fields.append(subtype_to_superfamily(fields[3]))
        fout.write("\t".join(fields) + "\n")