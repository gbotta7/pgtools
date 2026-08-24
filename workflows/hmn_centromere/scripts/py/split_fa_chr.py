import sys
import gzip
from Bio import SeqIO

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

fa_in = snakemake.input["fa_file"]
alias_file = snakemake.input["chromalias"]
fa_out = snakemake.output["fa_file"]
chrom_target = snakemake.params["chrom"]

contigs = set()
with open(alias_file) as f:
    for line in f:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        assembly_contig = fields[0].strip()
        ucsc_chr = fields[1].strip()
        chr_prefix = ucsc_chr.split("_")[0]
        if chr_prefix == chrom_target:
            contigs.add(assembly_contig)

print(f"[split_fasta_by_chr] found {len(contigs)} contig(s) for {chrom_target}: {sorted(contigs)}")

if not contigs:
    print(f"[split_fasta_by_chr] WARNING: no contigs matched {chrom_target}. Output will be empty.")

n_written = 0
with gzip.open(fa_in, "rt") as fin, gzip.open(fa_out, "wt") as fout:
    for record in SeqIO.parse(fin, "fasta"):
        if record.id in contigs:
            SeqIO.write(record, fout, "fasta")
            n_written += 1

print(f"[split_fasta_by_chr] wrote {n_written} record(s) to {fa_out}")

if n_written != len(contigs):
    print(
        f"[split_fasta_by_chr] WARNING: expected {len(contigs)} contigs but wrote {n_written}. "
        f"Check that assembly IDs in chromAlias match fasta headers exactly."
    )