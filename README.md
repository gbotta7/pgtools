# User's guide

**pgtools** is a fast, multi-threaded C toolkit for pangenome SNP-mer analysis. Given a set of genome assemblies (FASTA files), it identifies *SNP-mers* — k-mers whose central base varies across genomes — and builds a compact sparse matrix recording how often each SNP-mer allele appears in each genome.
Typical use cases can be:
- **Pangenome SNP calling** — efficiently identify genome-wide SNPs across large collections of bacterial, fungal, or viral assemblies without whole-genome alignment.
- **Population genomics** — characterize allele frequencies and SNP distributions across hundreds to thousands of genomes.

---

## How it works

pgtools operates in two stages:

**Stage 1 — k-mer counting and SNP-mer discovery (`pg_count`)**

Each genome is streamed in chunks and its k-mers are extracted and inserted into a partitioned hash table. Every time a number of genomes is processed, low-frequency k-mers are filtered out. Once all genomes are processed, a final filter retains only *SNP-mers* that appear in at least a user-defined fraction of genomes, indicating a single-nucleotide polymorphism at that position (only bi-allelic SNPs are retained).

**Stage 2 — Per-genome SNP-mer counting (`pg_findsnp`)**

The SNP-mer hash table is indexed and converted into a sparse matrix (CSR-like format). Each genome is then rescanned in parallel, and the occurrence counts of each SNP-mer allele per genome are stored. The result is written to a tab-separated output file.

---


## Installation

**Dependencies:** `zlib`, `pthreads`, a C compiler supporting C11 (GCC or Clang).

```bash
git clone <repo>
cd pgtools
make
```

---

## General usage

```
pgtools count [options] <in1.fa> [in2.fa [...]]
```

To identify SNP-mers across N bacterial assemblies, using 16 threads
```bash
pgtools count -k 31 -m 0.9 -t 16 -o snpmers.tsv genome_*.fa
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-k INT` | 31 | k-mer size (must be odd and ≤ 31) |
| `-m FLOAT` | 0.95 | Minimum fraction of genomes a k-mer must appear in to be retained |
| `-p INT` | 10 | Number of bits used for hash table partitioning (higher = more partitions) |
| `-f INT` | 2 | Filter type for intermediate k-mer pruning |
| `-t INT` | 4 | Number of worker threads |
| `-K INT` | 1.9g | Chunk size for streaming input |
| `-v` | off | Verbose logging |
| `-o FILE` | — | Output file (TSV) for per-genome SNP-mer counts |

---

## Output format

The output TSV (specified with `-o`) contains one row per SNP-mer observation, with columns for the SNP-mer ID, genome ID, and allele counts for each of the two observed central bases.

---

## File overview

| File | Purpose |
|------|---------|
| `main.c` | Entry point and CLI parsing |
| `count.c` | k-mer counting, SNP-mer discovery, and per-genome counting pipelines |
| `htab.c/h` | Partitioned hash table: insert, count, filter, index, serialize |
| `utils.c/h` | k-mer hashing, nucleotide tables, option initialization |
| `kthread.c/h` | Thread pool primitives (`kt_for`, `kt_pipeline`) |
| `parser.c/h` | FASTA/FASTQ streaming via kseq |
| `kseq.h` | Single-header FASTA/FASTQ parser |
| `khashl.h` | Generic open-addressing hash table |
| `ketopt.h` | Lightweight command-line option parser |
| `sys.c/h` | Wall-clock time, CPU time, and peak RSS reporting |
