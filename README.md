## User's guide

**pgkmc** is a fast, multi-threaded C toolkit for pangenome k-mer and SNP-mer analysis. Given a set of genome assemblies, it can count k-mers or identify and count *SNP-mers* — k-mers whose central base varies across genomes. It produces a TSV file with three columns:
1. k-mer or SNP-mer sequence
2. Counts
3. Optionally, the positions of each k-mer or SNP-mer in that sample

Typical use cases can be:
- **Pangenome SNP calling** — efficiently identify genome-wide SNPs across large collections of bacterial, fungal, or viral assemblies without whole-genome alignment and also in repetitive regions.
- **Population genomics** — characterize allele frequencies and SNP distributions across hundreds to thousands of genomes.

---

## How it works

When counting SNP-mers, pgkmc operates in two stages:

**Stage 1 — SNP-mer discovery (`pgkmc detect`)**

Each genome is streamed in chunks and its k-mers are extracted and inserted into a partitioned hash table. Every time a number of genomes is processed, low-frequency k-mers are filtered out. Once all genomes are processed, a final filter retains only *SNP-mers* that satisfy the user's or the default condition (--msf, --maf, --mko), indicating a single-nucleotide polymorphism at that position (only bi-allelic SNPs are retained).

**Stage 2 — Per-genome SNP-mer counting (`pgkmc count`)**

The genome is streamed again in chunks together with a set of SNP-mers, output of pgkmc detect. SNP-mers repopulates the partitioned hash table, per-genome allele counts are performed. Optionally, chromosome name and genomic positions can be tracked. For each genome, this information is written to a genome-specific TSV file. The counts for each file can be merged downstream with a Python code, since the order of SNP-mers is the same across the files (if the SNP-mers file used is always the same).

When counting k-mers, pgkmc operates in one stage:

The genome is streamed in chunks and its k-mers are extracted and inserted into a partitioned hash table, where they get counted.

---


## Installation

**Dependencies:** `zlib`, `pthreads`, a C compiler supporting C11 (GCC or Clang).

```bash
git clone https://github.com/gbotta7/pgkmc.git
cd pgkmc/src
make
```

---

## General usage

```bash
pgkmc detect [options] <in1.fa> [in2.fa [...]]
pgkmc count [options] <in.fa>
```

To know more about options:
```bash
pgkmc detect
pgkmc count
```

To identify unique SNP-mers present in at least 90% of the genomes in the mtb152 dataset, and count them in a specific file, using 12 threads:
```bash
pgkmc detect --snp -k21 -f2 --msf 0.9 -t12 -o mtb152_snpmers.txt mtb152_asm/*.fa;
pgkmc count --snp -k21 -t12 --kmers mtb152_snpmers.txt -o mtb152_snpmers.tsv mtb152_asm/fileN.fa;
```

To additionally keep mapping information of the SNP-mers:
```bash
pgkmc count --snp -w -k21 -t12 --kmers mtb152_snpmers.txt -o mtb152_snpmers.tsv mtb152_asm/fileN.fa;
```

To count k-mers in a specific file, using 3 threads:
```bash
pgkmc count -k21 -t3 -o mtb152_kmers.tsv mtb152_asm/fileN.fa;
```

To identify all centromeric SNP-mers present in at least 50% of the genomes, with MAF > 0.05, in the HPRC dataset, and count them in a specific file, using 3 threads:
```bash
pgkmc detect --snp -k31 -f0 --msf 0.5 --maf 0.05 -t3 \ 
        -b hprc_cent_anno_files.txt \ 
        -o hprc_cent_snpmers.txt \ 
        hprc_asm/*.fa;
pgkmc count --snp -k31 -t12 --kmers hprc_cent_snpmers.txt \ 
        -b hprc_cent_anno/fileN.bed \ 
        -o hprc_cent_snpmers.tsv hprc_asm/fileN.fa;
```    


#### SNP filtering (`-f/--filt_type`)

Controls which SNP-mers are retained:

| Value         | Description                                                                                                  |
| ------------- | ------------------------------------------------------------------------------------------------------------ |
| `2` (default) | Keep only unique SNP-mers.                                                                                   |
| `1`           | Keep SNPs that occur multiple times, provided that each genome contains only one allele at the SNP positions. |
| `0`           | Keep all SNPs without filtering.                                                                             |

Use `-f 1` or `-f 0` if you want to retain non-unique SNP-mers.

#### Restricting to specific regions (`-b`)

`-b/--bed_list` takes a single text file listing the paths of all the BED files in pgkmc detect, one per fasta file passed. Therefore, the list must have exactly as many entries as there are input files. In pgkmc count, `-b/--bed` takes a single bed file with the regions of interest in the input fasta file.

---

## File overview

| File | Purpose |
|------|---------|
| `main.c` | Entry point and CLI parsing |
| `bed.c/h` | BED files handling |
| `count.c` | K-mer and SNP-mer counting, SNP-mer discovery |
| `khtab.c/h` | Partitioned k-mer hash table: insert, count, filter |
| `shtab.c/h` | Partitioned SNP-mer hash table: insert, count, filter |
| `utils.c/h` | K-mer hashing, nucleotide tables, option initialization |
| `kthread.c/h` | Thread pool primitives (`kt_for`, `kt_pipeline`) |
| `parser.c/h` | FASTA/FASTQ streaming via kseq |
| `kseq.h` | Single-header FASTA/FASTQ parser |
| `khashl.h` | Generic open-addressing hash table |
| `ketopt.h` | Lightweight command-line option parser |
| `sys.c/h` | Wall-clock time, CPU time, and peak RSS reporting |
