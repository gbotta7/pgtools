## User's guide

**pgtools** is a fast, multi-threaded C toolkit for pangenome SNP-mer analysis. Given a set of genome assemblies, it identifies *SNP-mers* — k-mers whose central base varies across genomes — and builds a VCF file recording:
1. The reference and alternate alleles for each SNP, along with the SNP-mers it originated from.
2. Per-genome allele counts for each SNP.
3. AC, AN, and AF fields.
5. Optionally, all hits of each SNP-mer across all input assemblies (contig|strand|position) in the INFO field.

Typical use cases can be:
- **Pangenome SNP calling** — efficiently identify genome-wide SNPs across large collections of bacterial, fungal, or viral assemblies without whole-genome alignment.
- **Population genomics** — characterize allele frequencies and SNP distributions across hundreds to thousands of genomes.

---

## How it works

pgtools operates in two stages:

**Stage 1 — k-mer counting and SNP-mer discovery (`pg_count_k`)**

Each genome is streamed in chunks and its k-mers are extracted and inserted into a partitioned hash table. Every time a number of genomes is processed, low-frequency k-mers are filtered out. Once all genomes are processed, a final filter retains only *SNP-mers* that appear in at least a user-defined fraction of genomes (-m), indicating a single-nucleotide polymorphism at that position (only bi-allelic SNPs are retained).

**Stage 2 — Per-genome SNP-mer counting (`pg_count_snp`)**

Each genome is streamed again in chunks, and only its SNP-mers are inserted into a partitioned hash table, where per-genome allele counts - as well as chromosome name, genomic positions, and strand information (if required by the user) - are tracked. For each genome, this information is written to a genome-specific VCF file. After all genomes have been processed, the individual VCF files are merged into a single final VCF containing all SNPs and their corresponding per-genome information. This second pass is also parallelized across genomes, as the memory required to maintain multiple SNP-mer hash tables is substantially lower. Multiple genomes can therefore be processed concurrently, with each genome typically using 3–4 threads.

---


## Installation

**Dependencies:** `zlib`, `pthreads`, a C compiler supporting C11 (GCC or Clang).

```bash
git clone https://github.com/gbotta7/pgtools.git
cd pgtools/src
make
```

---

## General usage

```bash
pgtools count [options] <in1.fa> [in2.fa [...]]
```

To know more about options:
```bash
pgtools count
```

To identify unique SNP-mers across 90% of the genomes in the mtb152 dataset against the mtb reference genome, using 18 threads:
```bash
pgtools count -k21 -m0.9 -t18 \
        -o mtb152_snpmers.vcf \
        -r mtb152_asm/H37Rv.fa.gz \
        mtb152_asm/*.fa
```

To identify unique SNP-mers across 98% of the genomes in the hmn579 dataset against the reference (added to the assemblies), using 12 threads:
```bash
pgtools count -k31 -m0.98 -t12 \
        -o hmn580_snpmers.vcf \
        -r hmn_ref_asm/GRCh38_genomic.fna.gz \
        hmn_ref_asm/GRCh38_genomic.fna.gz hmn579_asm/*.fa
```

To identify unique SNP-mers across 98% of the genomes in the hmn579 dataset against the reference (added to the assemblies) and keep mapping information across the pangenome, using 12 threads:
```bash
pgtools count -k31 -m0.98 -t12 \
        -o hmn580_snpmers.vcf \
        -r hmn_ref_asm/GRCh38_genomic.fna.gz \
        -w \
        hmn_ref_asm/GRCh38_genomic.fna.gz hmn579_asm/*.fa
```

### Options

| Short      | Long               | Default | Description                                                                            |
| ---------- | ------------------ | ------- | -------------------------------------------------------------------------------------- |
| `-k INT`   | `--kmer INT`       | `31`    | k-mer size (must be odd and ≤ 31)                                                      |
| `-m FLOAT` | `--min-freq FLOAT` | `0.95`  | Minimum fraction of genomes a k-mer must appear in to be retained                      |
| `-p INT`   | `--pre INT`        | `10`    | Number of bits used for hash table partitioning (higher values create more partitions) |
| `-f INT`   | `--filt-type INT`  | `2`     | Filter type for k-mer counting                                                         |
| `-K INT`   | `--chunk-size INT` | `1.9g`  | Input chunk size used for streaming genomes                                            |
| `-t INT`   | `--threads INT`    | `4`     | Number of worker threads                                                               |
| `-w`       | `--write_info`     | off     | Write positions of all pangenome hits to the INFO field of the output VCF              |
| `-r FILE`  | `--ref FILE`       | —       | Reference genome used to define
SNP-mers                                         |
| `-o FILE`  | `--output FILE`    | —       | Output VCF containing genome-specific SNPs and all their hits in the pangenome (if -w is set)                                            |
| `-v`       | `--verbose`        | off     | Enable verbose logging                                                                 |


#### SNP filtering (`-f`)

Controls which SNP-mers are retained:

| Value         | Description                                                                                                  |
| ------------- | ------------------------------------------------------------------------------------------------------------ |
| `2` (default) | Keep only unique SNP-mers.                                                                                   |
| `1`           | Keep SNPs that occur multiple times, provided that each genome contains only one allele at the SNP position. |
| `0`           | Keep all SNPs without filtering.                                                                             |

Use `-f 1` or `-f 0` if you want to retain non-unique SNP-mers.

---

## File overview

| File | Purpose |
|------|---------|
| `main.c` | Entry point and CLI parsing |
| `count.c` | k-mer counting, SNP-mer discovery, and per-genome counting pipelines |
| `htab.c/h` | Partitioned hash table: insert, count, filter |
| `utils.c/h` | k-mer hashing, nucleotide tables, option initialization |
| `kthread.c/h` | Thread pool primitives (`kt_for`, `kt_pipeline`) |
| `parser.c/h` | FASTA/FASTQ streaming via kseq |
| `kseq.h` | Single-header FASTA/FASTQ parser |
| `khashl.h` | Generic open-addressing hash table |
| `ketopt.h` | Lightweight command-line option parser |
| `sys.c/h` | Wall-clock time, CPU time, and peak RSS reporting |
