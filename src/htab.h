#ifndef HTAB_H
#define HTAB_H

#include <pthread.h>
#include <stdint.h>
#include <stdbool.h>

#define PG_MAGIC "PGKM"
#define VAL_INFO_BITS 6
#define COUNTER_BITS (32 - VAL_INFO_BITS)
#define COUNTER_MAX ((1U << COUNTER_BITS/2) - 1)

#define val_count1(v) (((v) >> VAL_INFO_BITS) & ((1U << COUNTER_BITS/2) - 1))
#define val_count2(v) ((v) >> (VAL_INFO_BITS + COUNTER_BITS/2))
#define val_cb1(v) ((v) & 0x3U)
#define val_cb2(v) (((v) >> 2) & 0x3U)
#define val_snp1(v) (((v) >> 4) & 0x1U)
#define val_snp2(v) (((v) >> 5) & 0x1U)
#define val_pack(cnt2, cnt1, snp2, snp1, cb2, cb1) (((cnt2) << (VAL_INFO_BITS + COUNTER_BITS/2)) | ((cnt1) << VAL_INFO_BITS) | ((snp2) << 5) | ((snp1) << 4) | ((cb2) << 2) | (cb1))

typedef struct __attribute__((packed)) {
	uint64_t h_flanks;
	uint8_t cb; // stores central base of k-mer and SNP information
} ch_seq_t;

typedef struct { // terminal options
    int64_t chunk_size;
    int64_t n_threads;
    double min_freq;
	int32_t k;
    int32_t pre; // number of bits for partitioning.
    int verbose;
} pg_opt_t;

struct pg_ht_t; // see khashl.h and htab.c for the definition of pg_ht_t.

typedef struct { 
    struct pg_ht_t *h; // hash table for each bucket.
    pthread_mutex_t lock;
} pg_ht1_t;

typedef struct {
    pg_ht1_t *h; // array of partitions (size = 1 << pre).
    int n_ins_tot; // stores the total k-mer insertions in pangenome.
    int32_t k; // k-mers length.
    int32_t pre; // stores the k-mer flanking sequences, encoded as 2 bits per base (used for partitioning).
} pg_mht_t;

typedef struct {
    uint64_t flanks;   // hashed flanks (without genome ID bits)
    uint8_t cb1;
    uint8_t cb2;
} pg_snp_t;

typedef struct {
    int n_rows;
    int n_cols;

    pg_snp_t *snps;    // length n_snps
    uint64_t *mat;     // dimensions: (2*n_snps) x n_genomes
} pg_mtx_t;


pg_mht_t *pg_mht_init(int k, int pre);
void pg_mht_destroy(pg_mht_t *h);
int pg_mht_insert_list(pg_mht_t *h, int n, const ch_seq_t *a, uint64_t gnm_id);
int pg_mht_filter(pg_mht_t *h, int n_proc, int n_tot, double min_freq);
void pg_mht_tighten(pg_mht_t *h);
// void pg_mht_dump(const pg_mht_t *h, const char *fn, const char **fns);
pg_mtx_t *pg_mht2mtx(const pg_mht_t *h, int n_fns, int n_kmers);
void pg_mtx_dump(const char *fn, const pg_mht_t *h, const char **fns, const pg_mtx_t *m);

pg_mht_t *pg_count(const char **fns, int n_fns, const pg_opt_t *opt, int filt, pg_mht_t *h_init);

#endif // HTAB_H