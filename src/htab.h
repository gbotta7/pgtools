#ifndef HTAB_H
#define HTAB_H

#include <pthread.h>
#include <stdint.h>

#define PG_MAGIC "PGKM"
#define COUNTER_BITS 10
#define COUNTER_MAX ((1U << COUNTER_BITS) - 1)

typedef struct { // terminal options
    int64_t chunk_size;
    int64_t n_threads;
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
    uint64_t counts; // stores the total k-mer counts in pangenome.
    int32_t k; // k-mers length.
    int32_t pre; // stores the k-mer flanking sequences, encoded as 2 bits per base (used for partitioning).
} pg_mht_t;

// typedef struct __attribute__((packed)) {
//     uint64_t kmer; // stores the k-mer flanking sequences, encoded as 2 bits per base, and the central base (biallelic SNPs).
//     uint16_t count; // stores the count of k-mer with central base a and central base b.
// } pg_snpmer_t;

pg_mht_t *pg_mht_init(int k, int pre);
void pg_mht_destroy(pg_mht_t *h);
int pg_mht_insert_list(pg_mht_t *h, int n, const uint64_t *a);
// void pg_mht_filter(pg_mht_t *h, int min_count);
void pg_mht_tighten(pg_mht_t *h);
pg_mht_t *pg_count(const char **fns, int n_fns, const pg_opt_t *opt);
int pg_mht_dump(const pg_mht_t *h, const char *fn);

#endif // HTAB_H