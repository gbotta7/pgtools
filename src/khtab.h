#ifndef KHTAB_H
#define KHTAB_H

#include <stdint.h>

#include "shtab.h"
#include "utils.h"

// define the k-mer hash table value
#define K_COUNTER_BITS 10
#define K_COUNTER_MAX ((1 << K_COUNTER_BITS)-1)

// define the k-mer info hash table value
#define K_POS_BITS 32
#define K_POS_MAX ((1U << K_POS_BITS) - 1)

#define pg_k_eq(a, b) ((a) >> K_COUNTER_BITS == (b) >> K_COUNTER_BITS)
#define pg_k_hash(a) ((a) >> K_COUNTER_BITS)

typedef struct {
    uint32_t pos;
    uint16_t seq_idx;
} kci_entry_t;

typedef struct {
    kci_entry_t *i;
    int n, m;
} kci_t; // k-mer count and info

KHASHL_SET_INIT(KH_LOCAL, pg_kht_t, pg_kht, uint64_t, pg_k_hash, pg_k_eq)
KHASHL_MAP_INIT(KH_LOCAL, pg_kiht_t, pg_kiht, uint64_t, kci_t, pg_k_hash, pg_k_eq)

struct pg_kht_t;

typedef struct {
	struct pg_kht_t *h;
} pg_kht1_t;

typedef struct {
	struct pg_kiht_t *ih;
} pg_kiht1_t;

typedef struct {
    pthread_mutex_t mutex;
    cnames_t cnames;        // contig names
    pg_kht1_t *h;           // array of partitions (size = 1 << pre)
    pg_kiht1_t *ih;          // pg_kht1_t with info
    int64_t n_ins_tot;      // stores the total k-mer insertions in pangenome
    int64_t n_del_tot;      // stores the total filtered k-mers in the filter stage
    int32_t k;              // k-mers length
    int32_t pre;            // partition bits
} pg_mkht_t;

pg_mkht_t *pg_mkht_init(int k, int pre, int w);
void pg_mkht_destroy(pg_mkht_t *h, int w);
int pg_mkht_insert_list(pg_mkht_t *h, int n, const seq_t *a, seq_info_t *b, int w);
void pg_mkht_tighten(pg_mkht_t *h);

void write_kmer_tsv(const char *out_fn, pg_mkht_t *h, const char *gnm_fn, int w, int w_mko);

#endif // KHTAB_H