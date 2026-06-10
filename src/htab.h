#ifndef HTAB_H
#define HTAB_H

#include <pthread.h>
#include <stdint.h>
#include <stdbool.h>

#define PG_MAGIC "PGKM"

// define the kmer hash table value
#define K_VAL_INFO_BITS 8
#define K_PGNM_COUNTER_BITS 20
#define K_GNM_COUNTER_BITS (32 - (K_VAL_INFO_BITS + K_PGNM_COUNTER_BITS))
#define K_COUNTER_MAX ((1U << K_GNM_COUNTER_BITS/2) - 1)

#define k_val_gnm_count2(v) ((v) >> ((K_VAL_INFO_BITS + K_PGNM_COUNTER_BITS) + K_GNM_COUNTER_BITS/2))
#define k_val_gnm_count1(v) (((v) >> (K_VAL_INFO_BITS + K_PGNM_COUNTER_BITS)) & ((1U << K_GNM_COUNTER_BITS/2) - 1))
#define k_val_filt(v) (((v) >> (6 + K_PGNM_COUNTER_BITS)) & 0x1U)
#define k_val_snp2(v) (((v) >> (5 + K_PGNM_COUNTER_BITS)) & 0x1U)
#define k_val_snp1(v) (((v) >> (4 + K_PGNM_COUNTER_BITS)) & 0x1U)
#define k_val_cb2(v) (((v) >> (2 + K_PGNM_COUNTER_BITS)) & 0x3U)
#define k_val_cb1(v) (((v) >> K_PGNM_COUNTER_BITS) & 0x3U)
#define k_val_pgnm_count2(v) (((v) >> (K_PGNM_COUNTER_BITS/2)) & ((1U << K_PGNM_COUNTER_BITS/2) - 1)) 
#define k_val_pgnm_count1(v) ((v) & ((1U << K_PGNM_COUNTER_BITS/2) - 1))

#define k_val_pack(cnt2, cnt1, filt, snp2, snp1, cb2, cb1, pgnm_cnt2, pgnm_cnt1) (((cnt2) << (K_VAL_INFO_BITS + K_PGNM_COUNTER_BITS + K_GNM_COUNTER_BITS/2))| ((cnt1) << (K_VAL_INFO_BITS + K_PGNM_COUNTER_BITS)) | ((filt) << (K_PGNM_COUNTER_BITS + 6)) | ((snp2) << (K_PGNM_COUNTER_BITS + 5)) | ((snp1) << (K_PGNM_COUNTER_BITS + 4)) | ((cb2) << (K_PGNM_COUNTER_BITS + 2)) | ((cb1) << K_PGNM_COUNTER_BITS) | ((pgnm_cnt2) << (K_PGNM_COUNTER_BITS/2)) | (pgnm_cnt1))

// define the snpmer hash table value
#define S_VAL_INFO_BITS 6
#define S_COUNTER_BITS (32 - S_VAL_INFO_BITS)
#define S_COUNTER_MAX ((1U << S_COUNTER_BITS/2) - 1)

#define s_val_count2(v) ((v) >> (S_VAL_INFO_BITS + S_COUNTER_BITS/2))
#define s_val_count1(v) (((v) >> S_VAL_INFO_BITS) & ((1U << S_COUNTER_BITS/2) - 1))
#define s_val_snp2(v) (((v) >> 5) & 0x1U)
#define s_val_snp1(v) (((v) >> 4) & 0x1U)
#define s_val_cb2(v) (((v) >> 2) & 0x3U)
#define s_val_cb1(v) ((v) & 0x3U)

#define s_val_pack(cnt2, cnt1, snp2, snp1, cb2, cb1) (((cnt2) << (S_VAL_INFO_BITS + S_COUNTER_BITS/2))| ((cnt1) << S_VAL_INFO_BITS) | ((snp2) << 5) | ((snp1) << 4) | ((cb2) << 2) | (cb1));

#define M_COUNTER_BITS 16
#define M_COUNTER_MAX ((1U << M_COUNTER_BITS) - 1)

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
    int filt_type;
    int verbose;
} pg_opt_t;

struct pg_ht_t; // see khashl.h and htab.c for the definition of pg_ht_t.

typedef struct { 
    struct pg_ht_t *h; // hash table for each bucket.
} pg_ht1_t;

typedef struct {
    pg_ht1_t *h;            // array of partitions (size = 1 << pre).
    int64_t n_ins_tot;      // stores the total k-mer insertions in pangenome.
    int64_t n_del_tot;      // stores the total filtered k-mers in the filter stage
    int n_done;             // number of processed files
    int32_t k;              // k-mers length.
    int32_t pre;            // stores the k-mer flanking sequences, encoded as 2 bits per base (used for partitioning).
} pg_mht_t;

typedef struct {
    uint64_t flanks;   // hashed flanks (without genome ID bits)
    uint8_t cb1;
    uint8_t cb2;
} pg_snp_t;

// typedef struct {
//     int n_rows;
//     int n_cols;

//     pg_snp_t *snps;    // length n_snps
//     uint64_t *mat;     // dimensions: (2*n_snps) x n_genomes
// } pg_mtx_t;

typedef struct __attribute__((packed)) {
    uint32_t row_id;
    uint16_t col_id;
    uint16_t cnt1;
    uint16_t cnt2;
} pg_entry_t;

typedef struct {
    pthread_mutex_t mutex;  // mutex to count SNPmers in multiple genomes concurrently (when pldat_t has snp=1)
    pg_entry_t *entries;
    pg_snp_t *snpmer;
    int64_t n, m;
    uint32_t n_snps;
    uint32_t n_fns;
} pg_csr_t;

typedef struct {
    uint32_t *ids;
    uint32_t  n;
} pg_id_map_t;


pg_mht_t *pg_mht_init(int k, int pre);
pg_mht_t *pg_mht_copy(const pg_mht_t *src);
void pg_mht_destroy(pg_mht_t *h);
int64_t pg_mht_insert_list(pg_mht_t *h, int n, const ch_seq_t *a, int f);
void pg_mht_count_list(pg_mht_t *h, int n, const ch_seq_t *a);
void pg_mht_clear_k(pg_mht_t *h, long i, int f);
void pg_mht_clear_s(pg_mht_t *h, long i);
int64_t pg_mht_filter(pg_mht_t *h, int n_proc, int n_tot, double min_freq, int ff);
void pg_mht_tighten(pg_mht_t *h);
void pg_mht_rearrange(pg_mht_t *h, long i);
// void pg_mht_dump(const pg_mht_t *h, const char *fn);
pg_id_map_t *pg_mht_idx(pg_mht_t *h);

pg_csr_t *pg_csr_init(int n_snps, int n_fns, pg_mht_t *h, pg_id_map_t *id_maps);
void pg_csr_insert(pg_csr_t *csr, pg_mht_t *h, pg_id_map_t *id_maps, int gnm_id);
void pg_csr_dump(const char *fn, const pg_mht_t *h, const char **fns, const pg_csr_t *csr);

pg_mht_t *pg_count(const char **fns, const int n_fns, const pg_opt_t *opt);
pg_csr_t *pg_findsnp(const char **fns, const int n_fns, int n_snps, const pg_opt_t *opt, pg_mht_t *h);

#endif // HTAB_H