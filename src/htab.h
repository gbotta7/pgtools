#ifndef HTAB_H
#define HTAB_H

#include <pthread.h>
#include <stdint.h>
#include <stdbool.h>

#include "parser.h"

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
    int idx;    // index of the contig
    int pos;
	uint8_t cb; // stores central base of k-mer and SNP information
    char strand;
} s_ch_seq_t;

typedef struct __attribute__((packed)) {
	uint64_t h_flanks;
	uint8_t cb; // stores central base of k-mer and SNP information
} k_ch_seq_t;

typedef struct { // terminal options
    int64_t chunk_size;
    int64_t n_threads;
    double min_freq;
	int32_t k;
    int32_t pre; // number of bits for partitioning.
    int filt_type;
    int verbose;
} pg_opt_t;

typedef struct {
    int pos;
    int seq_idx;
    char strand;
} kinfo_entry_t;

typedef struct {
    kinfo_entry_t *i;
    int n, m;
} kinfo_t;

KHASHL_MAP_INIT(KH_LOCAL, pg_ht_t, pg_ht, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)
KHASHL_MAP_INIT(KH_LOCAL, pg_im_t, pg_im, uint64_t, kinfo_t, kh_hash_uint64, kh_eq_generic)

// see khashl.h for the definition of pg_ht_t.
typedef struct { 
    struct pg_ht_t *h; // count hash table for each bucket.
} pg_ht1_t;

// see khashl.h and htab.c for the definition of pg_im_t.
typedef struct { 
    struct pg_im_t *m; // info hash table for each bucket.
} pg_im1_t;

typedef struct {
    char **names;   // interned contig name strings
    int n;
    int m;
} cnames_t;

typedef struct {
    cnames_t cnames;        // contig names
    pg_ht1_t *h;            // array of partitions (size = 1 << pre)
    pg_im1_t *m;            // k-mers info
    int64_t n_ins_tot;      // stores the total k-mer insertions in pangenome
    int64_t n_del_tot;      // stores the total filtered k-mers in the filter stage
    int n_done;             // number of processed files
    int32_t k;              // k-mers length
    int32_t pre;            // stores the k-mer flanking sequences, encoded as 2 bits per base (used for partitioning)
} pg_mht_t;


pg_mht_t *pg_mht_init(int k, int pre);
pg_mht_t *pg_mht_copy(const pg_mht_t *src);
void pg_mht_destroy(pg_mht_t *h);
int64_t pg_mht_insert_list(pg_mht_t *h, int n, const k_ch_seq_t *a, int f);
void pg_mht_count_list(pg_mht_t *h, int n, const s_ch_seq_t *a);
void pg_mht_clear_k(pg_mht_t *h, long i, int f);
void pg_mht_clear_s(pg_mht_t *h, long i);
int64_t pg_mht_filter(pg_mht_t *h, int n_proc, int n_tot, double min_freq, int ff);
void pg_mht_tighten(pg_mht_t *h);
void pg_dump_snps(const char *fn, pg_mht_t *h);

void pg_mht_rearrange(pg_mht_t *h, long i);
void write_vcf(const char *out_fn, pg_mht_t *h, pg_mht_t *ref_h, char *gnm_fn);
void merge_vcfs(const char *out_fn, const char *tmpdir, int n_fns, int n_snps);

pg_mht_t *pg_count(const char **fns, const int n_fns, const pg_opt_t *opt, const char *out);
void pg_findsnp(const char **fns, const int n_fns, int64_t n_snps, const pg_opt_t *opt, pg_mht_t *h, const char *ref_fn, const char *out_fn);

#endif // HTAB_H