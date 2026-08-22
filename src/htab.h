#ifndef HTAB_H
#define HTAB_H

#include <pthread.h>
#include <stdint.h>
#include <stdbool.h>

#include "parser.h"

#define PG_MAGIC "PGKM"

// define the first pass kmer hash table value
#define F_VAL_INFO_BITS 8
#define F_PGNM_COUNTER_BITS 20
#define F_GNM_COUNTER_BITS (32 - (F_VAL_INFO_BITS + F_PGNM_COUNTER_BITS))
#define F_COUNTER_MAX ((1U << F_GNM_COUNTER_BITS/2) - 1)

#define f_val_gnm_count2(v) ((v) >> ((F_VAL_INFO_BITS + F_PGNM_COUNTER_BITS) + F_GNM_COUNTER_BITS/2))
#define f_val_gnm_count1(v) (((v) >> (F_VAL_INFO_BITS + F_PGNM_COUNTER_BITS)) & ((1U << F_GNM_COUNTER_BITS/2) - 1))
#define f_val_filt2(v) (((v) >> (7 + F_PGNM_COUNTER_BITS)) & 0x1U)
#define f_val_filt1(v) (((v) >> (6 + F_PGNM_COUNTER_BITS)) & 0x1U)
#define f_val_snp2(v) (((v) >> (5 + F_PGNM_COUNTER_BITS)) & 0x1U)
#define f_val_snp1(v) (((v) >> (4 + F_PGNM_COUNTER_BITS)) & 0x1U)
#define f_val_cb2(v) (((v) >> (2 + F_PGNM_COUNTER_BITS)) & 0x3U)
#define f_val_cb1(v) (((v) >> F_PGNM_COUNTER_BITS) & 0x3U)
#define f_val_pgnm_count2(v) (((v) >> (F_PGNM_COUNTER_BITS/2)) & ((1U << F_PGNM_COUNTER_BITS/2) - 1)) 
#define f_val_pgnm_count1(v) ((v) & ((1U << F_PGNM_COUNTER_BITS/2) - 1))

#define f_val_pack(cnt2, cnt1, filt2, filt1, snp2, snp1, cb2, cb1, pgnm_cnt2, pgnm_cnt1) (((cnt2) << (F_VAL_INFO_BITS + F_PGNM_COUNTER_BITS + F_GNM_COUNTER_BITS/2)) | ((cnt1) << (F_VAL_INFO_BITS + F_PGNM_COUNTER_BITS)) | ((filt2) << (F_PGNM_COUNTER_BITS + 7)) | ((filt1) << (F_PGNM_COUNTER_BITS + 6)) | ((snp2) << (F_PGNM_COUNTER_BITS + 5)) | ((snp1) << (F_PGNM_COUNTER_BITS + 4)) | ((cb2) << (F_PGNM_COUNTER_BITS + 2)) | ((cb1) << F_PGNM_COUNTER_BITS) | ((pgnm_cnt2) << (F_PGNM_COUNTER_BITS/2)) | (pgnm_cnt1))

// define the second pass kmer hash table value
#define S_VAL_INFO_BITS 6
#define S_COUNTER_BITS (32 - (S_VAL_INFO_BITS))
#define S_COUNTER_MAX ((1U << S_COUNTER_BITS/2) - 1)

#define s_val_count2(v) ((v) >> (S_VAL_INFO_BITS + S_COUNTER_BITS/2))
#define s_val_count1(v) (((v) >> S_VAL_INFO_BITS) & ((1U << S_COUNTER_BITS/2) - 1))
#define s_val_filt2(v) (((v) >> 5) & 0x1U)
#define s_val_filt1(v) (((v) >> 4) & 0x1U)
#define s_val_cb2(v) (((v) >> 2) & 0x3U)
#define s_val_cb1(v) ((v) & 0x3U)

#define s_val_pack(cnt2, cnt1, filt2, filt1, cb2, cb1) (((cnt2) << (S_VAL_INFO_BITS + S_COUNTER_BITS/2))| ((cnt1) << S_VAL_INFO_BITS) | ((filt2) << 5) | ((filt1) << 4) | ((cb2) << 2) | (cb1));

// define the second pass kmer info hash table value
#define M_POS_BITS 31
#define M_POS_MAX ((1U << M_POS_BITS) - 1)

#define m_val_pos(v) (((v) >> 1) & 0x7FFFFFFFU)
#define m_val_strand(v) ((v) & 0x1U)

#define m_postrand_pack(pos, strand) (((pos) << 1) | (strand))

typedef struct {
    uint32_t pos;
    uint16_t idx;       // index of the contig
    uint8_t strand;
} ch_info_t;

typedef struct __attribute__((packed)) {
	uint64_t h_flanks;
	uint8_t cb; // stores central base of k-mer and SNP information
} ch_seq_t;

typedef struct { // terminal options
    int64_t chunk_size;
    int64_t n_threads;
    double msf;
    double maf;
    int snp;
    int mko;
	int32_t k;
    int32_t pre; // number of bits for partitioning.
    int write_info;
    int write_mko;
    int filt_type;
    int verbose;
} pg_opt_t;

typedef struct {
    uint32_t postrand; // 1 bit for strand, and 31 for pos
    uint16_t seq_idx;
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
    char **names;           // interned contig name strings
    int n;
    int m;
} cnames_t;

typedef struct {
    pthread_mutex_t mutex;
    cnames_t cnames;        // contig names
    pg_ht1_t *h;            // array of partitions (size = 1 << pre)
    pg_im1_t *m;            // k-mers info
    int64_t n_ins_tot;      // stores the total k-mer insertions in pangenome
    int64_t n_del_tot;      // stores the total filtered k-mers in the filter stage
    int32_t k;              // k-mers length
    int32_t pre;            // stores the k-mer flanking sequences, encoded as 2 bits per base (used for partitioning)
} pg_mht_t;


pg_mht_t *pg_mht_init(int k, int pre);
pg_mht_t *pg_mht_copy(const pg_mht_t *src);
void pg_mht_destroy(pg_mht_t *h);
int64_t pg_mht_insert_list(pg_mht_t *h, int n, const ch_seq_t *a, int f);
void pg_mht_count_list(pg_mht_t *h, int n, const ch_seq_t *a, ch_info_t *b);
void pg_mht_clear1(pg_mht_t *h, long i, int f, int max_occ);
void pg_mht_clear2(pg_mht_t *h, long i);
int64_t pg_mht_filter(pg_mht_t *h, long i, int n_proc, int n_tot, int ff, pg_opt_t *opt);
void pg_mht_tighten(pg_mht_t *h);
pg_mht_t *pg_mht_repopulate(const char *kmer_file, pg_opt_t *opt);
void pg_mht_rearrange(pg_mht_t *h, long i);

pg_mht_t *pg_detect(const char **fa_fns, const char **bed_fns, const int n_fns, const pg_opt_t *opt, const char *out_fn);
void pg_count(const char **fa_fns, const char **bed_fns, const int n_fns, int64_t n_snps, const pg_opt_t *opt, pg_mht_t *h, const char *out_fn);

void pg_dump_kmers(const char *fn, pg_mht_t *h, int snp);
void write_tsv(const char *out_fn, pg_mht_t *h, const char *gnm_fn, int snp);
void merge_tsvs(const char *out_fn, const char *tmpdir, const char **fa_fns, int n_fns, int n_rows);

#endif // HTAB_H