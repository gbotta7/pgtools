#ifndef HTAB_H
#define HTAB_H

#include <pthread.h>
#include <stdint.h>
#include <stdbool.h>

#include "parser.h"

// define the first pass SNP-mer hash table value
#define F_PGNM_COUNTER_BITS 26
#define F_PGNM_COUNTER_MAX ((1U << F_PGNM_COUNTER_BITS/2) - 1)
#define F_GNM_COUNTER_BITS (32 - F_PGNM_COUNTER_BITS)
#define F_GNM_COUNTER_MAX ((1U << F_GNM_COUNTER_BITS/2) - 1)

#define f_val_gnm_count2(v) ((v) >> (F_PGNM_COUNTER_BITS + F_GNM_COUNTER_BITS/2))
#define f_val_gnm_count1(v) (((v) >> F_PGNM_COUNTER_BITS) & ((1U << F_GNM_COUNTER_BITS/2) - 1))
#define f_val_pgnm_count2(v) (((v) >> (F_PGNM_COUNTER_BITS/2)) & ((1U << F_PGNM_COUNTER_BITS/2) - 1)) 
#define f_val_pgnm_count1(v) ((v) & ((1U << F_PGNM_COUNTER_BITS/2) - 1))

#define F_VAL_INFO_BITS 7
#define F_VAL_MAX (((uint64_t)1 << F_VAL_INFO_BITS) - 1)
#define f_val_filt(v) (((v) >> 6) & 0x1U)
#define f_val_snp2(v) (((v) >> 5) & 0x1U)
#define f_val_snp1(v) (((v) >> 4) & 0x1U)
#define f_val_cb2(v) (((v) >> 2) & 0x3U)
#define f_val_cb1(v) ((v) & 0x3U)

#define f_key_pack(filt, snp2, snp1, cb2, cb1) (((filt) << 6) | ((snp2) << 5) | ((snp1) << 4) | ((cb2) << 2) | (cb1))
#define f_val_pack(cnt2, cnt1, pgnm_cnt2, pgnm_cnt1) (((cnt2) << (F_PGNM_COUNTER_BITS + F_GNM_COUNTER_BITS/2)) | ((cnt1) << (F_PGNM_COUNTER_BITS)) | ((pgnm_cnt2) << (F_PGNM_COUNTER_BITS/2)) | (pgnm_cnt1))

// define the second pass SNP-mer hash table value
#define S_VAL_INFO_BITS 4
#define S_VAL_MAX (((uint64_t)1 << S_VAL_INFO_BITS) - 1)
#define S_COUNTER_BITS 32
#define S_COUNTER_MAX ((1U << S_COUNTER_BITS/2) - 1)

#define s_val_count2(v) ((v) >> S_COUNTER_BITS/2)
#define s_val_count1(v) ((v) & ((1U << S_COUNTER_BITS/2) - 1))
#define s_val_cb2(v) (((v) >> 2) & 0x3U)
#define s_val_cb1(v) ((v) & 0x3U)

#define s_key_pack(cb2, cb1) (((cb2) << 2) | (cb1))
#define s_val_pack(cnt2, cnt1) (((cnt2) << S_COUNTER_BITS/2) | (cnt1))

// define the second pass SNP-mer info hash table value
#define I_POS_BITS 30
#define I_POS_MAX ((1U << I_POS_BITS) - 1)

#define i_val_pos(v) (((v) >> 2) & I_POS_MAX)
#define i_val_allele(v) ((v) & 0x3U)

#define i_posallele_pack(pos, allele) (((pos) << 2) | (allele))

typedef struct {
    uint32_t pos;
    uint16_t idx;       // index of the contig
    uint8_t allele;
} seq_info_t;

typedef struct __attribute__((packed)) {
	uint64_t h_flanks;
	uint8_t cb; // stores central base of k-mer and SNP information
} seq_t;

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
    uint32_t posallele; // 2 bit for allele, and 30 for pos
    uint16_t seq_idx;
} si_entry_t;

typedef struct {
    si_entry_t *i;
    int n, m;
    uint32_t cnt;
} sci_t; // snp-mer count and info

#define pg_f_eq(a, b) ((a) >> F_VAL_INFO_BITS == (b) >> F_VAL_INFO_BITS)
#define pg_f_hash(a) ((a) >> F_VAL_INFO_BITS)

KHASHL_MAP_INIT(KH_LOCAL, pg_sht_t, pg_sht, uint64_t, uint32_t, pg_f_hash, pg_f_eq)
KHASHL_MAP_INIT(KH_LOCAL, pg_siht_t, pg_siht, uint64_t, sci_t, pg_f_hash, pg_f_eq)

// see khashl.h for the definition of pg_sht_t.
typedef struct { 
    struct pg_sht_t *h; // count hash table for each bucket.
} pg_sht1_t;

// see khashl.h and htab.c for the definition of pg_siht_t.
typedef struct { 
    struct pg_siht_t *ih; // info hash table for each bucket.
} pg_siht1_t;

typedef struct {
    char **names;           // interned contig name strings
    int n;
    int m;
} cnames_t;

typedef struct {
    pthread_mutex_t mutex;
    cnames_t cnames;        // contig names
    pg_sht1_t *h;           // array of partitions (size = 1 << pre)
    pg_siht1_t *ih;          // pg_sht1_t with info
    int64_t n_ins_tot;      // stores the total SNP-mer insertions in pangenome
    int64_t n_del_tot;      // stores the total filtered SNP-mers in the filter stage
    int32_t k;              // SNP-mers length
    int32_t pre;            // partition bits
} pg_msht_t;


pg_msht_t *pg_msht_init(int k, int pre, int w);
void pg_msht_destroy(pg_msht_t *h, int w);
int64_t pg_msht_insert_list(pg_msht_t *h, int n, const seq_t *a, int f);
void pg_msht_count_list(pg_msht_t *h, int n, const seq_t *a, seq_info_t *b);
void pg_msht_clear1(pg_msht_t *h, long i, int f, int max_occ);
void pg_msht_clear2(pg_msht_t *h, long i, int w);
int64_t pg_msht_filter(pg_msht_t *h, long i, int n_proc, int n_tot, int ff, pg_opt_t *opt);
void pg_msht_tighten(pg_msht_t *h);
pg_msht_t *pg_msht_repopulate(const char *kmer_file, pg_opt_t *opt);
void pg_msht_rearrange(pg_msht_t *h, long i);

pg_msht_t *pg_detect(const char **fa_fns, const char **bed_fns, const int n_fns, const pg_opt_t *opt, const char *out_fn);
void pg_count(const char *fa_fn, const char *bed_fn, const pg_opt_t *opt, pg_msht_t *h, const char *out_fn);

void pg_dump_snpmers(const char *fn, pg_msht_t *h);
void write_snpmer_tsv(const char *out_fn, pg_msht_t *h, const char *gnm_fn, int w, int w_mko);
// void merge_tsvs(const char *out_fn, const char *tmpdir, const char **fa_fns, int n_fns, int n_rows);

#endif // HTAB_H