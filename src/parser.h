#ifndef BFC_BSEQ_H
#define BFC_BSEQ_H

#include <stdint.h>

#include "khashl.h"

#define MAX_NAME_LEN  32
#define LINE_BUF_SIZE 65536

struct bseq_file_s;
typedef struct bseq_file_s bseq_file_t;

typedef struct {
	int l_seq;
	uint32_t aux, aux2;
	char *name, *comment, *seq, *qual;
} bseq1_t;

extern unsigned char seq_nt6_table[256];

bseq_file_t *bseq_open(const char *fn);
void bseq_close(bseq_file_t *fp);
bseq1_t *bseq_read(bseq_file_t *fp, int64_t chunk_size, int keep_comment, int *n_);

// typedef struct {
//     char chrom_name[MAX_NAME_LEN];
//     long target_pos;
// } paf_rec_t;

// typedef struct {
//     off_t offset;    // byte offset in file
//     char chrom_name[MAX_NAME_LEN];
//     char qname[MAX_NAME_LEN];
// } paf_entry_t;

// struct paf_idx_t; // see khashl.h and htab.c for the definition of pg_ht_t.

// KHASHL_MAP_INIT(KH_LOCAL, paf_idx_t, paf_idx, uint64_t, paf_entry_t, kh_hash_uint64, kh_eq_generic)

// paf_idx_t *idx_paf(const char *fn, int32_t k);
// paf_rec_t fetch_paf_rec(paf_idx_t *m, FILE *fp, uint64_t hflanks, int32_t k);

#endif