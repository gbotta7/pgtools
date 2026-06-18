#ifndef BFC_BSEQ_H
#define BFC_BSEQ_H

#include <stdint.h>

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

typedef struct {
    char chrom_name[MAX_NAME_LEN];
    long target_pos;
} paf_rec_t;

paf_rec_t *parse_paf(const char *fn, int n_snps, int32_t k);

#endif