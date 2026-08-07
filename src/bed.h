#ifndef BED_H
#define BED_H

#include <stdint.h>

#include "khashl.h"

typedef struct {
    int32_t start, end;
} bed_entry_t;   // 0-based, half-open [start,end)

typedef struct {
	int32_t n, m;
	bed_entry_t *e;
} bed_ctg_t;

KHASHL_MAP_INIT(KH_LOCAL, bedmap_t, bedmap, const char *, bed_ctg_t, kh_hash_str, kh_eq_str)

typedef struct {
	bedmap_t *h;
} bedmap1_t;


bedmap1_t *bed_read(const char *fn);
void   bed_destroy(bedmap1_t *b);

const bed_ctg_t *bed_get(const bedmap1_t *b, const char *name);
// int64_t bed_covered_len(const bed_ctg_t *c, int32_t l);
int64_t bed_nk(const bed_ctg_t *c, int32_t l, int32_t k);
void mask_fa(char *seq, int32_t l, const bed_ctg_t *c);

#endif