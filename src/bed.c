#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bed.h"
#include "khashl.h"
#include "utils.h"

const bed_ctg_t *bed_get(const bedmap1_t *b, const char *name)
{
	khint_t k;
	if (b == 0) return 0;
	k = bedmap_get(b->h, name);
	return k == kh_end(b->h) ? 0 : &kh_val(b->h, k);
}

// int64_t bed_covered_len(const bed_ctg_t *c, int32_t l)
// {
// 	int64_t cov = 0; int32_t i;
// 	for (i = 0; i < c->n; ++i) {
// 		int32_t st = c->e[i].start, en = c->e[i].end;
// 		if (st < 0) st = 0;
// 		if (en > l) en = l; // BED may overrun the FASTA
// 		if (en > st) cov += en - st;
// 	}
// 	return cov;
// }

int64_t bed_nk(const bed_ctg_t *c, int32_t l, int32_t k)
{
	int64_t nk = 0; int32_t i;
	for (i = 0; i < c->n; ++i) {
		int32_t st = c->e[i].start, en = c->e[i].end;
        // if (st < 0) st = 0; // guard against extended flanks at the start of the sequence
		// if (en > l) en = l; // guard against extended flanks at the end of the sequence
		if (en - st >= k) nk += en - st - k + 1; // clamp at 0: short intervals contribute nothing
	}
	return nk;
}

void mask_fa(char *seq, int32_t l, const bed_ctg_t *c)
{
	int32_t i, prev = 0;
	for (i = 0; i < c->n; ++i) {
		int32_t st = c->e[i].start, en = c->e[i].end;
		// if (st < 0) st = 0; // guard against extended flanks at the start of the sequence
		// if (en > l) en = l; // guard against extended flanks at the end of the sequence
		if (st > prev) memset(seq + prev, 'N', st - prev);
		if (en > prev) prev = en;
	}
	if (prev < l) memset(seq + prev, 'N', l - prev);
}

static int entries_cmp(const void *a, const void *b)
{
	const bed_entry_t *x = (const bed_entry_t*)a, *y = (const bed_entry_t*)b;
	return x->start != y->start ? (x->start < y->start ? -1 : 1)
	     : x->end   != y->end   ? (x->end   < y->end   ? -1 : 1) : 0;
}

static void ctg_sort_merge(bed_ctg_t *c)
{
	int32_t i, j;
	if (c->n < 1) return;
	qsort(c->e, c->n, sizeof(bed_entry_t), entries_cmp);
	for (i = 1, j = 0; i < c->n; ++i) {
		if (c->e[i].start <= c->e[j].end) {
			if (c->e[i].end > c->e[j].end) c->e[j].end = c->e[i].end;
		} else c->e[++j] = c->e[i];
	}
	c->n = j + 1;
}

bedmap1_t *bed_read(const char *fn)
{
	FILE *fp;
	char *line = 0;
	size_t sz = 0;
	bedmap1_t *b;
	int32_t i;
    khint_t kit;

	fp = fopen(fn, "r");
	if (fp == 0) return 0;

	CALLOC(b, 1);
	b->h = bedmap_init();

	while (getline(&line, &sz, fp) > 0) {
		char *name, *p, *q, *save = 0;
		int32_t start, end, idx;
		int absent;
		bed_ctg_t *c;

        // skip lines that are not records
		if (line[0] == '#' || line[0] == '\n') continue;
		if (strncmp(line, "track", 5) == 0 || strncmp(line, "browser", 7) == 0) continue;

        // parse the first three columns of the BED record
		name = strtok_r(line, "\t \n", &save);
		p = strtok_r(0, "\t \n", &save);
		q = strtok_r(0, "\t \n", &save);
		if (name == 0 || p == 0 || q == 0) {
            fprintf(stderr, "[E::%s] malformed BED record: %s\n", __func__, line);
            continue;
        };
		start = (int32_t)strtol(p, 0, 10);
		end = (int32_t)strtol(q, 0, 10);
		if (end <= start) {
            fprintf(stderr, "[E::%s] malformed BED record: %s\n", __func__, line);
            continue;
        };

		kit = bedmap_put(b->h, name, &absent);
        c = &kh_val(b->h, kit); 
		if (absent) {
            kh_key(b->h, kit) = strdup(name);
            memset(c, 0, sizeof(bed_ctg_t));
        }

		if (c->n == c->m) {
            c->m = c->m < 8? 8 : c->m << 1;
            REALLOC(c->e, c->m);
        }
        c->e[c->n].start = start;
        c->e[c->n++].end = end;
	}
	free(line);
	fclose(fp);

	for (kit = 0; kit < kh_end(b->h); ++kit)
		if (kh_exist(b->h, kit))
			ctg_sort_merge(&kh_val(b->h, kit));

	return b;
}

void bed_destroy(bedmap1_t *b)
{
	khint_t kit;
	if (b == 0) return;
	for (kit = 0; kit < kh_end(b->h); ++kit) {
		if (!kh_exist(b->h, kit)) continue;
		free((char*)kh_key(b->h, kit));
		free(kh_val(b->h, kit).e);
	}
	bedmap_destroy(b->h);
	free(b);
}