#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

# include "htab.h"
#include "khashl.h" // hash table
#include "utils.h"

// #define pg_h_eq(a, b) ((a)>>COUNTER_BITS == (b)>>COUNTER_BITS) // lower COUNTER_BITS for counts; higher bits for SNPs (compare only them).
// #define pg_h(a) ((a)>>COUNTER_BITS) // hash only the k-mer part (higher bits).
// KHASHL_SET_INIT(, pg_ht_t, pg_ht, uint96_t, pg_h, pg_h_eq)
KHASHL_MAP_INIT(, pg_ht_t, pg_ht, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

// Operations on hash tables and bloom filter.

pg_mht_t *pg_mht_init(int k, int pre)
{
	pg_mht_t *h;
	int i;
	if (pre < COUNTER_BITS) return 0; // pre should be preceeded by the count bits.
	CALLOC(h, 1);
	h->k = k, h->pre = pre;
	CALLOC(h->h, 1<<h->pre); // allocate the array of partitions.
	for (i = 0; i < 1<<h->pre; ++i) {
		h->h[i].h = pg_ht_init(); // initialize hash table for each bucket.

		pthread_mutex_init(&h->h[i].lock, NULL); // initialize lock for each partition
	}
	return h;
}

void pg_mht_destroy(pg_mht_t *h)
{
	int i;
	if (h == 0) return;
	for (i = 0; i < 1<<h->pre; ++i) {
		pthread_mutex_destroy(&h->h[i].lock);
		pg_ht_destroy(h->h[i].h); // destroy hash table for each bucket.
	}
	free(h->h); free(h);
}

// void pg_mht_filter(pg_mht_t *h, int min_count)
// {
//     int i, n = 1 << h->pre;
//     for (i = 0; i < n; ++i) {
//         pg_ht1_t *g = h->h[i].h;
//         khint_t k;
//         for (k = 0; k < kh_end(g1); ++k)
//             if (kh_exist(g1, k) && (kh_key(g1, k) & COUNTER_MAX) < min_count)
//                 pg_ht1_del(g1, k);
//         for (k = 0; k < kh_end(g2); ++k)
//             if (kh_exist(g2, k) && (kh_key(g2, k) & COUNTER_MAX) < min_count)
//                 pg_ht2_del(g2, k);
//     }
// }

int pg_mht_insert_list(pg_mht_t *h, int n, const ch_seq_t *a)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	pg_ht1_t *g;
	if (n == 0) return 0;
	g = &h->h[a[0].h_flanks&mask]; // get hash table partition for the first (and all) k-mers.
	pthread_mutex_lock(&g->lock);
	
	for (j = 0; j < n; ++j) {
		int absent;
		uint64_t key;
		uint32_t cb = a[j].cb;
		key = a[j].h_flanks >> h->pre;
		if ((a[j].h_flanks & mask) != (a[0].h_flanks & mask)) continue;

		khint_t k = pg_ht_put(g->h, key, &absent);
		if (absent) {
			++n_ins;
			kh_val(g->h, k) = val_pack(1, 0, cb);  // first occurrence, SNP unknown
		} else {
			uint32_t v = kh_val(g->h, k);
			uint32_t cnt = val_count(v);
			uint32_t snp = val_snp(v) | (cb != val_cb(v));  // set SNP if CB differs
			kh_val(g->h, k) = val_pack(cnt < COUNTER_MAX ? cnt + 1 : cnt, snp, val_cb(v));
		}
	}

	pthread_mutex_unlock(&g->lock);
	return n_ins;
}


void pg_mht_tighten(pg_mht_t *h)
{
	int i;
	for (i = 0; i < 1<<h->pre; ++i) {
		pg_ht_t *g = h->h[i].h;
		if (kh_size(g) * 3 < kh_capacity(g))
			pg_ht_m_resize(g, kh_size(g) * 3);
	}
}


int pg_mht_dump(const pg_mht_t *h, const char *fn)
{
	FILE *fp;
	uint32_t t[3];
	int i;
	if ((fp = strcmp(fn, "-")? fopen(fn, "wb") : stdout) == 0) return -1;
	fwrite(PG_MAGIC, 1, 4, fp);
	t[0] = h->k, t[1] = h->pre, t[2] = COUNTER_BITS;
	fwrite(t, 4, 3, fp);
	for (i = 0; i < 1<<h->pre; ++i) {
		pg_ht_t *g = h->h[i].h;
		khint_t k;
		t[0] = kh_capacity(g), t[1] = kh_size(g);
		fwrite(t, 4, 2, fp);
		for (k = 0; k < kh_end(g); ++k)
			if (kh_exist(g, k))
				fwrite(&kh_key(g, k), 8, 1, fp);
	}
	fprintf(stderr, "[M::%s] Dumpped the hash table to file '%s'.\n", __func__, fn);
	fclose(fp);
	return 0;
}