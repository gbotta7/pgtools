#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "htab.h"
#include "khashl.h" // hash table
#include "utils.h"

KHASHL_MAP_INIT(, pg_ht_t, pg_ht, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

// Operations on hash tables and bloom filter.

pg_mht_t *pg_mht_init(int k, int pre)
{
	pg_mht_t *h;
	int i;
	CALLOC(h, 1);
	h->k = k;
	h->pre = pre;
	h->n_done = 0;
	CALLOC(h->h, 1<<h->pre); // allocate the array of partitions.
	for (i = 0; i < 1<<h->pre; ++i) {
		h->h[i].h = pg_ht_init(); // initialize hash table for each bucket.
	}
	return h;
}

// pg_ht_t *pg_ht_copy(const pg_ht_t *src) {
//     if (!src) return NULL;
//     pg_ht_t *dst = (pg_ht_t*)kcalloc(1, sizeof(pg_ht_t));
//     dst->bits  = src->bits;
//     dst->count = src->count;
//     khint_t n_buckets = kh_capacity(src);
//     if (n_buckets) {
//         // copy used flags
//         size_t flag_size = __kh_fsize(n_buckets) * sizeof(khint32_t);
//         dst->used = (khint32_t*)kmalloc(flag_size);
//         memcpy(dst->used, src->used, flag_size);
//         // copy keys+vals (bucket array)
//         dst->keys = kmalloc(n_buckets * sizeof(*src->keys));
//         memcpy(dst->keys, src->keys, n_buckets * sizeof(*src->keys));
//     }
//     return dst;
// }

pg_ht_t *pg_ht_copy(const pg_ht_t *src) {
    if (!src) return NULL;
    pg_ht_t *dst = (pg_ht_t*)kcalloc(1, sizeof(pg_ht_t));
    dst->bits  = src->bits;
    dst->count = src->count;
    khint_t n_buckets = kh_capacity(src);
    if (n_buckets) {
        size_t flag_size = __kh_fsize(n_buckets) * sizeof(khint32_t);
        dst->used = (khint32_t*)kmalloc(flag_size);
        memcpy(dst->used, src->used, flag_size);
        size_t bucket_size = n_buckets * sizeof(pg_ht_t_m_bucket_t);
        dst->keys = kmalloc(bucket_size);
        memcpy(dst->keys, src->keys, bucket_size);
    }
    return dst;
}

pg_mht_t *pg_mht_copy(const pg_mht_t *src) {
    pg_mht_t *dst = (pg_mht_t*)calloc(1, sizeof(pg_mht_t));
    dst->k = src->k;
    dst->pre = src->pre;
    dst->n_ins_tot = src->n_ins_tot;
    dst->n_del_tot = src->n_del_tot;
    dst->n_done = src->n_done;
    int n = 1 << src->pre;
    dst->h = (pg_ht1_t*)calloc(n, sizeof(pg_ht1_t));
    for (int i = 0; i < n; ++i)
        dst->h[i].h = pg_ht_copy(src->h[i].h);
    return dst;
}

void pg_mht_destroy(pg_mht_t *h)
{
	int i;
	if (h == 0) return;
	for (i = 0; i < 1<<h->pre; ++i) {
		pg_ht_destroy(h->h[i].h); // destroy hash table for each bucket.
	}
	free(h->h); free(h);
}

int64_t pg_mht_filter(pg_mht_t *h, int n_proc, int n_tot, double min_freq, int ff)
{	
	int64_t n_del = 0;
    int i, n = 1 << h->pre;
	int cond;

    for (i = 0; i < n; ++i) {
		// store entries to delete
        pg_ht1_t *g = &h->h[i];
		uint64_t *del_part = malloc(kh_size(g->h) * sizeof(uint64_t));
		int64_t n_del_part = 0;
        khint_t k;
        for (k = 0; k < kh_end(g->h); ++k) {
			if (!kh_exist(g->h, k)) continue;
			uint32_t v = kh_val(g->h, k);
			if (ff) { // final filter
				cond = (double)(n_proc - (k_val_pgnm_count1(v) + k_val_pgnm_count2(v))) / n_tot > (1.0 - min_freq + 1e-9) || (!k_val_snp1(v) || k_val_snp2(v));
			} else {
				cond = (double)(n_proc - (k_val_pgnm_count1(v) + k_val_pgnm_count2(v))) / n_tot > (1.0 - min_freq + 1e-9);
			}

			if (cond) {
				del_part[n_del_part++] = kh_key(g->h, k);
			}
		}

		// delete entries
		for (int d = 0; d < n_del_part; ++d) {
			k = pg_ht_get(g->h, del_part[d]);
    		if (k != kh_end(g->h)) {
				pg_ht_del(g->h, k);
			}
		}
		n_del += n_del_part;
		free(del_part);
    }

	return n_del;
}

int64_t pg_mht_insert_list(pg_mht_t *h, int n, const ch_seq_t *a, int f)
{
	int j, mask = (1<<h->pre) - 1;
	int64_t n_ins = 0;
	pg_ht1_t *g;
	if (n == 0) return 0;

	uint32_t gnm_cnt1, gnm_cnt2, cb1, cb2, snp1, snp2, pgnm_cnt1, pgnm_cnt2, v;

	g = &h->h[a[0].h_flanks & mask]; // get hash table partition for the first (and all) k-mers.
	
	for (j = 0; j < n; ++j) {
		int absent;
		uint32_t cb = a[j].cb;
		uint64_t key;
		key = (a[j].h_flanks >> h->pre);
		
		// k-mers pass
		khint_t k;
		if (!f) {
			k = pg_ht_put(g->h, key, &absent);
		} else {
			k = pg_ht_get(g->h, key);
			if (k == kh_end(g->h)) continue; // k-mer not found until now, skip
			absent = 0;
		}
		if (absent) { // first occurrence, SNP unknown
			++n_ins;
			gnm_cnt1 = 1; gnm_cnt2 = 0; snp1 = 0; snp2 = 0; cb1 = cb; cb2 = 0; pgnm_cnt1 = 0; pgnm_cnt2 = 0;
			kh_val(g->h, k) = k_val_pack(gnm_cnt2, gnm_cnt1, snp2, snp1, cb2, cb1, pgnm_cnt2, pgnm_cnt1);
		} else {
			v = kh_val(g->h, k);
			gnm_cnt1 = k_val_gnm_count1(v);
			gnm_cnt2 = k_val_gnm_count2(v);
			cb1 = k_val_cb1(v);
			cb2 = k_val_cb2(v);
			snp1 = k_val_snp1(v);
			snp2 = k_val_snp2(v);
			pgnm_cnt1 = k_val_pgnm_count1(v);
			pgnm_cnt2 = k_val_pgnm_count2(v);
			if (snp1 ^ snp2) { // already known as SNP, check if it is multi-allelic
				snp1 = 1;
				if (cb != cb1 && cb != cb2) {
					snp2 = 1; // multi-allelic SNP (do not count)
				} else {
					snp2 = 0; // bi-allelic SNP
					if (cb == cb1) {
						if (gnm_cnt1 < K_COUNTER_MAX) ++gnm_cnt1;
					}
					else  {
						if (gnm_cnt2 < K_COUNTER_MAX) {
							++gnm_cnt2;
						}
					}
				}
			} else if (snp1 & snp2) {
				snp1 = 1; snp2 = 1; // already known as multi-allelic SNP
				if (cb == cb1) {
					if (gnm_cnt1 < K_COUNTER_MAX) {
						gnm_cnt1 += 1;
					}
				}
				else if (cb == cb2) {
					if (gnm_cnt2 < K_COUNTER_MAX) {
						gnm_cnt2 += 1;
					} 
				} else {
						continue;
				}
			} else if (cb != cb1) { // newly identified SNP
					snp1 = 1; snp2 = 0;
					if (gnm_cnt2 < K_COUNTER_MAX) {
						++gnm_cnt2;
					}
					cb2 = cb; // store the second central base
			} else { // still non-SNP
				snp1 = 0; snp2 = 0;
				if (gnm_cnt1 < K_COUNTER_MAX) {
					++gnm_cnt1;
				}
			}

			kh_val(g->h, k) = k_val_pack(gnm_cnt2, gnm_cnt1, snp2, snp1, cb2, cb1, pgnm_cnt2, pgnm_cnt1);
		}
	}
	
	return n_ins;
}


void pg_mht_clear_k(pg_mht_t *h, long i)
{
	// store entries to delete
	pg_ht1_t *g = &h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g->h); ++k) {
		if (!kh_exist(g->h, k)) continue;
		uint32_t v = kh_val(g->h, k);
		if (k_val_gnm_count1(v) == 1 && k_val_gnm_count2(v) == 0) {
			kh_val(g->h, k) = k_val_pack(0, 0, k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v), k_val_pgnm_count1(v)+1);
		} else if (k_val_gnm_count1(v) == 0 && k_val_gnm_count2(v) == 1) {
			kh_val(g->h, k) = k_val_pack(0, 0, k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v)+1, k_val_pgnm_count1(v));
		} else {
			kh_val(g->h, k) = k_val_pack(0, 0, k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v), k_val_pgnm_count1(v));
		}
	}
}

void pg_mht_clear_s(pg_mht_t *h, long i)
{
	// store entries to delete
	pg_ht1_t *g = &h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g->h); ++k) {
		if (!kh_exist(g->h, k)) continue;
		uint32_t v = kh_val(g->h, k);
		kh_val(g->h, k) = s_val_pack(0, 0, s_val_snp2(v), s_val_snp1(v), s_val_cb2(v), s_val_cb1(v));
	}
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


void pg_mht_rearrange(pg_mht_t *h, long i)
{
	// store entries to delete
	pg_ht1_t *g = &h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g->h); ++k) {
		if (!kh_exist(g->h, k)) continue;
		uint32_t v = kh_val(g->h, k);

		uint32_t cb1 = k_val_cb1(v);
		uint32_t cb2 = k_val_cb2(v);
		uint32_t snp1 = k_val_snp1(v);
		uint32_t snp2 = k_val_snp2(v);

		kh_val(g->h, k) = s_val_pack(0, 0, snp2, snp1, cb2, cb1);
	}
}


void pg_mht_count_list(pg_mht_t *h, int n, const ch_seq_t *a)
{
	int j, mask = (1<<h->pre) - 1;
	pg_ht1_t *g;
	if (n == 0) return;

	uint32_t cnt1, cnt2, cb1, cb2, snp1, snp2, v;

	g = &h->h[a[0].h_flanks & mask]; // get hash table partition for the first (and all) k-mers.
	
	for (j = 0; j < n; ++j) {
		int absent;
		uint32_t cb = a[j].cb;
		uint64_t key;
		key = (a[j].h_flanks >> h->pre);
		
		// snp-mers pass
		khint_t k = pg_ht_get(g->h, key);
		if (k == kh_end(g->h)) continue; // not a SNP-mer, skip

		v = kh_val(g->h, k);
		cnt1 = s_val_count1(v);
		cnt2 = s_val_count2(v);
		cb1 = s_val_cb1(v);
		cb2 = s_val_cb2(v);
		snp1 = s_val_snp1(v);
		snp2 = s_val_snp2(v);
		if (snp1 ^ snp2) { // already known as SNP, check if it is multi-allelic
			snp1 = 1;
			if (cb != cb1 && cb != cb2) {
				snp2 = 1; // multi-allelic SNP (do not count)
			} else {
				snp2 = 0; // bi-allelic SNP
				if (cb == cb1) {
					if (cnt1 < S_COUNTER_MAX) ++cnt1;
				}
				else  {
					if (cnt2 < S_COUNTER_MAX) {
						++cnt2;
					}
				}
			}
		} else if (snp1 & snp2) {
			snp1 = 1; snp2 = 1; // already known as multi-allelic SNP
			if (cb == cb1) {
				if (cnt1 < S_COUNTER_MAX) {
					cnt1 += 1;
				}
			}
			else if (cb == cb2) {
				if (cnt2 < S_COUNTER_MAX) {
					cnt2 += 1;
				} 
			} else {
					continue;
			}
		} else if (cb != cb1) { // newly identified SNP
				snp1 = 1; snp2 = 0;
				if (cnt2 < S_COUNTER_MAX) {
					++cnt2;
				}
				cb2 = cb; // store the second central base
		} else { // still non-SNP
			snp1 = 0; snp2 = 0;
			if (cnt1 < S_COUNTER_MAX) {
				++cnt1;
			}
		}

		kh_val(g->h, k) = s_val_pack(cnt2, cnt1, snp2, snp1, cb2, cb1);
	}
}

pg_id_map_t *pg_mht_idx(pg_mht_t *h)
{
	pg_id_map_t *id_maps = calloc(1<<h->pre, sizeof(pg_id_map_t));
	uint32_t n_snps = 0;

	for (int i = 0; i < 1<<h->pre; ++i) {
		pg_ht1_t *g = &h->h[i];
		uint32_t cap = kh_capacity(g->h);
		id_maps[i].n = cap;
		id_maps[i].ids = malloc(cap * sizeof(uint32_t));
		memset(id_maps[i].ids, 0xFF, cap * sizeof(uint32_t)); // UINT32_MAX = not a SNP
		for (khint_t k = 0; k < kh_end(g->h); ++k) {
			if (!kh_exist(g->h, k)) continue;
			uint32_t v = kh_val(g->h, k);
			id_maps[i].ids[k] = n_snps++;
		}
	}

	return id_maps;
}


void pg_csr_insert(pg_csr_t *csr, pg_mht_t *h, pg_id_map_t *id_maps, int gnm_id)
{	
	for (int i = 0; i < 1<<h->pre; ++i) {
		pg_ht1_t *g = &h->h[i];
		for (khint_t k = 0; k < kh_end(g->h); ++k) {
			if (!kh_exist(g->h, k)) continue;
			// if (id_maps[i].ids[k] == UINT32_MAX) continue; // not a SNP-mer

			uint32_t v = kh_val(g->h, k);
			uint32_t cnt1 = s_val_count1(v);
			uint32_t cnt2 = s_val_count2(v);
			if (cnt1 == 0 && cnt2 == 0) continue; // not seen in this genome

			if (csr->n == csr->m) {
				csr->m = csr->m < 1024 ? 1024 : csr->m + (csr->m >> 1);
				REALLOC(csr->entries, csr->m);
			}
			csr->entries[csr->n].row_id = id_maps[i].ids[k];
			csr->entries[csr->n].col_id = gnm_id;
			csr->entries[csr->n].cnt1 = cnt1 > M_COUNTER_MAX ? M_COUNTER_MAX : (uint16_t)cnt1;
			csr->entries[csr->n].cnt2 = cnt2 > M_COUNTER_MAX ? M_COUNTER_MAX : (uint16_t)cnt2;
			csr->n++;
		}
	}
}

pg_csr_t *pg_csr_init(int n_snps, int n_fns, pg_mht_t *h, pg_id_map_t *id_maps)
{
    pg_csr_t *csr;
    CALLOC(csr, 1);
    csr->n_snps = n_snps;
    csr->n_fns  = n_fns;

    // fill snpmer metadata once
    csr->snpmer = calloc(n_snps, sizeof(pg_snp_t));
    for (int i = 0; i < 1<<h->pre; ++i) {
        pg_ht1_t *g = &h->h[i];
        for (khint_t k = 0; k < kh_end(g->h); ++k) {
            if (!kh_exist(g->h, k)) continue;
            uint32_t id = id_maps[i].ids[k];
            // if (id == UINT32_MAX) continue;
            uint32_t v = kh_val(g->h, k);
            csr->snpmer[id].cb1 = (uint8_t)s_val_cb1(v);
            csr->snpmer[id].cb2 = (uint8_t)s_val_cb2(v);
            csr->snpmer[id].flanks = kh_key(g->h, k);
        }
    }
	return csr;
}

static int csr_entry_cmp(const void *a, const void *b) {
    const pg_entry_t *ea = (const pg_entry_t *)a;
    const pg_entry_t *eb = (const pg_entry_t *)b;
    if (ea->row_id != eb->row_id) return (ea->row_id > eb->row_id) - (ea->row_id < eb->row_id);
    return (ea->col_id > eb->col_id) - (ea->col_id < eb->col_id);
}

void pg_csr_dump(const char *fn, const pg_mht_t *h, const char **fns, const pg_csr_t *csr)
{
    FILE *fp;
    if ((fp = strcmp(fn, "-") ? fopen(fn, "w") : stdout) == 0)
        return;

    fprintf(fp, "kmer\tSNP");
    for (int j = 0; j < csr->n_fns; ++j)
        fprintf(fp, "\t%s", fns[j]);
    fprintf(fp, "\n");

    uint64_t hash_mask = (1ULL << ((h->k - 1) * 2)) - 1;
    char seq[h->k + 1];
    int mid = h->k >> 1;

    uint16_t *row1 = calloc(csr->n_fns, sizeof(uint16_t));
    uint16_t *row2 = calloc(csr->n_fns, sizeof(uint16_t));

    qsort(csr->entries, csr->n, sizeof(pg_entry_t), csr_entry_cmp);

    int64_t e = 0;
    for (uint32_t s = 0; s < csr->n_snps; ++s) {

        memset(row1, 0, csr->n_fns * sizeof(uint16_t));
        memset(row2, 0, csr->n_fns * sizeof(uint16_t));

        while (e < csr->n && csr->entries[e].row_id == s) {
            uint16_t c = csr->entries[e].col_id;
            row1[c] = csr->entries[e].cnt1;
            row2[c] = csr->entries[e].cnt2;
            e++;
        }

        uint64_t flanks = pg_hash64_inv(csr->snpmer[s].flanks << h->pre, hash_mask);
        for (int j = 0; j < mid; ++j)
            seq[h->k - 1 - j] = nt4_seq_table[(flanks >> (j * 2)) & 3];
        for (int j = 0; j < mid; ++j)
            seq[mid - 1 - j] = nt4_seq_table[(flanks >> ((mid + j) * 2)) & 3];
        seq[h->k] = '\0';

        seq[mid] = nt4_seq_table[csr->snpmer[s].cb1];
        fprintf(fp, "%s\t%c/%c", seq,
                nt4_seq_table[csr->snpmer[s].cb1],
                nt4_seq_table[csr->snpmer[s].cb2]);
        for (int j = 0; j < csr->n_fns; ++j)
            fprintf(fp, "\t%u", row1[j]);
        fprintf(fp, "\n");

        seq[mid] = nt4_seq_table[csr->snpmer[s].cb2];
        fprintf(fp, "%s\t%c/%c", seq,
                nt4_seq_table[csr->snpmer[s].cb1],
                nt4_seq_table[csr->snpmer[s].cb2]);
        for (int j = 0; j < csr->n_fns; ++j)
            fprintf(fp, "\t%u", row2[j]);
        fprintf(fp, "\n");
    }

    free(row1); free(row2);
    if (fp != stdout) fclose(fp);
}
