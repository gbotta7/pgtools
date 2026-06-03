#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

# include "htab.h"
#include "khashl.h" // hash table
#include "utils.h"

KHASHL_MAP_INIT(, pg_ht_t, pg_ht, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

// Operations on hash tables and bloom filter.

pg_mht_t *pg_mht_init(int k, int pre)
{
	pg_mht_t *h;
	int i;
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

int pg_mht_filter(pg_mht_t *h, int n_proc, int n_tot, double min_freq)
{	
	int n_del = 0;
    int i, n = 1 << h->pre;

    for (i = 0; i < n; ++i) {
		// store entries to delete
        pg_ht1_t *g = &h->h[i];
		uint64_t *del_part = malloc(kh_size(g->h) * sizeof(uint64_t));
		int n_del_part = 0;
        khint_t k;
        // pthread_mutex_lock(&g->lock);
        for (k = 0; k < kh_end(g->h); ++k) {
			if (!kh_exist(g->h, k)) continue;
			uint32_t v = kh_val(g->h, k);
			if ((double)(n_proc - (val_count1(v)) + (val_count2(v))) / n_tot > (1.0 - min_freq) || (!val_snp1(v) || val_snp2(v))) {
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
        // pthread_mutex_unlock(&g->lock);
    }

	return n_del;
}

int pg_mht_insert_list(pg_mht_t *h, int n, const ch_seq_t *a, uint64_t gnm_id)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	pg_ht1_t *g;
	if (n == 0) return 0;

	g = &h->h[a[0].h_flanks & mask]; // get hash table partition for the first (and all) k-mers.
	pthread_mutex_lock(&g->lock);
	
	for (j = 0; j < n; ++j) {
		int absent;
		uint32_t cb = a[j].cb;
		uint64_t key;
		key = (a[j].h_flanks >> h->pre) | (gnm_id << (8*sizeof(uint64_t) - h->pre));

		uint32_t cnt1 = 0, cnt2 = 0, cb1, cb2, snp1, snp2, v;

		khint_t k = pg_ht_get(g->h, a[j].h_flanks >> h->pre);

		// SNP-mers pass
		if (gnm_id && k == kh_end(g->h)) 
			continue; // this k-mer does not exist in the pangenome hash table (with gnm_id=0)
		else if (gnm_id && k != kh_end(g->h)) { // this k-mer exists in the pangenome hash table, check if it exists in this genome and get the counts
			uint32_t v_pgnm = kh_val(g->h, k);
			// Inherit SNP state and cb assignments from pangenome
			snp1 = val_snp1(v_pgnm);
			snp2 = val_snp2(v_pgnm);
			cb1  = val_cb1(v_pgnm);
			cb2  = val_cb2(v_pgnm);

			k = pg_ht_put(g->h, key, &absent);

			if (absent) {
				++n_ins;
				if (cb == cb1) {
					cnt1 = 1;
				}
				else {
					cnt2 = 1;
				}
				kh_val(g->h, k) = val_pack(cnt2, cnt1, snp2, snp1, cb2, cb1); // first occurrence in this genome, SNP unknown
			}
			else {
				v = kh_val(g->h, k);
				cnt1 = val_count1(v);
				cnt2 = val_count2(v);

				if (cb == cb1) {
					if (cnt1 < COUNTER_MAX) ++cnt1;
				} else {
					if (cnt2 < COUNTER_MAX) ++cnt2;
				}
			
				kh_val(g->h, k) = val_pack(cnt2, cnt1, snp2, snp1, cb2, cb1);
			}
		}

		// k-mers pass
		else {
			k = pg_ht_put(g->h, key, &absent);
			if (absent) {
				++n_ins;
				kh_val(g->h, k) = val_pack(0, 1, 0, 0, 0, cb);  // first occurrence, SNP unknown
			} else {
				v = kh_val(g->h, k);
				cb1 = val_cb1(v);
				cb2 = val_cb2(v);
				cnt1 = val_count1(v);
				cnt2 = val_count2(v);
				if (val_snp1(v) ^ val_snp2(v)) { // already known as SNP, check if it is multi-allelic
					snp1 = 1;
					if (cb != val_cb1(v) && cb != val_cb2(v)) {
						snp2 = 1; // multi-allelic SNP (do not count)
					} else {
						snp2 = 0; // bi-allelic SNP
						if (cb == val_cb1(v))
							if (cnt1 < COUNTER_MAX) ++cnt1;
						else 
							if (cnt2 < COUNTER_MAX) ++cnt2;
					}
				} else {
					if (val_snp1(v) & val_snp2(v)) { // already known as multi-allelic SNP
						snp1 = 1; snp2 = 1; // already known as multi-allelic SNP (do not count)
						if (cb == val_cb1(v))
							if (cnt1 < COUNTER_MAX) cnt1 += 1;
						else if (cb == val_cb2(v)) 
							if (cnt2 < COUNTER_MAX) cnt2 += 1;
						else continue;
					} else if (cb != val_cb1(v)) { // newly identified SNP
						snp1 = 1; snp2 = 0;
						if (cnt2 < COUNTER_MAX) ++cnt2;
						cb2 = cb; // store the second central base
					} else { // still non-SNP
						snp1 = 0; snp2 = 0;
						if (cnt1 < COUNTER_MAX) ++cnt1;
					}
				}

				kh_val(g->h, k) = val_pack(cnt2, cnt1, snp2, snp1, cb2, cb1);
			}
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


pg_mtx_t *pg_mht2mtx(const pg_mht_t *h, int n_fns, int n_kmers)
{
    pg_mtx_t *m;
    CALLOC(m, 1);

    m->n_rows = 2 * n_kmers;
    m->n_cols = n_fns;

    m->snps = (pg_snp_t *)calloc(n_kmers, sizeof(pg_snp_t));
    m->mat  = (uint64_t *)calloc(m->n_rows * m->n_cols, sizeof(uint64_t));

    pg_mht_t *h_copy = pg_mht_init(h->k, h->pre);

    int idx = 0;

    for (int i = 0; i < (1 << h->pre); ++i) {
        pg_ht_t *g = h->h[i].h;
        pg_ht_t *g_copy = h_copy->h[i].h;

        for (khint_t k = 0; k < kh_end(g); ++k) {
            if (!kh_exist(g, k))
                continue;

            uint64_t key = kh_key(g, k) & ((1ULL << (64 - h->pre)) - 1);
			khint_t k_pgnm = pg_ht_get(g, key);  // key has gnm_id bits = 0
			if (k_pgnm == kh_end(g)) continue;
			uint32_t v_pgnm = kh_val(g, k_pgnm);

            int absent;
            khint_t k_copy = pg_ht_put(g_copy, key, &absent);

            if (!absent)
                continue;

            m->snps[idx].flanks = key;
            m->snps[idx].cb1 = val_cb1(v_pgnm);
            m->snps[idx].cb2 = val_cb2(v_pgnm);

            for (int gnm_id = 1; gnm_id <= n_fns; ++gnm_id) {
                uint64_t key_gnm =
                    ((uint64_t)gnm_id << (64 - h->pre)) | key;

                khint_t k_gnm = pg_ht_get(g, key_gnm);

                if (k_gnm == kh_end(g))
                    continue;

                uint32_t v = kh_val(g, k_gnm);

                m->mat[(2 * idx + 0) * n_fns + (gnm_id - 1)] =
                    val_count1(v);

                m->mat[(2 * idx + 1) * n_fns + (gnm_id - 1)] =
                    val_count2(v);
            }

            ++idx;
        }
    }

    pg_mht_destroy(h_copy);

    return m;
}

void pg_mtx_dump(const char *fn, const pg_mht_t *h, const char **fns, const pg_mtx_t *m)
{
    FILE *fp;

    if ((fp = strcmp(fn, "-") ? fopen(fn, "w") : stdout) == 0)
        return;

    fprintf(fp, "kmer\tSNP");

    for (int j = 0; j < m->n_cols; ++j)
        fprintf(fp, "\t%s", fns[j]);

    fprintf(fp, "\n");

    uint64_t hash_mask =
        (1ULL << ((h->k - 1) * 2)) - 1;

    char seq[h->k + 1];
    int mid = h->k >> 1;

    for (int s = 0; s < m->n_rows/2; ++s) {

        uint64_t flanks =
            pg_hash64_inv(m->snps[s].flanks << h->pre,
                          hash_mask);

        for (int j = 0; j < mid; ++j)
            seq[h->k - 1 - j] =
                nt4_seq_table[(flanks >> (j * 2)) & 3];

        for (int j = 0; j < mid; ++j)
            seq[mid - 1 - j] =
                nt4_seq_table[(flanks >> ((mid + j) * 2)) & 3];

        seq[h->k] = '\0';
		
        seq[mid] = nt4_seq_table[m->snps[s].cb1];

        fprintf(fp, "%s\t%c/%c",
                seq,
                nt4_seq_table[m->snps[s].cb1],
                nt4_seq_table[m->snps[s].cb2]);

        for (int j = 0; j < m->n_cols; ++j)
            fprintf(fp, "\t%llu",
                    (unsigned long long)
                    m->mat[(2 * s + 0) * m->n_cols + j]);

        fprintf(fp, "\n");

        seq[mid] = nt4_seq_table[m->snps[s].cb2];

        fprintf(fp, "%s\t%c/%c",
                seq,
                nt4_seq_table[m->snps[s].cb1],
                nt4_seq_table[m->snps[s].cb2]);

        for (int j = 0; j < m->n_cols; ++j)
            fprintf(fp, "\t%llu",
                    (unsigned long long)
                    m->mat[(2 * s + 1) * m->n_cols + j]);

        fprintf(fp, "\n");
    }

    fclose(fp);
}