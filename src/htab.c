#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include "htab.h"
#include "khashl.h" // hash table
#include "kseq.h"
#include "utils.h"

KHASHL_SET_INIT(, strset_t, strset, const char *, kh_hash_str, kh_eq_str)

// Operations on hash tables

pg_mht_t *pg_mht_init(int k, int pre)
{
	pg_mht_t *h;
	int i;
	CALLOC(h, 1);
	h->k = k;
	h->pre = pre;
	pthread_mutex_init(&h->mutex, 0);
	CALLOC(h->h, 1<<h->pre); // allocate the array of partitions.
	CALLOC(h->m, 1<<h->pre);
	for (i = 0; i < 1<<h->pre; ++i) {
		h->h[i].h = pg_ht_init(); // initialize hash table for each bucket.
		h->m[i].m = pg_im_init();
	}
	return h;
}

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
	pthread_mutex_init(&dst->mutex, 0); // dst has own mutex
    int n = 1 << src->pre;
    dst->h = (pg_ht1_t*)calloc(n, sizeof(pg_ht1_t));
	dst->m = (pg_im1_t*)calloc(n, sizeof(pg_im1_t));
    for (int i = 0; i < n; ++i) {
		dst->h[i].h = pg_ht_copy(src->h[i].h);
		dst->m[i].m = pg_im_init();
	}
    return dst;
}

void pg_im1_destroy(pg_im1_t *m)
{
    khint_t k;
    for (k = 0; k < kh_end(m->m); ++k) {
        if (!kh_exist(m->m, k)) continue;
        kinfo_t *info = &kh_val(m->m, k);
        free(info->i);
    }
    pg_im_destroy(m->m);
}

void pg_mht_destroy(pg_mht_t *h)
{
	int i;
	if (h == 0) return;
	pthread_mutex_destroy(&h->mutex);
	for (i = 0; i < 1<<h->pre; ++i) {
		pg_ht_destroy(h->h[i].h); // destroy hash table for each bucket.
		pg_im1_destroy(&h->m[i]);
	}
	free(h->h); free(h->m); free(h);
}

int64_t pg_mht_filter(pg_mht_t *h, long i, int n_proc, int n_tot, int ff, pg_opt_t *opt)
{	
	int64_t n_del = 0;
	int cond = 0, cond_msf = 0, cond_msf1, cond_msf2, cond_maf = 0, cond_snp = 0;

	pg_ht1_t *g = &h->h[i];
	uint64_t *del_part = malloc(kh_size(g->h) * sizeof(uint64_t));
	int64_t n_del_part = 0;
	khint_t k;
	uint32_t v, filt1, filt2;
	for (k = 0; k < kh_end(g->h); ++k) {
		if (!kh_exist(g->h, k)) continue;
		v = kh_val(g->h, k);
		filt1 = f_val_filt1(v);
		filt2 = f_val_filt2(v);
		if (opt->snp) { // SNP-mer filters
			cond_msf = (double)(n_proc - (f_val_pgnm_count1(v) + f_val_pgnm_count2(v))) / n_tot > (1.0 - opt->msf + 1e-9);		// check that one of the two alleles is present in at least opt->msf*N files
			if (ff) cond_maf = (f_val_pgnm_count1(v) / (double)(f_val_pgnm_count1(v) + f_val_pgnm_count2(v))) < opt->maf || (f_val_pgnm_count2(v) / (double)(f_val_pgnm_count1(v) + f_val_pgnm_count2(v))) < opt->maf;		// check minimum allelic frequency > opt->maf (only in ff)
			if (ff) cond_snp = !f_val_snp1(v) || f_val_snp2(v);										// remove multi-allelic and non SNP-mers
			cond = cond_msf || cond_maf || cond_snp || (f_val_filt1(v) || f_val_filt2(v));			// just one of the filt bits needs to be on in SNP-mers
		} else { // k-mer filters
			cond_msf1 = (double)(n_proc - f_val_pgnm_count1(v)) / n_tot > (1.0 - opt->msf + 1e-9);		// check that each kmer is present in at least opt->msf*N files
			cond_msf2 = (double)(n_proc - f_val_pgnm_count2(v)) / n_tot > (1.0 - opt->msf + 1e-9);		// check that each kmer is present in at least opt->msf*N files
			if (cond_msf1) {
				filt1 = 1;
			}
			if (cond_msf2) {
				filt2 = 1;
			}
			if (cond_msf1 && cond_msf2) { // filter the hash table entry if both k-mers should be filtered
				cond_msf = 1;
			} else { // flag the k-mer
				cond_msf = 0;
				kh_val(g->h, k) = f_val_pack(f_val_gnm_count2(v),f_val_gnm_count1(v), filt2, filt1, f_val_snp2(v), f_val_snp1(v), f_val_cb2(v), f_val_cb1(v), f_val_pgnm_count2(v), f_val_pgnm_count1(v));
			}
			cond = cond_msf || (f_val_filt1(v) && f_val_filt2(v));	// both filt bits need to be on in k-mers
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

	return n_del;
}

int64_t pg_mht_insert_list(pg_mht_t *h, int n, const ch_seq_t *a, int f)
{
	int j, mask = (1<<h->pre) - 1;
	int64_t n_ins = 0;
	pg_ht1_t *g = &h->h[a[0].h_flanks & mask]; // get hash table partition for the first (and all) k-mers
	if (n == 0) return 0;

	uint32_t gnm_cnt1, gnm_cnt2, cb1, cb2, snp1, snp2, pgnm_cnt1, pgnm_cnt2, filter1, filter2, v;
	
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
		if (absent) { // first occurrence, k-mer unknown
			++n_ins;
			gnm_cnt1 = 1; gnm_cnt2 = 0; filter1 = 0; filter2 = 0; snp1 = 0; snp2 = 0; cb1 = cb; cb2 = 0; pgnm_cnt1 = 0; pgnm_cnt2 = 0;
			kh_val(g->h, k) = f_val_pack(gnm_cnt2, gnm_cnt1, filter2, filter1, snp2, snp1, cb2, cb1, pgnm_cnt2, pgnm_cnt1);
		} else {
			v = kh_val(g->h, k);
			gnm_cnt1 = f_val_gnm_count1(v);
			gnm_cnt2 = f_val_gnm_count2(v);
			cb1 = f_val_cb1(v);
			cb2 = f_val_cb2(v);
			snp1 = f_val_snp1(v);
			snp2 = f_val_snp2(v);
			filter1 = f_val_filt1(v);
			filter2 = f_val_filt2(v);
			pgnm_cnt1 = f_val_pgnm_count1(v);
			pgnm_cnt2 = f_val_pgnm_count2(v);

			if (snp1 ^ snp2) { // already known as SNP, check if it is multi-allelic
				snp1 = 1;
				if (cb != cb1 && cb != cb2) {
					snp2 = 1; // multi-allelic SNP (do not count)
				} else {
					snp2 = 0; // bi-allelic SNP
					if (cb == cb1) {
						if (gnm_cnt1 < F_COUNTER_MAX) ++gnm_cnt1;
					}
					else  {
						if (gnm_cnt2 < F_COUNTER_MAX) {
							++gnm_cnt2;
						}
					}
				}
			} else if (snp1 & snp2) {
				snp1 = 1; snp2 = 1; // already known as multi-allelic SNP
				if (cb == cb1) {
					if (gnm_cnt1 < F_COUNTER_MAX) {
						gnm_cnt1 += 1;
					}
				}
				else if (cb == cb2) {
					if (gnm_cnt2 < F_COUNTER_MAX) {
						gnm_cnt2 += 1;
					} 
				} else {
						continue;
				}
			} else if (cb != cb1) { // newly identified SNP
					snp1 = 1; snp2 = 0;
					if (gnm_cnt2 < F_COUNTER_MAX) {
						++gnm_cnt2;
					}
					cb2 = cb; // store the second central base
			} else { // still non-SNP
				snp1 = 0; snp2 = 0;
				if (gnm_cnt1 < F_COUNTER_MAX) {
					++gnm_cnt1;
				}
			}

			kh_val(g->h, k) = f_val_pack(gnm_cnt2, gnm_cnt1, filter2, filter1, snp2, snp1, cb2, cb1, pgnm_cnt2, pgnm_cnt1);
		}
	}
	
	return n_ins;
}


void pg_mht_clear1(pg_mht_t *h, long i, int f, int max_occ)
{
	// store entries to delete
	pg_ht1_t *g = &h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g->h); ++k) {
		if (!kh_exist(g->h, k)) continue;
		uint32_t v = kh_val(g->h, k);
		uint32_t gnm_cnt1 = f_val_gnm_count1(v);
		uint32_t gnm_cnt2 = f_val_gnm_count2(v);
		uint32_t pgnm_cnt1 = f_val_pgnm_count1(v);
		uint32_t pgnm_cnt2 = f_val_pgnm_count2(v);
		uint32_t filt1 = f_val_filt1(v);
		uint32_t filt2 = f_val_filt2(v);
		uint32_t snp1 = f_val_snp1(v);
		uint32_t snp2 = f_val_snp2(v);
		uint32_t cb1 = f_val_cb1(v);
		uint32_t cb2 = f_val_cb2(v);

		// filters list
		if (f == 0) { // the mildest filter, keep everything that has counts larger than 0 and minimum than max_occ if passed
			if (gnm_cnt1 > 0 && gnm_cnt1 <= max_occ) {
				pgnm_cnt1++;
			}
			if (gnm_cnt2 > 0 && gnm_cnt2 <= max_occ) {
				pgnm_cnt2++;
			}
			if (gnm_cnt1 > max_occ) {
				filt1 = 1;
			}
			if (gnm_cnt2 > max_occ) {
				filt2 = 1;
			}
		} else if (f == 1) {
			if (gnm_cnt1 > 0 && gnm_cnt2 > 0) {
				filt1 = 1;
				filt2 = 1;
			} else if (gnm_cnt1 > 0 && gnm_cnt1 <= max_occ && gnm_cnt2 == 0) {
				pgnm_cnt1++;
			} else if (gnm_cnt1 == 0 && gnm_cnt2 > 0 && gnm_cnt2 <= max_occ) {
				pgnm_cnt2++;
			} 
			if (gnm_cnt1 > max_occ) {
				filt1 = 1;
			}
			if (gnm_cnt2 > max_occ) {
				filt2 = 1;
			}
		} else if (f == 2) { // the strictest filter, keeps only unikmers
			if (f_val_gnm_count1(v) == 1 && f_val_gnm_count2(v) == 0) {
				pgnm_cnt1++;
			} else if (f_val_gnm_count1(v) == 0 && f_val_gnm_count2(v) == 1) {
				pgnm_cnt2++;
			} else if (f_val_gnm_count1(v) == 0 && f_val_gnm_count2(v) == 0) {
				; // do nothing
			} else {
				filt1 = 1;
				filt2 = 1;
			}
		}
		kh_val(g->h, k) = f_val_pack(0, 0, filt2, filt1, snp2, snp1, cb2, cb1, pgnm_cnt2, pgnm_cnt1);
	}
}

void pg_mht_clear2(pg_mht_t *h, long i)
{
	// store entries to delete
	pg_ht1_t *g = &h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g->h); ++k) {
		if (!kh_exist(g->h, k)) continue;
		uint32_t v = kh_val(g->h, k);
		kh_val(g->h, k) = s_val_pack(0, 0, s_val_filt2(v), s_val_filt1(v), s_val_cb2(v), s_val_cb1(v));
	}

	// clear info map
    pg_im1_t *m = &h->m[i];
    for (k = 0; k < kh_end(m->m); ++k) {
        if (!kh_exist(m->m, k)) continue;
        kinfo_t *info = &kh_val(m->m, k);
        free(info->i);
    }
    pg_im_destroy(m->m);
	m->m = pg_im_init();
}


void pg_mht_tighten(pg_mht_t *h)
{
	int i;
	for (i = 0; i < 1<<h->pre; ++i) {
		pg_ht_t *g = h->h[i].h;
		uint32_t sz = kh_size(g);
		if (sz == 0) {
            pg_ht_destroy(g);
            h->h[i].h = pg_ht_init();
		}
		else if (sz * 3 < kh_capacity(g))
			pg_ht_m_resize(g, sz * 3);
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

		uint32_t cb1 = f_val_cb1(v);
		uint32_t cb2 = f_val_cb2(v);
		uint32_t filt1 = f_val_filt1(v);
		uint32_t filt2 = f_val_filt2(v);

		kh_val(g->h, k) = s_val_pack(0, 0, filt2, filt1, cb2, cb1);
	}
}


void pg_mht_count_list(pg_mht_t *h, int n, const ch_seq_t *a, ch_info_t *b)
{
	int j, mask = (1<<h->pre) - 1;
	pg_ht1_t *g;
	pg_im1_t *m;
	if (n == 0) return;

	uint32_t cnt1, cnt2, cb1, cb2, filt1, filt2, v;

	g = &h->h[a[0].h_flanks & mask]; // get hash table partition for the first (and all) k-mers.
	m = &h->m[a[0].h_flanks & mask];
	
	for (j = 0; j < n; ++j) {
		uint32_t cb = a[j].cb;
		uint64_t key;
		key = (a[j].h_flanks >> h->pre);
		
		// snp-mers pass
		khint_t k = pg_ht_get(g->h, key);
		if (k == kh_end(g->h)) continue; // not a SNP-mer, skip (in theory can be removed given lookup in ch_insert_buf)

		// add counts
		v = kh_val(g->h, k);
		cnt1 = s_val_count1(v);
		cnt2 = s_val_count2(v);
		filt1 = s_val_filt1(v);
		filt2 = s_val_filt2(v);
		cb1 = s_val_cb1(v);
		cb2 = s_val_cb2(v);

		if (cb == cb1) {
			if (cnt1 < S_COUNTER_MAX) ++cnt1;
		} else if (cb == cb2) {
			if (cnt2 < S_COUNTER_MAX) {
				cnt2 += 1;
			} 
		}

		kh_val(g->h, k) = s_val_pack(cnt2, cnt1, filt2, filt1, cb2, cb1);

		// add info
		if (b) {
			int absent;
			khint_t i = pg_im_put(m->m, key, &absent);
			kinfo_t *info = &kh_val(m->m, i);;
			if (absent) {
				memset(info, 0, sizeof(kinfo_t));
				info->n = 1;
				info->m = 1;
				MALLOC(info->i, 1);
				info->i[0].postrand = m_postrand_pack(b[j].pos, b[j].strand);
				info->i[0].seq_idx = b[j].idx;
			} else {
				if (info->n == info->m) {
					info->m = info->m < 2 ? info->m + 1 : info->m + (info->m >> 1); // grow by 1.5x, but always by >=1
					REALLOC(info->i, info->m);
				}
				int n = info->n++;
				info->i[n].postrand = m_postrand_pack(b[j].pos, b[j].strand);
				info->i[n].seq_idx = b[j].idx;
			}
		}
	}
}

pg_mht_t *pg_mht_repopulate(const char *kmer_file, pg_opt_t *opt)
{	
	pg_mht_t *h = pg_mht_init(opt->k, opt->pre);
	FILE *fp;
	char *line = NULL;
	size_t line_cap = 0;
	ssize_t len;
	int half = opt->k/2, i;
	char *left, *right;
	uint64_t x[2], mask = (1ULL<<opt->k*2) - 1, shift = (opt->k - 1) * 2;
	uint64_t hash_mask = (1ULL<<((opt->k-1)*2)) - 1; // to hash only the flanks
	int64_t n_ins = 0, n_skipped = 0;

	fp = fopen(kmer_file, "r");
	if (!fp) {
		fprintf(stderr, "[E::%s] failed to open '%s'\n", __func__, kmer_file);
		return -1;
	}

	left = (char*)malloc(half + 1);
	right = (char*)malloc(half + 1);

	while ((len = getline(&line, &line_cap, fp)) >= 0) {
		// strip trailing newline/carriage return
		while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
			line[--len] = '\0';
		char a1, a2;
		int bucket, absent;
		uint64_t h_flanks, key;
		uint32_t snp1, cb1, cb2, filt2, v;
		khint_t k;

		// parse the k-mer
		int n = (int)strlen(line);

		// check lentgh
		if (opt->snp) { // expects SNP-mers
			if (n != opt->k + 4) {
				fprintf(stderr, "[E::%s] failed to parse SNP-mer: expected %d-mer, but got %d-mer\n", __func__, opt->k, n-4);
				return -1;
			}
		} else {
			if (n != opt->k) {
				fprintf(stderr, "[E::%s] failed to parse k-mer: expected %d-mer, but got %d-mer\n", __func__, opt->k, n);
				if (n == opt->k + 4) {
					fprintf(stderr, "[E::%s] use the --snp option if passing SNP-mers\n", __func__, opt->k, n);
				}
				return -1;
			}
		}

		// left flank
		memcpy(left, line, half);
		left[half] = '\0';
		// right flank
		memcpy(right, line + half + (opt->snp ? 5 : 1), half);
		right[half] = '\0';
		if (opt->snp) { // expects SNP-mers
			// a1/a2 alleles
			char *b = line + half;
			if (b[0] != '[' || b[2] != '/' || b[4] != ']') {
				fprintf(stderr, "[E::%s] failed to parse SNP-mer: alleles not in the right format\n", __func__);
				return -1;
			}
			a1 = b[1];
			a2 = b[3];
		}
		else { // expects normal k-mers
			// a1 allele
			char *b = line + half;
			a1 = b[0];
		}

		// create k-mers hash table
		for (i = 0; i < half; ++i) {
			int c = seq_nt4_table[(uint8_t)left[i]];
			if (c >= 4) {
				fprintf(stderr, "[E::%s] found a base that is not A/a, T/t, C/c, or G/g\n", __func__);
				return -1;
			}
			x[0] = (x[0] << 2 | c) & mask;                  								// forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  								// reverse strand
		}
		{
			int c = seq_nt4_table[(uint8_t)a1];
			if (c >= 4) {
				fprintf(stderr, "[E::%s] found a base that is not A/a, T/t, C/c, or G/g\n", __func__);
				return -1;
			}
			x[0] = (x[0] << 2 | c) & mask;                  								// forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  								// reverse strand
		}
		for (i = 0; i < half; ++i) {
			int c = seq_nt4_table[(uint8_t)right[i]];
			if (c >= 4) {
				fprintf(stderr, "[E::%s] found a base that is not A/a, T/t, C/c, or G/g\n", __func__);
				return -1;
			}
			x[0] = (x[0] << 2 | c) & mask;                  								// forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  								// reverse strand
		}
		
		uint64_t y = x[0] < x[1] ? x[0] : x[1];
		uint64_t y_rev = x[0] < x[1] ? x[1] : x[0];
		uint64_t flanks = (y & ((1ULL<<(opt->k/2)*2)-1))          				// right flank from raw y
						| ((y >> ((opt->k/2+1)*2)) << ((opt->k/2)*2)); 		// left flank from raw y
		uint64_t rev_flanks = (y_rev & ((1ULL<<(opt->k/2)*2)-1))          		// right flank from raw y
						| ((y_rev >> ((opt->k/2+1)*2)) << ((opt->k/2)*2)); 	// left flank from raw y

		if (flanks == rev_flanks) {
			if (opt->verbose) {
				if(opt->snp) {
					fprintf(stderr, "[E::%s] skipped SNP-mer %s%c%s because it is palindromic\n", __func__, left, a1, right);
				} else {
					fprintf(stderr, "[E::%s] skipped k-mer %s%c%s because it is palindromic\n", __func__, left, a1, right);
				}
			}
			continue;
		}
		h_flanks = pg_hash64(flanks, hash_mask);
		bucket = h_flanks & ((1<<opt->pre) - 1);
		pg_ht1_t *g = &h->h[bucket]; // get hash table partition for the k-mer
		key = h_flanks >> opt->pre;

		k = pg_ht_put(g->h, key, &absent);
		if (absent) {
			++n_ins;
			snp1 = opt->snp ? 1:0;
			cb1 = seq_nt4_table[a1];
			cb2 = opt->snp ? seq_nt4_table[a2] : 0;
			filt2 = opt->snp ? 0 : 1;	// filt2 is always 1 to start for k-mers, not for SNP-mers
			kh_val(g->h, k) = f_val_pack(0, 0, filt2, 0, 0, snp1, cb2, cb1, 0, 0);
		} else {
			v = kh_val(g->h, k);
			cb1 = f_val_cb1(v);

			if (seq_nt4_table[a1] == cb1) {
				if (opt->verbose) {
					if(opt->snp) {
						fprintf(stderr, "[E::%s] skipped duplicated SNP-mer %s%c%s, just kept once\n", __func__, left, a1, right);
					} else {
						fprintf(stderr, "[E::%s] skipped duplicated k-mer %s%c%s, just kept once\n", __func__, left, a1, right);
					}
				}
			} else { // add two k-mers to the same key as a SNP-mer, it will be handled downstream if you are counting k-mers
				kh_val(g->h, k) = f_val_pack(0, 0, 0, 0, 0, 1, seq_nt4_table[a1], cb1, 0, 0);
			}
		}
	}

	free(line); free(left); free(right);
	fclose(fp);

	h->n_ins_tot = n_ins;

	fprintf(stderr, "[M::%s] loaded %d SNP-mers from '%s' (%d skipped)\n", __func__, n_ins, kmer_file, n_skipped);

	return h;
}

// WRITE FILES
void pg_dump_kmers(const char *fn, pg_mht_t *h, int snp)
{
	FILE *fp;
    if ((fp = strcmp(fn, "-") ? fopen(fn, "w") : stdout) == 0)
        return;

	uint64_t hash_mask = (1ULL << ((h->k - 1) * 2)) - 1;
	int mid = h->k >> 1;
	int right_off = snp ? mid + 5 : mid + 1;
	char seq[64];

	for (int i = 0; i < 1 << h->pre; ++i) {
		pg_ht1_t *g = &h->h[i];
		for (khint_t k = 0; k < kh_end(g->h); ++k) {
			if (!kh_exist(g->h, k)) continue;

			uint64_t flanks = pg_hash64_inv(((uint64_t)kh_key(g->h, k) << h->pre) | (uint64_t)i, hash_mask);
			uint32_t v = kh_val(g->h, k);
			uint32_t cb1 = f_val_cb1(v);
			uint32_t cb2 = f_val_cb2(v);

			// left flank
			for (int j = 0; j < mid; ++j)
				seq[mid - 1 - j] = nt4_seq_table[(flanks >> ((mid + j) * 2)) & 3];

			// right flank
			for (int j = 0; j < mid; ++j)
				seq[right_off + (mid - 1 - j)] = nt4_seq_table[(flanks >> (j * 2)) & 3];

			if (snp) {
				seq[mid] = '[';
				seq[mid + 1] = nt4_seq_table[cb1];
				seq[mid + 2] = '/';
				seq[mid + 3] = nt4_seq_table[cb2];
				seq[mid + 4] = ']';
				seq[h->k + 4] = '\0';
				fprintf(fp, "%s\n", seq); // dump the SNP-mer
			} else {
				if (f_val_pgnm_count1(v) > 0 && !f_val_filt1(v)) { // dump the first k-mer
					seq[mid] = nt4_seq_table[cb1];
					seq[h->k] = '\0';
					fprintf(fp, "%s\n", seq);
				}
				if (f_val_pgnm_count2(v) > 0 && !f_val_filt2(v)) { // dump the second k-mer
					seq[mid] = nt4_seq_table[cb2];
					fprintf(fp, "%s\n", seq);
				}
			}
		}
	}

	if (fp != stdout) fclose(fp);
}

void write_tsv(const char *out_fn, pg_mht_t *h, const char *gnm_fn, int snp)
{
	FILE *fp = fopen(out_fn, "w");
	if (!fp) {
		fprintf(stderr, "[M::%s] Failed to open '%s'\n", __func__, out_fn);
		return;
	}

	// get sample name
	const char *bname = strrchr(gnm_fn, '/');
	bname = bname ? bname + 1 : gnm_fn;
	char sample_name[256];
	strncpy(sample_name, bname, sizeof(sample_name) - 1);
	sample_name[sizeof(sample_name) - 1] = '\0';
	// strip .gz if present
	size_t len = strlen(sample_name);
	if (len > 3 && strcmp(sample_name + len - 3, ".gz") == 0)
		sample_name[len - 3] = '\0';
	// strip .fa, .fna, or .fasta
	static const char *fa_exts[] = { ".fasta", ".fna", ".fa", NULL };
	for (int i = 0; fa_exts[i]; i++) {
		len = strlen(sample_name);
		size_t elen = strlen(fa_exts[i]);
		if (len > elen && strcmp(sample_name + len - elen, fa_exts[i]) == 0) {
			sample_name[len - elen] = '\0';
			break;
		}
	}
	if (snp)
		fprintf(fp, "snpmer\t%s\n", sample_name);
	else
		fprintf(fp, "kmer\t%s\n", sample_name);

	int mid = h->k >> 1;
	uint64_t hash_mask = (1ULL << ((h->k - 1) * 2)) - 1;
	char seq[64];

	for (int i = 0; i < 1 << h->pre; ++i) {
		pg_ht1_t *g = &h->h[i];
		for (khint_t k = 0; k < kh_end(g->h); ++k) {
			if (!kh_exist(g->h, k)) continue;

			uint64_t key = kh_key(g->h, k);
			uint32_t v = kh_val(g->h, k);
			uint32_t cb1 = s_val_cb1(v);
			uint32_t cb2 = s_val_cb2(v);
			uint64_t flanks = pg_hash64_inv(((uint64_t)key << h->pre) | (uint64_t)i, hash_mask);
			
			if (snp) { // write SNP-mer
				// rebuild "left[a1/a2]right"
				for (int j = 0; j < mid; ++j)
					seq[mid - 1 - j] = nt4_seq_table[(flanks >> ((mid + j) * 2)) & 3]; // left
				for (int j = 0; j < mid; ++j)
					seq[mid + 5 + j] = nt4_seq_table[(flanks >> (j * 2)) & 3];         // right
				seq[mid] = '[';
				seq[mid + 1] = nt4_seq_table[cb1];
				seq[mid + 2] = '/';
				seq[mid + 3] = nt4_seq_table[cb2];
				seq[mid + 4] = ']';
				seq[h->k + 4] = '\0';

				fprintf(fp, "%s\t%u,%u\n", seq, s_val_count1(v), s_val_count2(v));
			} else { // write k-mers
				// rebuild plain k-mer: left flank[a1]right flank
				for (int j = 0; j < mid; ++j)
					seq[mid - 1 - j] = nt4_seq_table[(flanks >> ((mid + j) * 2)) & 3]; // left
				for (int j = 0; j < mid; ++j)
					seq[mid + 1 + j] = nt4_seq_table[(flanks >> (j * 2)) & 3];         // right
				seq[mid] = nt4_seq_table[cb1];
				seq[h->k] = '\0';

				fprintf(fp, "%s\t%u\n", seq, s_val_count1(v));

				// print the second k-mer if it was in the passed k-mers
				if (!s_val_filt2(v)) { // if k-mers list is scanned against new genomes and the other allele is found, but not originally in the list (so discard)
					seq[mid] = nt4_seq_table[cb2];
					fprintf(fp, "%s\t%u\n", seq, s_val_count2(v));
				}
			}
		}
	}

	fclose(fp);
}

void merge_tsvs(const char *out_fn, const char *tmpdir, const char **fa_fns, int n_fns, int n_rows)
{
	FILE **fps = malloc(n_fns * sizeof(FILE*));
	for (int i = 0; i < n_fns; ++i) {
		char p[4096];
		snprintf(p, sizeof p, "%s/gnm.%d.tsv", tmpdir, i);
		fps[i] = fopen(p, "r");
		if (!fps[i])
			fprintf(stderr, "[M::%s] Failed to open '%s'\n", __func__, p);
	}

	FILE *out;
	if ((out = strcmp(out_fn, "-") ? fopen(out_fn, "wb") : stdout) == 0)
		return;

	// header: row-label column from file 0's header, sample names from fa_fns
	char line[65536];
	for (int i = 0; i < n_fns; ++i) {
		if (fps[i]) fgets(line, sizeof line, fps[i]); // discard header line, just advance past it
		if (i == 0) {
			line[strcspn(line, "\n")] = '\0';
			char *tab = strchr(line, '\t');
			if (tab) *tab = '\0';
			fprintf(out, "%s", line); // "kmer" or "snpmer"
		}

		const char *bname = strrchr(fa_fns[i], '/');
		bname = bname ? bname + 1 : fa_fns[i];
		fprintf(out, "\t%s", bname);
	}
	fprintf(out, "\n");

	// content: one row per k-mer/SNPmer, values from each file's matching line
	for (int r = 0; r < n_rows; ++r) {
		for (int i = 0; i < n_fns; ++i) {
			if (!fps[i] || !fgets(line, 65536, fps[i])) {
				fprintf(out, "%s.", i ? "\t" : "");
				continue;
			}
			line[strcspn(line, "\n")] = '\0';

			char *tab = strchr(line, '\t');
			if (!tab) { fprintf(out, "%s.", i ? "\t" : ""); continue; }
			*tab = '\0';

			if (i == 0) fprintf(out, "%s", line); // row key, from file 0 only
			fprintf(out, "\t%s", tab + 1);        // this sample's value
		}
		fprintf(out, "\n");
	}

	for (int i = 0; i < n_fns; ++i)
		if (fps[i]) fclose(fps[i]);
	free(fps);

	if (out && out != stdout)
		fclose(out);
}