#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include "shtab.h"
#include "khashl.h" // hash table
#include "kseq.h"
#include "utils.h"

// Operations on hash tables

pg_msht_t *pg_msht_init(int k, int pre, int w)
{
	pg_msht_t *h;
	int i;
	CALLOC(h, 1);
	h->k = k;
	h->pre = pre;
	pthread_mutex_init(&h->mutex, 0);
	// allocate the array of partitions
	if(w) {
		CALLOC(h->ih, 1<<h->pre);
	} else {
		CALLOC(h->h, 1<<h->pre);
	}
	for (i = 0; i < 1<<h->pre; ++i) {
		// initialize hash table for each bucket
		if(w) {
			h->ih[i].ih = pg_siht_init();
		} else {
			h->h[i].h = pg_sht_init();
		}
	}
	return h;
}

void pg_siht1_destroy(pg_siht1_t *ih)
{
    khint_t k;
    for (k = 0; k < kh_end(ih->ih); ++k) {
        if (!kh_exist(ih->ih, k)) continue;
        sci_t *info = &kh_val(ih->ih, k);
        free(info->i);
    }
    pg_siht_destroy(ih->ih);
}

void pg_msht_destroy(pg_msht_t *h, int w)
{
	int i;
	if (h == 0) return;
	pthread_mutex_destroy(&h->mutex);
	for (i = 0; i < 1<<h->pre; ++i) {
		// destroy hash table for each bucket.
		if (w) {
			pg_siht1_destroy(&h->ih[i]);
		} else {
			pg_sht_destroy(h->h[i].h);
		}
		
	}
	free(h->h); free(h->ih); free(h);
}

int64_t pg_msht_filter(pg_msht_t *h, long i, int n_proc, int n_tot, int ff, pg_opt_t *opt)
{	
	int64_t n_del = 0;
	int cond = 0, cond_msf = 0, cond_maf = 0, cond_snp = 0;

	pg_sht1_t *g = &h->h[i];
	uint64_t *del_part = malloc(kh_size(g->h) * sizeof(uint64_t)), kv;
	khint_t k;
	uint32_t v;
	for (k = 0; k < kh_end(g->h); ++k) {
		if (!kh_exist(g->h, k)) continue;
		kv = kh_key(g->h, k);
		v = kh_val(g->h, k);
		cond_msf = (double)(n_proc - (f_val_pgnm_count1(v) + f_val_pgnm_count2(v))) / n_tot > (1.0 - opt->msf + 1e-9);		// check that one of the two alleles is present in at least opt->msf*N files
		if (ff) cond_maf = (f_val_pgnm_count1(v) / (double)(f_val_pgnm_count1(v) + f_val_pgnm_count2(v))) < opt->maf || (f_val_pgnm_count2(v) / (double)(f_val_pgnm_count1(v) + f_val_pgnm_count2(v))) < opt->maf;		// check minimum allelic frequency > opt->maf (only in ff)
		if (ff) cond_snp = !f_val_snp1(kv) || f_val_snp2(kv);				// remove multi-allelic and non SNP-mers
		cond = cond_msf || cond_maf || cond_snp || f_val_filt(kv);			// just one of the filt bits needs to be on in SNP-mers
		if (cond) {
			del_part[n_del++] = kh_key(g->h, k);
		}
	}

	// delete entries
	for (int d = 0; d < n_del; ++d) {
		k = pg_sht_get(g->h, del_part[d]);
		if (k != kh_end(g->h)) {
			pg_sht_del(g->h, k);
		}
	}
	
	// fprintf(stderr, "n_del = %lld\n", n_del);
	free(del_part);

	return n_del;
}

int64_t pg_msht_insert_list(pg_msht_t *h, int n, const seq_t *a, int f)
{
	int j;
	uint64_t mask = (1<<h->pre) - 1;
	int64_t n_ins = 0;
	if (n == 0) return 0;
	pg_sht1_t *g = &h->h[a[0].h_flanks & mask]; // get hash table partition for the first (and all) SNP-mers

	uint32_t gnm_cnt1, gnm_cnt2, pgnm_cnt1, pgnm_cnt2, v;
	uint64_t filter, cb1, cb2, snp1, snp2, kv;
	
	for (j = 0; j < n; ++j) {
		int absent;
		uint32_t cb = a[j].cb;
		uint64_t key = (a[j].h_flanks >> h->pre);
		if ((a[j].h_flanks & mask) != (a[0].h_flanks & mask)) continue; // skip if the partition bits are different (should not happen)
		
		khint_t k;
		if (!f) {
			k = pg_sht_put(g->h, key << F_VAL_INFO_BITS, &absent);
		} else {
			k = pg_sht_get(g->h, key << F_VAL_INFO_BITS);
			if (k == kh_end(g->h)) continue; // SNP-mer not found until now, skip
			absent = 0;
		}
		if (absent) { // first occurrence, SNP-mer unknown
			++n_ins;
			gnm_cnt1 = 1; gnm_cnt2 = 0; filter = 0; snp1 = 0; snp2 = 0; cb1 = cb; cb2 = 0; pgnm_cnt1 = 0; pgnm_cnt2 = 0;
			kh_key(g->h, k) = key << F_VAL_INFO_BITS | f_key_pack(filter, snp2, snp1, cb2, cb1);
			kh_val(g->h, k) = f_val_pack(gnm_cnt2, gnm_cnt1, pgnm_cnt2, pgnm_cnt1);
		} else {
			kv = kh_key(g->h, k);
			cb1 = f_val_cb1(kv);
			cb2 = f_val_cb2(kv);
			snp1 = f_val_snp1(kv);
			snp2 = f_val_snp2(kv);
			filter = f_val_filt(kv);
			v = kh_val(g->h, k);
			gnm_cnt1 = f_val_gnm_count1(v);
			gnm_cnt2 = f_val_gnm_count2(v);
			pgnm_cnt1 = f_val_pgnm_count1(v);
			pgnm_cnt2 = f_val_pgnm_count2(v);

			if (snp1 ^ snp2) { // already known as SNP, check if it is multi-allelic
				snp1 = 1;
				if (cb != cb1 && cb != cb2) {
					snp2 = 1; // multi-allelic SNP (do not count)
				} else {
					snp2 = 0; // bi-allelic SNP
					if (cb == cb1) {
						if (gnm_cnt1 < F_GNM_COUNTER_MAX) {
							gnm_cnt1++;
						}
					}
					else  {
						if (gnm_cnt2 < F_GNM_COUNTER_MAX) {
							gnm_cnt2++;
						}
					}
				}
			} else if (snp1 & snp2) {
				snp1 = 1; snp2 = 1; // already known as multi-allelic SNP
				if (cb == cb1) {
					if (gnm_cnt1 < F_GNM_COUNTER_MAX) {
						gnm_cnt1++;
					}
				}
				else if (cb == cb2) {
					if (gnm_cnt2 < F_GNM_COUNTER_MAX) {
						gnm_cnt2++;
					} 
				} else {
					continue;
				}
			} else if (cb != cb1) { // newly identified SNP
					snp1 = 1; snp2 = 0;
					if (gnm_cnt2 < F_GNM_COUNTER_MAX) {
						++gnm_cnt2;
					}
					cb2 = cb; // store the second central base
			} else { // still non-SNP
				snp1 = 0; snp2 = 0;
				if (gnm_cnt1 < F_GNM_COUNTER_MAX) {
					++gnm_cnt1;
				}
			}

			kh_key(g->h, k) = (kv & ~F_VAL_MAX) | f_key_pack(filter, snp2, snp1, cb2, cb1);
			kh_val(g->h, k) = f_val_pack(gnm_cnt2, gnm_cnt1, pgnm_cnt2, pgnm_cnt1);
		}
	}
	
	return n_ins;
}


void pg_msht_clear1(pg_msht_t *h, long i, int f, int max_occ) // first pass
{
	// store entries to delete
	pg_sht1_t *g = &h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g->h); ++k) {
		if (!kh_exist(g->h, k)) continue;
		uint64_t kv = kh_key(g->h, k);
		uint64_t filt = f_val_filt(kv);
		uint64_t snp1 = f_val_snp1(kv);
		uint64_t snp2 = f_val_snp2(kv);
		uint64_t cb1 = f_val_cb1(kv);
		uint64_t cb2 = f_val_cb2(kv);
		uint32_t v = kh_val(g->h, k);
		uint32_t gnm_cnt1 = f_val_gnm_count1(v);
		uint32_t gnm_cnt2 = f_val_gnm_count2(v);
		uint32_t pgnm_cnt1 = f_val_pgnm_count1(v);
		uint32_t pgnm_cnt2 = f_val_pgnm_count2(v);

		// filters list
		if (f == 0) { // the mildest filter, keep everything that has counts larger than 0 and minimum than max_occ if passed
			if (gnm_cnt1 > 0 && gnm_cnt1 <= max_occ) {
				if (pgnm_cnt1 < F_PGNM_COUNTER_MAX) {
					pgnm_cnt1++;
				}
			}
			if (gnm_cnt2 > 0 && gnm_cnt2 <= max_occ) {
				if (pgnm_cnt2 < F_PGNM_COUNTER_MAX) {
					pgnm_cnt2++;
				}
			}
			if (gnm_cnt1 > max_occ || gnm_cnt2 > max_occ) {
				filt = 1;
			}
		} else if (f == 1) {
			if (gnm_cnt1 > 0 && gnm_cnt2 > 0 || (gnm_cnt1 > max_occ || gnm_cnt2 > max_occ)) {
				filt = 1;
			} else if (gnm_cnt1 > 0 && gnm_cnt1 <= max_occ && gnm_cnt2 == 0) {
				if (pgnm_cnt1 < F_PGNM_COUNTER_MAX) {
					pgnm_cnt1++;
				}
			} else if (gnm_cnt1 == 0 && gnm_cnt2 > 0 && gnm_cnt2 <= max_occ) {
				if (pgnm_cnt2 < F_PGNM_COUNTER_MAX) {
					pgnm_cnt2++;
				}
			} 
		} else if (f == 2) { // the strictest filter, keeps only unikmers
			if (f_val_gnm_count1(v) == 1 && f_val_gnm_count2(v) == 0) {
				if (pgnm_cnt1 < F_PGNM_COUNTER_MAX) {
					pgnm_cnt1++;
				}
			} else if (f_val_gnm_count1(v) == 0 && f_val_gnm_count2(v) == 1) {
				if (pgnm_cnt2 < F_PGNM_COUNTER_MAX) {
					pgnm_cnt2++;
				}
			} else if (f_val_gnm_count1(v) == 0 && f_val_gnm_count2(v) == 0) {
				; // do nothing
			} else {
				filt = 1;
			}
		}
		kh_key(g->h, k) = (kv & ~F_VAL_MAX) | f_key_pack(filt, snp2, snp1, cb2, cb1);
		kh_val(g->h, k) = f_val_pack(0, 0, pgnm_cnt2, pgnm_cnt1);
	}
}

void pg_msht_clear2(pg_msht_t *h, long i, int w) // second pass
{
	// store entries to delete
	pg_sht1_t *g = &h->h[i];
	pg_siht1_t *ig = &h->ih[i];
	khint_t k;

	if (w) {
		for (k = 0; k < kh_end(ig->ih); ++k) {
			if (!kh_exist(ig->ih, k)) continue;
			sci_t *v = &kh_val(ig->ih, k);
			kh_val(ig->ih, k).cnt = 0;
			free(v->i);
			v->i = NULL;
			v->n = v->m = 0;
		}
		pg_siht_destroy(ig->ih);
		ig->ih = pg_siht_init();
	} else {
		for (k = 0; k < kh_end(g->h); ++k) {
			if (!kh_exist(g->h, k)) continue;
			kh_val(g->h, k) = 0;
		}
	}
}


void pg_msht_tighten(pg_msht_t *h)
{
	int i;
	for (i = 0; i < 1<<h->pre; ++i) {
		pg_sht_t *g = h->h[i].h;
		uint32_t sz = kh_size(g);
		if (sz == 0) {
            pg_sht_destroy(g);
            h->h[i].h = pg_sht_init();
		}
		else if (sz * 3 < kh_capacity(g))
			pg_sht_m_resize(g, sz * 3);
	}
}


void pg_msht_rearrange(pg_msht_t *h, long i)
{
	// store entries to delete
	pg_sht1_t *g = &h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g->h); ++k) {
		if (!kh_exist(g->h, k)) continue;
		uint64_t kv = kh_key(g->h, k);

		uint64_t cb1 = f_val_cb1(kv);
		uint64_t cb2 = f_val_cb2(kv);

		kh_key(g->h, k) = (kv & ~S_VAL_MAX) | s_key_pack(cb2, cb1);
		kh_val(g->h, k) = 0;
	}
}


void pg_msht_count_list(pg_msht_t *h, int n, const seq_t *a, seq_info_t *b)
{
	int j, mask = (1<<h->pre) - 1;
	pg_sht1_t *g;
	pg_siht1_t *ig;
	if (n == 0) return;
	uint32_t cnt1, cnt2, v;
	uint64_t cb1, cb2, kv;
	sci_t *vi;

	// get hash table partition for the first (and all) SNP-mers.
	if (b) {
		ig = &h->ih[a[0].h_flanks & mask];
	} else {
		g = &h->h[a[0].h_flanks & mask];
	}
	
	for (j = 0; j < n; ++j) {
		uint64_t cb = a[j].cb;
		uint64_t key;
		if ((a[j].h_flanks & mask) != (a[0].h_flanks & mask)) continue; // skip if the partition bits are different (should not happen)
		key = (a[j].h_flanks >> h->pre);
		
		if (b) { // with info
			khint_t k = pg_siht_get(ig->ih, key << F_VAL_INFO_BITS);
			if (k == kh_end(ig->ih)) continue; // not a SNP-mer, skip

			// add counts
			vi = &kh_val(ig->ih, k);
			kv = kh_key(ig->ih, k);
			
			cnt1 = s_val_count1(vi->cnt);
			cnt2 = s_val_count2(vi->cnt);
			cb1 = s_val_cb1(kv);
			cb2 = s_val_cb2(kv);

			if (cb == cb1) {
				if (cnt1 < S_COUNTER_MAX) {
					cnt1 += 1;
				}
			} else if (cb == cb2) {
				if (cnt2 < S_COUNTER_MAX) {
					cnt2 += 1;
				} 
			}
			
			kh_val(ig->ih, k).cnt = s_val_pack(cnt2, cnt1);

			// add info
			if (vi->n == vi->m) {
				vi->m = vi->m < 2 ? vi->m + 1 : vi->m + (vi->m >> 1); // grow by 1.5x, but always by >=1
				REALLOC(vi->i, vi->m);
			}
			int idx = vi->n++;
			vi->i[idx].posallele = i_posallele_pack(b[j].pos, b[j].allele);
			vi->i[idx].seq_idx = b[j].idx;
		} else { // without info
			khint_t k = pg_sht_get(g->h, key << F_VAL_INFO_BITS);
			if (k == kh_end(g->h)) continue; // not a SNP-mer, skip

			// add counts
			v = kh_val(g->h, k);
			kv = kh_key(g->h, k);
			
			cnt1 = s_val_count1(v);
			cnt2 = s_val_count2(v);
			cb1 = s_val_cb1(kv);
			cb2 = s_val_cb2(kv);

			if (cb == cb1) {
				if (cnt1 < S_COUNTER_MAX) {
					cnt1++;
				}
			} else if (cb == cb2) {
				if (cnt2 < S_COUNTER_MAX) {
					cnt2++;
				} 
			}

			kh_val(g->h, k) = s_val_pack(cnt2, cnt1);
		}
	}
}

pg_msht_t *pg_msht_repopulate(const char *kmer_file, pg_opt_t *opt)
{	
	pg_msht_t *h = pg_msht_init(opt->k, opt->pre, opt->write_info);
	FILE *fp;
	char *line = NULL;
	size_t line_cap = 0;
	ssize_t len;
	int half = opt->k/2, i;
	char *left, *right;
	uint64_t x[2] = {0, 0}, mask = (1ULL << opt->k*2) - 1, shift = (opt->k - 1) * 2;
	uint64_t hash_mask = (1ULL << ((opt->k-1)*2)) - 1; // to hash only the flanks
	int64_t n_ins = 0, n_skipped = 0;

	fp = fopen(kmer_file, "r");
	if (!fp) {
		fprintf(stderr, "[E::%s] failed to open '%s'\n", __func__, kmer_file);
		fclose(fp);
		return NULL;
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
		uint64_t cb1, cb2;
		khint_t k;

		// parse the k-mer
		int n = (int)strlen(line);

		// check lentgh
		if (n != opt->k + 4) {
			fprintf(stderr, "[E::%s] failed to parse SNP-mer: expected %d-mer, but got %d-mer\n", __func__, opt->k, n-4);
			free(line); free(left); free(right); fclose(fp);
			return NULL;
		}

		// left flank
		memcpy(left, line, half);
		left[half] = '\0';
		// right flank
		memcpy(right, line + half + 5, half);
		right[half] = '\0';
		// a1/a2 alleles
		char *b = line + half;
		if (b[0] != '[' || b[2] != '/' || b[4] != ']') {
			fprintf(stderr, "[E::%s] failed to parse SNP-mer: alleles not in the right format\n", __func__);
			free(line); free(left); free(right); fclose(fp);
			return NULL;
		}
		a1 = b[1];
		a2 = b[3];

		// create SNP-mers hash table
		for (i = 0; i < half; ++i) {
			int c = seq_nt4_table[(uint8_t)left[i]];
			if (c >= 4) {
				fprintf(stderr, "[E::%s] found a base that is not A/a, T/t, C/c, or G/g\n", __func__);
				free(line); free(left); free(right); fclose(fp);
				return NULL;
			}
			x[0] = (x[0] << 2 | c) & mask;                  								// forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  								// reverse strand
		}

		cb1 = seq_nt4_table[(uint8_t)a1];
		cb2 = seq_nt4_table[(uint8_t)a2];
		if (cb1 >= 4 || cb2 >= 4) {
			fprintf(stderr, "[E::%s] found a base that is not A/a, T/t, C/c, or G/g\n", __func__);
			free(line); free(left); free(right); fclose(fp);
			return NULL;
		}
		x[0] = (x[0] << 2 | cb1) & mask;                  								// forward strand
		x[1] = x[1] >> 2 | (uint64_t)(3 - cb1) << shift;  								// reverse strand

		for (i = 0; i < half; ++i) {
			int c = seq_nt4_table[(uint8_t)right[i]];
			if (c >= 4) {
				fprintf(stderr, "[E::%s] found a base that is not A/a, T/t, C/c, or G/g\n", __func__);
				free(line); free(left); free(right); fclose(fp);
				return NULL;
			}
			x[0] = (x[0] << 2 | c) & mask;                  								// forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  								// reverse strand
		}
		
		// re-do because user could load non-canonical SNP-mers
		uint64_t y = x[0] < x[1] ? x[0] : x[1];
		uint64_t y_rev = x[0] < x[1] ? x[1] : x[0];
		uint64_t flanks = (y & ((1ULL<<(opt->k/2)*2)-1))          				// right flank from raw y
						| ((y >> ((opt->k/2+1)*2)) << ((opt->k/2)*2)); 		// left flank from raw y
		uint64_t rev_flanks = (y_rev & ((1ULL<<(opt->k/2)*2)-1))          		// right flank from raw y
						| ((y_rev >> ((opt->k/2+1)*2)) << ((opt->k/2)*2)); 	// left flank from raw y

		if (flanks == rev_flanks) {
			if (opt->verbose) {
				fprintf(stderr, "[E::%s] skipped SNP-mer %s%c%s because it is palindromic\n", __func__, left, a1, right);
			}
			continue;
		}

		h_flanks = pg_hash64(flanks, hash_mask);
		bucket = h_flanks & ((1<<opt->pre) - 1);
		key = h_flanks >> opt->pre;

		if (opt->write_info) {
			pg_siht1_t *ig = &h->ih[bucket];
			k = pg_siht_put(ig->ih, key << F_VAL_INFO_BITS, &absent);
			if (absent) {
				++n_ins;
				kh_key(ig->ih, k) = key << F_VAL_INFO_BITS | s_key_pack(cb2, cb1);
				sci_t *v = &kh_val(ig->ih, k);
				v->i = NULL;
				v->n = 0;
				v->m = 0;
				v->cnt = 0;
			} else {
				n_skipped++;
				if (opt->verbose) {
					fprintf(stderr, "[E::%s] skipped duplicated SNP-mer %s%c%s, just kept once\n", __func__, left, a1, right);
				}
			}
		} else {
			pg_sht1_t *g = &h->h[bucket];
			k = pg_sht_put(g->h, key << F_VAL_INFO_BITS, &absent);
			if (absent) {
				++n_ins;
				kh_key(g->h, k) = key << F_VAL_INFO_BITS | f_key_pack(0, 0, 0, cb2, cb1);
				kh_val(g->h, k) = 0;
			} else {
				n_skipped++;
				if (opt->verbose) {
					fprintf(stderr, "[E::%s] skipped duplicated SNP-mer %s%c%s, just kept once\n", __func__, left, a1, right);
				}
			}
		}
	}

	free(line); free(left); free(right); fclose(fp);

	h->n_ins_tot = n_ins;

	fprintf(stderr, "[M::%s] loaded %lld SNP-mers from '%s' (%lld skipped)\n", __func__, n_ins, kmer_file, n_skipped);

	return h;
}

// WRITE FILES
void pg_dump_snpmers(const char *fn, pg_msht_t *h)
{
	FILE *fp;
    if ((fp = strcmp(fn, "-") ? fopen(fn, "w") : stdout) == 0)
        return;

	uint64_t hash_mask = (1ULL << ((h->k - 1) * 2)) - 1;
	int mid = h->k >> 1;
	int right_off = mid + 5;
	char seq[64];

	for (int i = 0; i < 1 << h->pre; ++i) {
		pg_sht1_t *g = &h->h[i];
		for (khint_t k = 0; k < kh_end(g->h); ++k) {
			if (!kh_exist(g->h, k)) continue;

			uint64_t key = kh_key(g->h, k);
			uint64_t flanks = pg_hash64_inv((((key >> F_VAL_INFO_BITS) << h->pre)) | (uint64_t)i, hash_mask);
			uint64_t cb1 = f_val_cb1(key);
			uint64_t cb2 = f_val_cb2(key);

			// left flank
			for (int j = 0; j < mid; ++j)
				seq[mid - 1 - j] = nt4_seq_table[(flanks >> ((mid + j) * 2)) & 3];

			// right flank
			for (int j = 0; j < mid; ++j)
				seq[right_off + (mid - 1 - j)] = nt4_seq_table[(flanks >> (j * 2)) & 3];

			seq[mid] = '[';
			seq[mid + 1] = nt4_seq_table[cb1];
			seq[mid + 2] = '/';
			seq[mid + 3] = nt4_seq_table[cb2];
			seq[mid + 4] = ']';
			seq[h->k + 4] = '\0';
			fprintf(fp, "%s\n", seq); // dump the SNP-mer
		}
	}

	if (fp != stdout) fclose(fp);
}

void write_snpmer_tsv(const char *out_fn, pg_msht_t *h, const char *gnm_fn, int w, int w_mko)
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
	if (w) {
		fprintf(fp, "snpmer\t%s\tpositions\n", sample_name);
	} else {
		fprintf(fp, "snpmer\t%s\n", sample_name);
	}
	

	pg_siht1_t *ig = NULL;
	pg_sht1_t *g = NULL;
	int mid = h->k >> 1;
	int right_off = mid + 5;
	uint64_t hash_mask = (1ULL << ((h->k - 1) * 2)) - 1, key, flanks, cb1, cb2;
	sci_t *vi;
	uint32_t cnt1, cnt2, v;
	char seq[64];

	for (int i = 0; i < 1 << h->pre; ++i) {
		if (w) {
			ig = &h->ih[i];
		} else {
			g = &h->h[i];
		}

		khint_t kend = w ? kh_end(ig->ih) : kh_end(g->h);
		for (khint_t k = 0; k < kend; ++k) {
			if (w) {
				if (!kh_exist(ig->ih, k)) continue;
			} else {
				if (!kh_exist(g->h, k)) continue;
			}

			if (w) {
				key = kh_key(ig->ih, k);
				flanks = pg_hash64_inv((((key >> F_VAL_INFO_BITS) << h->pre)) | (uint64_t)i, hash_mask);
				cb1 = s_val_cb1(key);
				cb2 = s_val_cb2(key);
				vi = &kh_val(ig->ih, k);
				cnt1 = s_val_count1(vi->cnt);
				cnt2 = s_val_count2(vi->cnt);
			} else {
				key = kh_key(g->h, k);
				flanks = pg_hash64_inv((((key >> F_VAL_INFO_BITS) << h->pre)) | (uint64_t)i, hash_mask);
				cb1 = s_val_cb1(key);
				cb2 = s_val_cb2(key);
				v = kh_val(g->h, k);
				cnt1 = s_val_count1(v);
				cnt2 = s_val_count2(v);
			}

			// write SNP-mer
			// rebuild "left[a1/a2]right"
			for (int j = 0; j < mid; ++j)
				seq[mid - 1 - j] = nt4_seq_table[(flanks >> ((mid + j) * 2)) & 3]; 			// left
			for (int j = 0; j < mid; ++j)
				seq[right_off + (mid - 1 - j)] = nt4_seq_table[(flanks >> (j * 2)) & 3];	// right
			seq[mid] = '[';
			seq[mid + 1] = nt4_seq_table[cb1];
			seq[mid + 2] = '/';
			seq[mid + 3] = nt4_seq_table[cb2];
			seq[mid + 4] = ']';
			seq[h->k + 4] = '\0';

			fprintf(fp, "%s\t%u,%u", seq, cnt1, cnt2);

			if (w) {
				int occ = vi->n < w_mko ? vi->n : w_mko;
				for (int j = 0; j < occ; ++j) {
					uint32_t posallele = vi->i[j].posallele;
					uint32_t pos = i_val_pos(posallele);
					uint32_t allele = i_val_allele(posallele);
					int32_t idx = vi->i[j].seq_idx;
					const char *cname = h->cnames.names[idx];

					fprintf(fp, "%c%s:%c:%u", j == 0 ? '\t' : ',', cname, nt4_seq_table[allele], pos);
				}
			}
			fprintf(fp, "\n");
		}
	}

	fclose(fp);
}

// void merge_tsvs(const char *out_fn, const char *tmpdir, const char **fa_fns, int n_fns, int n_rows)
// {
// 	FILE **fps = malloc(n_fns * sizeof(FILE*));
// 	for (int i = 0; i < n_fns; ++i) {
// 		char p[4096];
// 		snprintf(p, sizeof p, "%s/gnm.%d.tsv", tmpdir, i);
// 		fps[i] = fopen(p, "r");
// 		if (!fps[i])
// 			fprintf(stderr, "[M::%s] Failed to open '%s'\n", __func__, p);
// 	}

// 	FILE *out;
// 	if ((out = strcmp(out_fn, "-") ? fopen(out_fn, "wb") : stdout) == 0)
// 		return;

// 	// header: row-label column from file 0's header, sample names from fa_fns
// 	char line[65536];
// 	for (int i = 0; i < n_fns; ++i) {
// 		if (fps[i]) fgets(line, sizeof line, fps[i]); // discard header line, just advance past it
// 		if (i == 0) {
// 			line[strcspn(line, "\n")] = '\0';
// 			char *tab = strchr(line, '\t');
// 			if (tab) *tab = '\0';
// 			fprintf(out, "%s", line); // "kmer" or "snpmer"
// 		}

// 		const char *bname = strrchr(fa_fns[i], '/');
// 		bname = bname ? bname + 1 : fa_fns[i];
// 		fprintf(out, "\t%s", bname);
// 	}
// 	fprintf(out, "\n");

// 	// content: one row per k-mer/SNPmer, values from each file's matching line
// 	for (int r = 0; r < n_rows; ++r) {
// 		for (int i = 0; i < n_fns; ++i) {
// 			if (!fps[i] || !fgets(line, 65536, fps[i])) {
// 				fprintf(out, "%s.", i ? "\t" : "");
// 				continue;
// 			}
// 			line[strcspn(line, "\n")] = '\0';

// 			char *tab = strchr(line, '\t');
// 			if (!tab) { fprintf(out, "%s.", i ? "\t" : ""); continue; }
// 			*tab = '\0';

// 			if (i == 0) fprintf(out, "%s", line); // row key, from file 0 only
// 			fprintf(out, "\t%s", tab + 1);        // this sample's value
// 		}
// 		fprintf(out, "\n");
// 	}

// 	for (int i = 0; i < n_fns; ++i)
// 		if (fps[i]) fclose(fps[i]);
// 	free(fps);

// 	if (out && out != stdout)
// 		fclose(out);
// }