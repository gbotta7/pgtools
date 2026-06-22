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
    dst->n_done = src->n_done;
    int n = 1 << src->pre;
    dst->h = (pg_ht1_t*)calloc(n, sizeof(pg_ht1_t));
	dst->m = (pg_im1_t*)calloc(n, sizeof(pg_im1_t));
    for (int i = 0; i < n; ++i) {
		dst->h[i].h = pg_ht_copy(src->h[i].h);
		dst->m[i].m = pg_im_init();
	}
	// dst->cnames.names = NULL;
    // dst->cnames.n = 0;
    // dst->cnames.m = 0;
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
	for (i = 0; i < 1<<h->pre; ++i) {
		pg_ht_destroy(h->h[i].h); // destroy hash table for each bucket.
		pg_im1_destroy(&h->m[i]);
	}
	free(h->h); free(h->m); free(h);
}

// void *pg_mht_idx(pg_mht_t *h)
// {
// 	uint32_t n_snps = 0;

// 	for (int i = 0; i < 1<<h->pre; ++i) {
// 		pg_ht_t *g = &h->h[i];
// 		uint32_t cap = kh_capacity(g);
// 		h->id_map[i].n = cap;
// 		h->id_map[i].ids = malloc(cap * sizeof(uint32_t));
// 		memset(h->id_map[i].ids, 0xFF, cap * sizeof(uint32_t)); // UINT32_MAX = not a SNP
// 		for (khint_t k = 0; k < kh_end(g); ++k) {
// 			if (!kh_exist(g, k)) continue;
// 			uint32_t v = kh_val(g, k);
// 			h->id_map[i].ids[k] = n_snps++;
// 		}
// 	}
// }

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
				cond = (double)(n_proc - (k_val_pgnm_count1(v) + k_val_pgnm_count2(v))) / n_tot > (1.0 - min_freq + 1e-9) || (!k_val_snp1(v) || k_val_snp2(v) || k_val_filt(v));
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

int64_t pg_mht_insert_list(pg_mht_t *h, int n, const k_ch_seq_t *a, int f)
{
	int j, mask = (1<<h->pre) - 1;
	int64_t n_ins = 0;
	pg_ht1_t *g = &h->h[a[0].h_flanks & mask]; // get hash table partition for the first (and all) k-mers
	if (n == 0) return 0;

	uint32_t gnm_cnt1, gnm_cnt2, cb1, cb2, snp1, snp2, pgnm_cnt1, pgnm_cnt2, filter, v;
	
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
			gnm_cnt1 = 1; gnm_cnt2 = 0; filter = 0; snp1 = 0; snp2 = 0; cb1 = cb; cb2 = 0; pgnm_cnt1 = 0; pgnm_cnt2 = 0;
			kh_val(g->h, k) = k_val_pack(gnm_cnt2, gnm_cnt1, filter, snp2, snp1, cb2, cb1, pgnm_cnt2, pgnm_cnt1);
		} else {
			v = kh_val(g->h, k);
			gnm_cnt1 = k_val_gnm_count1(v);
			gnm_cnt2 = k_val_gnm_count2(v);
			cb1 = k_val_cb1(v);
			cb2 = k_val_cb2(v);
			snp1 = k_val_snp1(v);
			snp2 = k_val_snp2(v);
			filter = k_val_filt(v);
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

			kh_val(g->h, k) = k_val_pack(gnm_cnt2, gnm_cnt1, filter, snp2, snp1, cb2, cb1, pgnm_cnt2, pgnm_cnt1);
		}
	}
	
	return n_ins;
}


void pg_mht_clear_k(pg_mht_t *h, long i, int f)
{
	// store entries to delete
	pg_ht1_t *g = &h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g->h); ++k) {
		if (!kh_exist(g->h, k)) continue;
		uint32_t v = kh_val(g->h, k);
		uint32_t gnm_cnt1 = k_val_gnm_count1(v);
		uint32_t gnm_cnt2 = k_val_gnm_count2(v);
		// filters list
		if (f == 0) { // the mildest filter, keep everything that has counts larger than 0
			if (gnm_cnt1 > 0) {
				kh_val(g->h, k) = k_val_pack(0, 0, k_val_filt(v), k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v), k_val_pgnm_count1(v) + 1);
			}
			if (gnm_cnt2 > 0) {
				kh_val(g->h, k) = k_val_pack(0, 0, k_val_filt(v), k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v) + 1, k_val_pgnm_count1(v));
			}
			if (gnm_cnt1 == 0 && gnm_cnt2 == 0) {
				kh_val(g->h, k) = k_val_pack(0, 0, k_val_filt(v), k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v), k_val_pgnm_count1(v));
			}
		} else if (f == 1) {
			if (gnm_cnt1 > 0 && gnm_cnt2 > 0) {
				kh_val(g->h, k) = k_val_pack(0, 0, 1, k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v), k_val_pgnm_count1(v)); // set filt=1 so it gets deleted in the final filter
			} else if (gnm_cnt1 > 0 && gnm_cnt2 == 0) {
				kh_val(g->h, k) = k_val_pack(0, 0, k_val_filt(v), k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v), k_val_pgnm_count1(v) + 1);
			} else if (gnm_cnt1 == 0 && gnm_cnt2 > 0) {
				kh_val(g->h, k) = k_val_pack(0, 0, k_val_filt(v), k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v) + 1, k_val_pgnm_count1(v));
			} else {
				kh_val(g->h, k) = k_val_pack(0, 0, k_val_filt(v), k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v), k_val_pgnm_count1(v));
			}
		} else if (f == 2) { // the strictest filter, keeps only unikmers
			if (k_val_gnm_count1(v) == 1 && k_val_gnm_count2(v) == 0) {
				kh_val(g->h, k) = k_val_pack(0, 0, k_val_filt(v), k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v), k_val_pgnm_count1(v) + 1);
			} else if (k_val_gnm_count1(v) == 0 && k_val_gnm_count2(v) == 1) {
				kh_val(g->h, k) = k_val_pack(0, 0, k_val_filt(v), k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v) + 1, k_val_pgnm_count1(v));
			} else if (k_val_gnm_count1(v) == 0 && k_val_gnm_count2(v) == 0) {
				kh_val(g->h, k) = k_val_pack(0, 0, k_val_filt(v), k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v), k_val_pgnm_count1(v));
			} else {
				kh_val(g->h, k) = k_val_pack(0, 0, 1, k_val_snp2(v), k_val_snp1(v), k_val_cb2(v), k_val_cb1(v), k_val_pgnm_count2(v), k_val_pgnm_count1(v)); // set filt=1 so it gets deleted in the final filter
			}
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

		uint32_t cb1 = k_val_cb1(v);
		uint32_t cb2 = k_val_cb2(v);
		uint32_t snp1 = k_val_snp1(v);
		uint32_t snp2 = k_val_snp2(v);

		kh_val(g->h, k) = s_val_pack(0, 0, snp2, snp1, cb2, cb1);
	}
}


void pg_mht_count_list(pg_mht_t *h, int n, const s_ch_seq_t *a)
{
	int j, mask = (1<<h->pre) - 1;
	pg_ht1_t *g;
	pg_im1_t *m;
	if (n == 0) return;

	uint32_t cnt1, cnt2, cb1, cb2, snp1, snp2, v;

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
		cb1 = s_val_cb1(v);
		cb2 = s_val_cb2(v);
		snp1 = s_val_snp1(v);
		snp2 = s_val_snp2(v);

		if (cb == cb1) {
			if (cnt1 < S_COUNTER_MAX) ++cnt1;
		} else if (cb == cb2) {
			if (cnt2 < S_COUNTER_MAX) {
				cnt2 += 1;
			} 
		}

		kh_val(g->h, k) = s_val_pack(cnt2, cnt1, snp2, snp1, cb2, cb1);

		// add info
		int absent;
		khint_t i = pg_im_put(m->m, key, &absent);
		kinfo_t *info = &kh_val(m->m, i);;
		if (absent) {
			memset(info, 0, sizeof(kinfo_t));
			info->n = 1;
			info->m = 1;
			MALLOC(info->i, 1);
			info->i[0].pos = a[j].pos;
			info->i[0].strand = a[j].strand;
			info->i[0].seq_idx = a[j].idx;
		} else {
			if (info->n == info->m) {
				info->m = info->m + (info->m >> 1); // grow by 1.5x
				REALLOC(info->i, info->m);
			}
			int n = info->n++;
			info->i[n].pos = a[j].pos;
			info->i[n].strand = a[j].strand;
			info->i[n].seq_idx = a[j].idx;
		}

	}
}

void pg_dump_snps(const char *fn, pg_mht_t *h) {
	FILE *fp;
    if ((fp = strcmp(fn, "-") ? fopen(fn, "w") : stdout) == 0)
        return;

	uint64_t hash_mask = (1ULL << ((h->k - 1) * 2)) - 1;
    uint64_t flanks;
	uint32_t v, cb1, cb2;
	char seq[h->k + 1];
	int mid = h->k >> 1;
	
	for (int i = 0; i < 1<<h->pre; ++i) {
        pg_ht1_t *g = &h->h[i];
        for (khint_t k = 0; k < kh_end(g->h); ++k) {
            if (!kh_exist(g->h, k)) continue;
            flanks = pg_hash64_inv(((uint64_t)kh_key(g->h, k) << h->pre) | (uint64_t)i, hash_mask);
			v = kh_val(g->h, k);
			cb1 = k_val_cb1(v);
			cb2 = k_val_cb2(v);

			for (int j = 0; j < mid; ++j)
				seq[h->k - 1 - j] = nt4_seq_table[(flanks >> (j * 2)) & 3];
			for (int j = 0; j < mid; ++j)
				seq[mid - 1 - j] = nt4_seq_table[(flanks >> ((mid + j) * 2)) & 3];
			seq[h->k] = '\0';

			seq[mid] = nt4_seq_table[cb1];
			// fprintf(fp, ">seq%u_ref\n", h->id_map[i].ids[k]);
			// fprintf(fp, ">%s\n", seq);
			fprintf(fp, "%s\n", seq);

			seq[mid] = nt4_seq_table[cb2];
			// fprintf(fp, ">seq%u_alt\n", h->id_map[i].ids[k]);
			// fprintf(fp, ">%s\n", seq);
			fprintf(fp, "%s\n", seq);

		}
	}	
}


// VCF writer
void write_vcf(const char *out_fn, pg_mht_t *h, pg_mht_t *ref_h, char *gnm_fn)
{
    FILE *fp = fopen(out_fn, "w");
    if (!fp) {
		fprintf(stderr, "[M::%s] Failed to open '%s'\n", __func__, out_fn);
        return;
    }

    // VCF header
    fprintf(fp,
		"##fileformat=VCFv4.2\n"
		"##INFO=<ID=KPOS,Number=1,Type=String,Description=\"K-mer positions: genome|strand|pos,... per genome separated by ;\">\n"
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
		"##FORMAT=<ID=KC,Number=2,Type=Integer,Description=\"K-mer counts for REF and ALT alleles\">\n");
	

	// collect unique contigs with a hash table
	strset_t *seen = strset_init();

	for (int i = 0; i < 1<<h->pre; ++i) {
		pg_im1_t *ref_m = &ref_h->m[i];
        khint_t k;

		for (k = 0; k < kh_end(ref_m->m); ++k) {
			if (!kh_exist(ref_m->m, k)) continue;

			int absent;
			const char *ref_cname = ref_h->cnames.names[kh_val(ref_m->m, k).i[0].seq_idx];
			strset_put(seen, ref_cname, &absent);
			if (absent)
				fprintf(fp, "##contig=<ID=%s>\n", ref_cname);
		}
	}
	fprintf(fp, "##contig=<ID=.>\n"); // add contig not found in reference

	// strip path and extension from gnm_fn for sample name
	const char *bname = strrchr(gnm_fn, '/');
	bname = bname ? bname + 1 : gnm_fn;
	char sample_name[256];
	strncpy(sample_name, bname, sizeof(sample_name) - 1);
	sample_name[sizeof(sample_name) - 1] = '\0';
	char *dot = strrchr(sample_name, '.');
	if (dot) *dot = '\0';

	fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", sample_name);

	// VCF content
	char *chrom_name;
	int target_pos;
    for (int i = 0; i < 1<<h->pre; ++i) {
        pg_ht1_t *g = &h->h[i];
		pg_im1_t *m = &h->m[i];
		pg_im1_t *ref_m = &ref_h->m[i];
        khint_t k;

        for (k = 0; k < kh_end(g->h); ++k) {
            if (!kh_exist(g->h, k)) continue;

            uint32_t v = kh_val(g->h, k);
            uint32_t cb1 = s_val_cb1(v);
            uint32_t cb2 = s_val_cb2(v);
            uint32_t cnt1 = s_val_count1(v);
            uint32_t cnt2 = s_val_count2(v); /// ADD REF HERE SO YOU KNOW WHICH ONE IS THE REF

            char ref = nt4_seq_table[cb1];
            char alt = nt4_seq_table[cb2];

            const char *gt;
			if (cnt1 > 0 && cnt2 == 0) gt = "0";
			else if (cnt1 == 0 && cnt2 > 0) gt = "1";
			else gt = ".";

			// get CHROM/POS from reference info map
			uint64_t key = kh_key(g->h, k);
			chrom_name = ".";
			target_pos = 0;
			khint_t ref_info_k = pg_im_get(ref_m->m, key);
			if (ref_info_k != kh_end(ref_m->m)) {
				kinfo_t *ref_occ = &kh_val(ref_m->m, ref_info_k);
				if (ref_occ->n > 0) {
					chrom_name = ref_h->cnames.names[ref_occ->i[0].seq_idx];
					target_pos = ref_occ->i[0].pos;
				}
			}

			// build INFO field: seq_name,strand,pos;seq_name,strand,pos;...
			char info[65536];
			int info_len = 0;

			khint_t info_k = pg_im_get(m->m, key);
			if (info_k != kh_end(m->m)) {
				kinfo_t *occ = &kh_val(m->m, info_k);
				for (int j = 0; j < occ->n; ++j) {
					const char *name = h->cnames.names[occ->i[j].seq_idx];
					khint_t ref_seen_k = strset_get(seen, name);
					if (ref_seen_k != kh_end(seen) && occ->n == 1) {
						info[0] = '.'; info[1] = '\0'; info_len = 1;
					} else if (ref_seen_k != kh_end(seen) && occ->n > 1) {
						if (j) {
							info_len += snprintf(info + info_len, sizeof(info) - info_len,
												"%s|%c|%d,",
												name, occ->i[j].strand, occ->i[j].pos);
						}
					} else {
						info_len += snprintf(info + info_len, sizeof(info) - info_len,
											"%s|%c|%d,",
											name, occ->i[j].strand, occ->i[j].pos);
					}
				}
			}			
			if (info_len == 0) {
				info[0] = '.'; info[1] = '\0';
			} else {
				if (info[info_len - 1] == ',') info[--info_len] = '\0';
				// prepend KPOS=
				char tmp[65536];
				snprintf(tmp, sizeof(tmp), "KPOS=%s", info);
				memcpy(info, tmp, strlen(tmp) + 1);
			}

            fprintf(fp,
					"%s\t%d\t.\t%c\t%c\t.\t.\t%s\tGT:KC\t%s:%u,%u\n",
					chrom_name, target_pos, ref, alt,
					info, gt, cnt1, cnt2);
		}
    }

	strset_destroy(seen);
    fclose(fp);
}

void merge_vcfs(const char *out_fn, const char *tmpdir, int n_fns, int n_snps)
{
    FILE **fps = malloc(n_fns * sizeof(FILE*));
    for (int i = 0; i < n_fns; ++i) {
        char p[4096];
        snprintf(p, sizeof p, "%s/gnm.%d.vcf", tmpdir, i);
        fps[i] = fopen(p, "r");
        if (!fps[i])
            fprintf(stderr, "[M::%s] Failed to open '%s'\n", __func__, p);
    }

	FILE *out;
	if ((out = strcmp(out_fn, "-") ? fopen(out_fn, "wb") : stdout) == 0)
		return -1;

    // VCF header — meta-lines from file 0 only
    char **sample_names = malloc(n_fns * sizeof(char*));
    char line[65536];
    for (int i = 0; i < n_fns; ++i) {
        sample_names[i] = NULL;
        while (fgets(line, sizeof line, fps[i])) {
            if (strncmp(line, "#CHROM", 6) == 0) {
                char *tok = strrchr(line, '\t');
                if (tok) {
                    tok++;
                    tok[strcspn(tok, "\n")] = '\0';
                    sample_names[i] = strdup(tok);
                }
                break;
            }
            if (i == 0 && line[0] == '#')
                fputs(line, out);
        }
    }

    fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (int i = 0; i < n_fns; ++i)
        fprintf(out, "\t%s", sample_names[i] ? sample_names[i] : "UNKNOWN");
    fprintf(out, "\n");

    // VCF content
    char **lines = malloc(n_fns * sizeof(char*));
    for (int i = 0; i < n_fns; ++i)
        lines[i] = malloc(65536);

    for (int s = 0; s < n_snps; ++s) {
        for (int i = 0; i < n_fns; ++i) {
            if (!fps[i] || !fgets(lines[i], 65536, fps[i]))
                lines[i][0] = '\0';
            lines[i][strcspn(lines[i], "\n")] = '\0';
        }

        char *p = lines[0];
        char *fields[10] = {0};
        char *tmp = lines[0];
        int fc = 0;
        fields[fc++] = tmp;
        for (char *c = tmp; *c && fc < 10; ++c) {
            if (*c == '\t') {
                *c = '\0';
                fields[fc++] = c + 1;
            }
        }

        // merge INFO from all files, stripping KPOS= prefix and re-adding once
        char merged_info[65536];
        int info_len = 0;
        for (int i = 0; i < n_fns; ++i) {
            char *sp = lines[i];
            int tc = 0;
            char *info_start = NULL, *info_end = NULL;
            for (char *c = sp; *c; ++c) {
                if (*c == '\t') {
                    ++tc;
                    if (tc == 7) info_start = c + 1;
                    if (tc == 8) { info_end = c; break; }
                }
            }
            if (info_start && info_end && info_end > info_start) {
                // strip KPOS= prefix
                if (strncmp(info_start, "KPOS=", 5) == 0) info_start += 5;
                int len = info_end - info_start;
                if (!(len == 1 && *info_start == '.')) {
                    if (info_len > 0) merged_info[info_len++] = ',';
                    memcpy(merged_info + info_len, info_start, len);
                    info_len += len;
                }
            }
        }
        if (info_len == 0) { merged_info[0] = '.'; merged_info[1] = '\0'; }
        else merged_info[info_len] = '\0';

        // write output line with KPOS= prefix on merged INFO
        if (merged_info[0] == '.')
            fprintf(out, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t%s\t%s",
                    fields[0], fields[1], fields[2], fields[3],
                    fields[4], fields[5], fields[6],
                    fields[8], fields[9]);
        else
            fprintf(out, "%s\t%s\t%s\t%s\t%s\t%s\t%s\tKPOS=%s\t%s\t%s",
                    fields[0], fields[1], fields[2], fields[3],
                    fields[4], fields[5], fields[6],
                    merged_info, fields[8], fields[9]);

        // append sample columns from remaining files
        for (int i = 1; i < n_fns; ++i) {
            char *sp = lines[i];
            int tc = 0;
            for (char *c = sp; *c; ++c) {
                if (*c == '\t' && ++tc == 9) {
                    fprintf(out, "\t%s", c + 1);
                    break;
                }
            }
        }
        fprintf(out, "\n");
    }

    for (int i = 0; i < n_fns; ++i) {
        if (fps[i]) fclose(fps[i]);
        free(lines[i]);
        free(sample_names[i]);
    }
    free(fps); free(lines); free(sample_names);

	if (out && out != stdout)
        fclose(out);
    // fclose(out);
}