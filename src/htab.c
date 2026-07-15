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

int64_t pg_mht_insert_list(pg_mht_t *h, int n, const ch_seq_t *a, int f)
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


void pg_mht_count_list(pg_mht_t *h, int n, const ch_seq_t *a, ch_info_t *b)
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
					info->m = info->m + (info->m >> 1); // grow by 1.5x
					REALLOC(info->i, info->m);
				}
				int n = info->n++;
				info->i[n].postrand = m_postrand_pack(b[j].pos, b[j].strand);
				info->i[n].seq_idx = b[j].idx;
			}
		}
	}
}


// VCF writer
void write_vcf(const char *out_fn, pg_mht_t *h, pg_mht_t *ref_h, char *gnm_fn, int write_info)
{
    FILE *fp = fopen(out_fn, "w");
    if (!fp) {
		fprintf(stderr, "[M::%s] Failed to open '%s'\n", __func__, out_fn);
        return;
    }

    // VCF header
	fprintf(fp,
			"##fileformat=VCFv4.2\n"
			"##source=pgtools\n"
			"##INFO=<ID=KPOS,Number=1,Type=String,Description=\"K-mer hit positions across the pangenome, formatted as contig|strand|pos. Multiple hits across the pangenome are comma-separated. Strand is + or - relative to each genome, position is always reported relative to the + strand.\">\n"
			"##INFO=<ID=SNPMER,Number=1,Type=String,Description=\"SNP-mer context around the variant, relative to the + strand.\">\n"
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Haplotype genotype call based on k-mer presence: 0=REF allele only, 1=ALT allele only, .=ambiguous (both or neither allele k-mers found; inspect KC for counts).\">\n"
			"##FORMAT=<ID=KC,Number=2,Type=Integer,Description=\"K-mer counts for the REF and ALT alleles respectively. Non zero counts on both alleles (GT=.) may indicate a copy number variant or a repetitive region.\">\n");
	

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
	char *ref_chrom_name;
	int ref_pos;
	int ref_strand;
	uint32_t cnt_ref;
	uint32_t cnt_alt;
    for (int i = 0; i < 1<<h->pre; ++i) {
        pg_ht1_t *g = &h->h[i];
		pg_im1_t *m = &h->m[i];
		pg_ht1_t *ref_g = &ref_h->h[i];
		pg_im1_t *ref_m = &ref_h->m[i];
        khint_t k;

        for (k = 0; k < kh_end(g->h); ++k) {
            if (!kh_exist(g->h, k)) continue;

			uint64_t key = kh_key(g->h, k);

			// find the reference allele
			khint_t k_ref = pg_ht_get(ref_g->h, key);
			uint32_t v_ref = kh_val(ref_g->h, k_ref);
			uint32_t cb1_ref = s_val_cb1(v_ref);
			uint32_t cb2_ref = s_val_cb2(v_ref);
			uint32_t cnt1_ref = s_val_count1(v_ref);
            uint32_t cnt2_ref = s_val_count2(v_ref);

			uint32_t cb_ref;
			if (cnt1_ref > cnt2_ref) { // keep as a reference the one that appear more frequently if both appear
				cb_ref = cb1_ref;
			} else {
				cb_ref = cb2_ref;
			}

            uint32_t v = kh_val(g->h, k);
            uint32_t cb1 = s_val_cb1(v);
            uint32_t cb2 = s_val_cb2(v);
            uint32_t cnt1 = s_val_count1(v);
            uint32_t cnt2 = s_val_count2(v);
			
			// get CHROM/POS and strand from reference info map
			ref_chrom_name = ".";
			ref_pos = 0;
			ref_strand = 0;
			khint_t ref_info_k = pg_im_get(ref_m->m, key);
			if (ref_info_k != kh_end(ref_m->m)) {
				kinfo_t *ref_occ = &kh_val(ref_m->m, ref_info_k);
				if (ref_occ->n > 0) {
					ref_chrom_name = ref_h->cnames.names[ref_occ->i[0].seq_idx];
					ref_pos = m_val_pos(ref_occ->i[0].postrand) + 1; // VCF is 1-based
					ref_strand = m_val_strand(ref_occ->i[0].postrand);
				}
			}
			
			// get REF/ALT and build SNPMER INFO field based on ref_strand
			char snpmer[64];
			uint64_t hash_mask = (1ULL << ((h->k - 1) * 2)) - 1;
			uint64_t flanks = pg_hash64_inv(((uint64_t)kh_key(ref_g->h, k) << ref_h->pre) | (uint64_t)i, hash_mask);
			int mid = ref_h->k >> 1;
			char ref, alt;
			const char *gt;
			if (cb1 == cb_ref) {
				if (ref_strand) { // already on forward strand
					ref = nt4_seq_table[cb1];
					alt = nt4_seq_table[cb2];
					for (int j = 0; j < mid; ++j)
						snpmer[ref_h->k - 1 - j] = nt4_seq_table[(flanks >> (j * 2)) & 3];
					for (int j = 0; j < mid; ++j)
						snpmer[mid - 1 - j] = nt4_seq_table[(flanks >> ((mid + j) * 2)) & 3];
					snpmer[mid] = ref;
					snpmer[ref_h->k] = '\0';
				} else { // on reverse strand, change it for VCF standard
					ref = nt4_seq_table[3 - cb1];
					alt = nt4_seq_table[3 - cb2];
					for (int j = 0; j < mid; ++j)
						snpmer[j] = nt4_seq_table[3 - ((flanks >> (j * 2)) & 3)];
					for (int j = 0; j < mid; ++j)
						snpmer[mid + 1 + j] = nt4_seq_table[3 - ((flanks >> ((mid + j) * 2)) & 3)];
					snpmer[mid] = ref;
					snpmer[ref_h->k] = '\0';
				}
				cnt_ref = cnt1;
				cnt_alt = cnt2;
			} else {
				if (ref_strand) { // already on forward strand
					ref = nt4_seq_table[cb2];
					alt = nt4_seq_table[cb1];
					for (int j = 0; j < mid; ++j)
						snpmer[ref_h->k - 1 - j] = nt4_seq_table[(flanks >> (j * 2)) & 3];
					for (int j = 0; j < mid; ++j)
						snpmer[mid - 1 - j] = nt4_seq_table[(flanks >> ((mid + j) * 2)) & 3];
					snpmer[mid] = ref;
					snpmer[ref_h->k] = '\0';
				} else { // on reverse strand, change it for VCF standard
					ref = nt4_seq_table[3 - cb2];
					alt = nt4_seq_table[3 - cb1];
					for (int j = 0; j < mid; ++j)
						snpmer[j] = nt4_seq_table[3 - ((flanks >> (j * 2)) & 3)];
					for (int j = 0; j < mid; ++j)
						snpmer[mid + 1 + j] = nt4_seq_table[3 - ((flanks >> ((mid + j) * 2)) & 3)];
					snpmer[mid] = ref;
					snpmer[ref_h->k] = '\0';
				}
				cnt_ref = cnt2;
				cnt_alt = cnt1;
			}
			
			// build GT field
			if (cnt_ref > 0 && cnt_alt == 0) gt = "0";
				else if (cnt_ref == 0 && cnt_alt > 0) gt = "1";
				else gt = ".";

			// build INFO field: seq_name,strand,pos;seq_name,strand,pos;...
			char info[65536];
			int info_len = 0;
			char info_strand;

			if (write_info) {
				khint_t info_k = pg_im_get(m->m, key);
				if (info_k != kh_end(m->m)) {
					kinfo_t *occ = &kh_val(m->m, info_k);
					for (int j = 0; j < occ->n; ++j) {
						const char *name = h->cnames.names[occ->i[j].seq_idx];
						khint_t ref_seen_k = strset_get(seen, name);
						if (ref_seen_k != kh_end(seen) && occ->n == 1) { // reference genome, one occurrence
							info[0] = '.'; info[1] = '\0'; info_len = 1;
						} else if (ref_seen_k != kh_end(seen) && occ->n > 1) { // reference genome, multiple occurrences
							if (j) {
								if (ref_strand) { // already on forward strand
									info_strand = m_val_strand(occ->i[j].postrand) ? '+' : '-';
								} else { // on reverse strand, change it for VCF standard
									info_strand = !m_val_strand(occ->i[j].postrand) ? '+' : '-';
								}
								info_len += snprintf(info + info_len, sizeof(info) - info_len,
													"%s|%c|%d,",
													name, info_strand, m_val_pos(occ->i[j].postrand) + 1); // VCF is 1-based
							}
						} else if (ref_seen_k == kh_end(seen) && occ->n == 1) { // not reference genome, one occurrence
							if (ref_strand) { // already on forward strand
								info_strand = m_val_strand(occ->i[j].postrand) ? '+' : '-';
							} else { // on reverse strand, change it for VCF standard
								info_strand = !m_val_strand(occ->i[j].postrand) ? '+' : '-';
							}
							info_len += snprintf(info + info_len, sizeof(info) - info_len,
												"%s|%c|%d,",
												name, info_strand, m_val_pos(occ->i[j].postrand) + 1); // VCF is 1-based
						} else { // not reference genome, multiple occurrences
							if (ref_strand) { // already on forward strand
								info_strand = m_val_strand(occ->i[j].postrand) ? '+' : '-';
							} else { // on reverse strand, change it for VCF standard
								info_strand = !m_val_strand(occ->i[j].postrand) ? '+' : '-';
							}
							info_len += snprintf(info + info_len, sizeof(info) - info_len,
												"%s|%c|%d,",
												name, info_strand, m_val_pos(occ->i[j].postrand) + 1); // VCF is 1-based
						}
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

			// append SNPMER= to INFO
			{
				char tmp[65600];
				if (strcmp(info, ".") == 0)
					snprintf(tmp, sizeof(tmp), "SNPMER=%s", snpmer);
				else
					snprintf(tmp, sizeof(tmp), "%s;SNPMER=%s", info, snpmer);
				memcpy(info, tmp, strlen(tmp) + 1);
			}

            fprintf(fp,
					"%s\t%d\t.\t%c\t%c\t.\t.\t%s\tGT:KC\t%s:%u,%u\n",
					ref_chrom_name, ref_pos, ref, alt,
					info, gt, cnt_ref, cnt_alt);
		}
    }

	strset_destroy(seen);
    fclose(fp);
}


static void count_alleles(const char *gt, int *ac, int n_alt, int *an)
{
    const char *p = gt;
    while (*p) {
        if (*p == '.') {
            p++;
        } else if (isdigit((unsigned char)*p)) {
            int v = 0;
            while (isdigit((unsigned char)*p)) { v = v * 10 + (*p - '0'); p++; }
            (*an)++;
            if (v >= 1 && v <= n_alt) ac[v - 1]++;
        } else {
            p++;
            continue;
        }
		
        if (*p == '/' || *p == '|') p++;
        else break;
    }
}

void merge_vcfs(const char *out_fn, const char *tmpdir, int n_fns, int n_snps, int ref_idx)
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
        return;

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

    // INFO header lines for the fields we compute during merge
    fputs("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele\">\n", out);
    fputs("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n", out);
    fputs("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency, for each ALT allele\">\n", out);

    fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (int i = 0; i < n_fns; ++i)
        fprintf(out, "\t%s", sample_names[i] ? sample_names[i] : "UNKNOWN");
    fprintf(out, "\n");

    // VCF content
    char **lines = malloc(n_fns * sizeof(char*));
    for (int i = 0; i < n_fns; ++i)
        lines[i] = malloc(65536);

    char **sample_gt = malloc(n_fns * sizeof(char*));

    for (int s = 0; s < n_snps; ++s) {
        for (int i = 0; i < n_fns; ++i) {
            if (!fps[i] || !fgets(lines[i], 65536, fps[i]))
                lines[i][0] = '\0';
            lines[i][strcspn(lines[i], "\n")] = '\0';
        }

        // merge INFO (KPOS) from all files, and grab SNPMER from the first genome only (same SNP-mer everywhere)
        char merged_info[65536];
        int info_len = 0;
        char snpmer[64] = {0};
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
                // walk ';'-separated subfields within this INFO column
                char *p = info_start;
                while (p < info_end) {
                    char *sub_end = memchr(p, ';', info_end - p);
                    if (!sub_end || sub_end > info_end) sub_end = info_end;

                    if (i == 0 && strncmp(p, "SNPMER=", 7) == 0) {
                        int len = sub_end - (p + 7);
                        if (len > 0 && len < (int)sizeof(snpmer)) {
                            memcpy(snpmer, p + 7, len);
                            snpmer[len] = '\0';
                        }
                    } else if (strncmp(p, "KPOS=", 5) == 0) {
                        char *kp = p + 5;
                        int len = sub_end - kp;
                        if (!(len == 1 && *kp == '.')) {
                            if (info_len > 0) merged_info[info_len++] = ',';
                            memcpy(merged_info + info_len, kp, len);
                            info_len += len;
                        }
                    }
                    p = sub_end + 1;
                }
            }
        }
        if (info_len == 0) { merged_info[0] = '.'; merged_info[1] = '\0'; }
        else merged_info[info_len] = '\0';

        // extract each sample's GT pointer before lines[0] is destroyed by fields parsing
        for (int i = 0; i < n_fns; ++i) {
            sample_gt[i] = NULL;
            int tc = 0;
            for (char *c = lines[i]; *c; ++c) {
                if (*c == '\t' && ++tc == 9) {
                    sample_gt[i] = c + 1;
                    break;
                }
            }
        }

        // parse fixed fields from lines[0] (modifies lines[0] in place)
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

        // number of ALT alleles, from the ALT field
        int n_alt = 1;
        if (fields[4])
            for (char *c = fields[4]; *c; ++c)
                if (*c == ',') ++n_alt;

        // compute AC / AN / AF across all samples' genotypes
        int *ac = calloc(n_alt, sizeof(int));
        int an = 0;
        for (int i = 0; i < n_fns; ++i)
            if (sample_gt[i])
                count_alleles(sample_gt[i], ac, n_alt, &an);

        char ac_str[4096] = {0}, af_str[4096] = {0};
        int ac_len = 0, af_len = 0;
        for (int a = 0; a < n_alt; ++a) {
            double af = an > 0 ? (double)ac[a] / an : 0.0;
            ac_len += snprintf(ac_str + ac_len, sizeof(ac_str) - ac_len, "%s%d", a ? "," : "", ac[a]);
            af_len += snprintf(af_str + af_len, sizeof(af_str) - af_len, "%s%.4g", a ? "," : "", af);
        }
        free(ac);

        // write CHROM..FILTER, INFO (AC, AN, AF, SNPMER, then KPOS if present), FORMAT
        fprintf(out, "%s\t%s\t%s\t%s\t%s\t%s\t%s\tAC=%s;AN=%d;AF=%s",
                fields[0], fields[1], fields[2], fields[3],
                fields[4], fields[5], fields[6],
                ac_str, an, af_str);
        if (snpmer[0])
            fprintf(out, ";SNPMER=%s", snpmer);
        if (merged_info[0] != '.')
            fprintf(out, ";KPOS=%s", merged_info);
        fprintf(out, "\t%s", fields[8]);

        // sample columns, using saved GT pointers
        for (int i = 0; i < n_fns; ++i)
            fprintf(out, "\t%s", sample_gt[i] ? sample_gt[i] : ".");
        fprintf(out, "\n");
    }

    free(sample_gt);
    for (int i = 0; i < n_fns; ++i) {
        if (fps[i]) fclose(fps[i]);
        free(lines[i]);
        free(sample_names[i]);
    }
    free(fps); free(lines); free(sample_names);

    if (out && out != stdout)
        fclose(out);
}