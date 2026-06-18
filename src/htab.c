#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

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
	pg_ht1_t *g;
	if (n == 0) return 0;

	uint32_t gnm_cnt1, gnm_cnt2, cb1, cb2, snp1, snp2, pgnm_cnt1, pgnm_cnt2, filter, v;

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
			fprintf(fp, "%s\n", seq);

			seq[mid] = nt4_seq_table[cb2];
			fprintf(fp, "%s\n", seq);

		}
	}	
}


// VCF writer
void write_vcf(const char *out_fn, pg_mht_t *h, const paf_rec_t *recs, int n_snps, char *gnm_fn)
{
    FILE *fp = fopen(out_fn, "w");
    if (!fp) {
        fprintf(stderr, "[E::write_vcf] failed to open '%s'\n", out_fn);
        return;
    }

    // VCF header
    fprintf(fp,
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
		"##FORMAT=<ID=KC,Number=2,Type=Integer,Description=\"K-mer counts for REF and ALT alleles\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",
        gnm_fn);

	// VCF content
    int snp_idx = 0; // tracks position in paf_rec array
    for (int i = 0; i < 1<<h->pre; ++i) {
        pg_ht1_t *g = &h->h[i];
        khint_t k;

        for (k = 0; k < kh_end(g->h); ++k) {
            if (!kh_exist(g->h, k)) continue;

            uint32_t v = kh_val(g->h, k);
            uint32_t cb1 = s_val_cb1(v);
            uint32_t cb2 = s_val_cb2(v);
            uint32_t cnt1 = s_val_count1(v);
            uint32_t cnt2 = s_val_count2(v);

            const paf_rec_t *rec = &recs[snp_idx++];

            char ref = nt4_seq_table[cb1];
            char alt = nt4_seq_table[cb2];

            // Simple GT: if cnt2 > 0 and cnt1 > 0 → het (0/1), else hom
            const char *gt;
            if (cnt1 == 0 && cnt2 > 0) gt = "0/1";
            else if (cnt1 > 0 && cnt2 == 0) gt = "1/0";
			else if (cnt1 == 0 && cnt2 == 0) gt = "0/0";
            else gt = "1/1";

            fprintf(fp,
				"%s\t%ld\t.\t%c\t%c\t.\t.\t.\tGT:KC\t%s:%u/%u\n",
				rec->chrom_name,
				rec->target_pos,
				ref, alt,
				gt,
				cnt1, cnt2);
		}
    }

    fclose(fp);
}


void merge_vcfs(const char *out_fn, const char *tmpdir, int n_fns, int n_snps)
{
    FILE **fps = malloc(n_fns * sizeof(FILE*));
    for (int i = 0; i < n_fns; ++i) {
        char p[4096];
        snprintf(p, sizeof p, "%s/gnm.%d.vcf", tmpdir, i);
        fps[i] = fopen(p, "r");
        if (!fps[i]) {
            fprintf(stderr, "[E::merge_vcfs] failed to open '%s'\n", p);
            // handle error
        }
    }

    FILE *out = fopen(out_fn, "w");

    // VCF header
    char **sample_names = malloc(n_fns * sizeof(char*));
    char line[4096];
    for (int i = 0; i < n_fns; ++i) {
        sample_names[i] = NULL;
        while (fgets(line, sizeof line, fps[i])) {
            if (strncmp(line, "#CHROM", 6) == 0) {
                // extract sample name (last tab-separated token)
                char *tok = strrchr(line, '\t');
                if (tok) {
                    tok++; // skip the tab
                    // strip newline
                    tok[strcspn(tok, "\n")] = '\0';
                    sample_names[i] = strdup(tok);
                }
                break; // done with header for this file
            }
            // only write meta-lines from genome 0
            if (i == 0 && line[0] == '#') {
                fputs(line, out);
            }
        }
    }

    fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (int i = 0; i < n_fns; ++i)
        fprintf(out, "\t%s", sample_names[i] ? sample_names[i] : "UNKNOWN");
    fprintf(out, "\n");

	// VCF content
    char **lines = malloc(n_fns * sizeof(char*));
    for (int i = 0; i < n_fns; ++i)
        lines[i] = malloc(4096);

    for (int s = 0; s < n_snps; ++s) {
        // read one line from each file
        for (int i = 0; i < n_fns; ++i) {
            if (!fgets(lines[i], 4096, fps[i])) {
                fprintf(stderr, "[E::merge_vcfs] unexpected EOF in genome %d at SNP %d\n", i, s);
            }
            // strip newline
            lines[i][strcspn(lines[i], "\n")] = '\0';
        }

        // we want everything up to and including FORMAT for the first file
        char *p = lines[0];
        int tab_count = 0;
        char *sample_start = NULL;
        for (char *c = p; *c; ++c) {
            if (*c == '\t') {
                ++tab_count;
                if (tab_count == 9) { // after FORMAT field
                    *c = '\0';
                    sample_start = c + 1;
                    break;
                }
            }
        }
        fprintf(out, "%s\t%s", p, sample_start); // fixed fields + sample 0

        // append sample columns from following genomes
        for (int i = 1; i < n_fns; ++i) {
            // extract just the sample column (after the 9th tab)
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

    // cleanup
    for (int i = 0; i < n_fns; ++i) {
        fclose(fps[i]);
        free(lines[i]);
        free(sample_names[i]);
    }
    free(fps); free(lines); free(sample_names);
    fclose(out);
}