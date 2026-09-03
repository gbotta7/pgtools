#include "khtab.h"
#include "utils.h"

pg_mkht_t *pg_mkht_init(int k, int pre, int w)
{
	pg_mkht_t *h;
	int i;
	CALLOC(h, 1);
	h->k = k;
	h->pre = pre;
	pthread_mutex_init(&h->mutex, 0);
	if (w) {
		CALLOC(h->ih, 1<<h->pre);
	} else {
		CALLOC(h->h, 1<<h->pre);
	}
	for (i = 0; i < 1<<h->pre; ++i) {
		if (w) {
			h->ih[i].ih = pg_kiht_init();
		} else {
			h->h[i].h = pg_kht_init();
		}
	}

	return h;
}

void pg_mkht_destroy(pg_mkht_t *h, int w)
{
	int i;
	if (h == 0) return;
	pthread_mutex_destroy(&h->mutex);
	for (i = 0; i < 1<<h->pre; ++i) {
		// destroy hash table for each bucket
		if (w) {
			pg_kiht1_destroy(&h->ih[i]);
		} else {
			pg_kht_destroy(h->h[i].h);
		}
	}
	free(h->h); free(h->ih); free(h);
}

void pg_kiht1_destroy(pg_kiht1_t *ih)
{
    khint_t k;
    for (k = 0; k < kh_end(ih->ih); ++k) {
        if (!kh_exist(ih->ih, k)) continue;
        kci_t *info = &kh_val(ih->ih, k);
        free(info->i);
    }
    pg_kiht_destroy(ih->ih);
}

int pg_mkht_insert_list(pg_mkht_t *h, int n, const seq_t *a, seq_info_t *b, int w)
{
	int j, n_ins = 0;
	uint64_t mask = (1<<h->pre) - 1;
	pg_kht1_t *g;
	pg_kiht1_t *ig;
	if (n == 0) return 0;
	if (w) {
		ig = &h->ih[a[0].h_seq & mask];
	} else {
		g = &h->h[a[0].h_seq & mask];
	}
	for (j = 0; j < n; ++j) {
		int ins = 1, absent;
		uint64_t key = a[j].h_seq >> h->pre, kv, cnt;
		if ((a[j].h_seq & mask) != (a[0].h_seq & mask)) continue; // skip if the partition bits are different (should not happen)
		khint_t k;
		kci_t *vi;
		
		if (w) {
			k = pg_kiht_put(ig->ih, key << K_COUNTER_BITS, &absent);
			if (absent) { // first occurrence, k-mer unknown
				++n_ins;
				kh_key(ig->ih, k) = key << K_COUNTER_BITS | 0x1U;
				// add info
				vi = &kh_val(ig->ih, k);
				vi->i = NULL;
				vi->n = 1;
				vi->m = 1;
				REALLOC(vi->i, vi->m);
				vi->i[0].pos = b[j].pos;
				vi->i[0].seq_idx = b[j].idx;
			} else {
				k = pg_kiht_get(ig->ih, key << K_COUNTER_BITS);
				if (k == kh_end(ig->ih)) continue;

				kv = kh_key(ig->ih, k);
				cnt = kv & K_COUNTER_MAX;

				if (cnt < K_COUNTER_MAX) {
					cnt++;
				}
				kh_key(ig->ih, k) = (kv & ~K_COUNTER_MAX) | cnt;

				// add info
				vi = &kh_val(ig->ih, k);
				if (vi->n == vi->m) {
					vi->m = vi->m < 2 ? vi->m + 1 : vi->m + (vi->m >> 1); // grow by 1.5x, but always by >=1
					REALLOC(vi->i, vi->m);
				}
				int idx = vi->n++;
				vi->i[idx].pos = b[j].pos;
				vi->i[idx].seq_idx = b[j].idx;
			}
		} else {
			k = pg_kht_put(g->h, key << K_COUNTER_BITS, &absent);
			if (absent) { // first occurrence, k-mer unknown
				++n_ins;
				++kh_key(g->h, k);
			} else {
				k = pg_kht_get(g->h, key << K_COUNTER_BITS);
				if (k != kh_end(g->h) && (kh_key(g->h, k) & K_COUNTER_MAX) < K_COUNTER_MAX) {
					++kh_key(g->h, k);
				}
			}
		}
	}

	return n_ins;
}

void pg_mkht_tighten(pg_mkht_t *h)
{
	int i;
	for (i = 0; i < 1<<h->pre; ++i) {
		pg_kht_t *g = h->h[i].h;
		if (kh_size(g) * 3 < kh_capacity(g))
			pg_kht_resize(g, kh_size(g) * 3);
	}
}


void write_kmer_tsv(const char *out_fn, pg_mkht_t *h, const char *gnm_fn, int w, int w_mko)
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
		fprintf(fp, "kmer\t%s\tpositions\n", sample_name);
	} else {
		fprintf(fp, "kmer\t%s\n", sample_name);
	}
	

	pg_kiht1_t *ig = NULL;
	pg_kht1_t *g = NULL;
	uint64_t hash_mask = (1ULL << (h->k * 2)) - 1, key, kmer, cnt;
	kci_t *vi;
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
				kmer = pg_hash64_inv((((key >> K_COUNTER_BITS) << h->pre)) | (uint64_t)i, hash_mask);
				vi = &kh_val(ig->ih, k);
				cnt = key & K_COUNTER_MAX;
			} else {
				key = kh_key(g->h, k);
				kmer = pg_hash64_inv((((key >> K_COUNTER_BITS) << h->pre)) | (uint64_t)i, hash_mask);
				cnt = key & K_COUNTER_MAX;
			}

			// write k-mer
			for (int j = 0; j < h->k; ++j)
				seq[j] = nt4_seq_table[(kmer >> ((h->k - 1 - j) * 2)) & 3];
			seq[h->k] = '\0';

			fprintf(fp, "%s\t%u", seq, cnt);

			if (w) {
				int occ = vi->n < w_mko ? vi->n : w_mko;

				fputc('\t', fp);
				for (int j = 0; j < occ; ++j) {
					uint32_t pos = vi->i[j].pos;
					uint16_t idx = vi->i[j].seq_idx;
					const char *cname = h->cnames.names[idx];

					if (j) fputc(',', fp);
					fprintf(fp, "%s:%u", cname, pos);
				}
			}
			fprintf(fp, "\n");
		}
	}

	fclose(fp);
}