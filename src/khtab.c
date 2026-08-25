#include "khtab.h"

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
			h->ih[i].ih = pg_kht_init();
		}
	}

	return h;
}

void pg_mkht_destroy(pg_mkht_t *h)
{
	int i;
	if (h == 0) return;
	pthread_mutex_destroy(&h->mutex);
	for (i = 0; i < 1<<h->pre; ++i)
		pg_kht_destroy(h->h[i].h);
	free(h->h); free(h);
}

int pg_mkht_insert_list(pg_mkht_t *h, int n, const uint64_t *a)
{
	int j, n_ins = 0;
	uint64_t mask = (1<<h->pre) - 1;
	pg_kht1_t *g;
	if (n == 0) return 0;
	g = &h->h[a[0] & mask];
	for (j = 0; j < n; ++j) {
		int ins = 1, absent;
		uint64_t key = a[j] >> h->pre;
		if ((a[j] & mask) != (a[0] & mask)) continue; // skip if the partition bits are different (should not happen)
		khint_t k;
		
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

pg_mkht_t *pg_mkht_repopulate(const char *kmer_file, pg_opt_t *opt)
{	
	pg_mkht_t *h = pg_mkht_init(opt->k, opt->pre, opt->write_info);
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
		uint32_t snp1, cb1, cb2, v;
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
		pg_sht1_t *g = &h->h[bucket]; // get hash table partition for the k-mer
		key = h_flanks >> opt->pre;

		k = pg_sht_put(g->h, key, &absent);
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


void pg_dump_kmers(const char *fn, pg_msht_t *h, int snp)
{
	FILE *fp;
    if ((fp = strcmp(fn, "-") ? fopen(fn, "w") : stdout) == 0)
        return;

	uint64_t hash_mask = (1ULL << ((h->k - 1) * 2)) - 1;
	int mid = h->k >> 1;
	int right_off = snp ? mid + 5 : mid + 1;
	char seq[64];

	for (int i = 0; i < 1 << h->pre; ++i) {
		pg_sht1_t *g = &h->h[i];
		for (khint_t k = 0; k < kh_end(g->h); ++k) {
			if (!kh_exist(g->h, k)) continue;

			uint64_t flanks = pg_hash64_inv(((uint64_t)kh_key(g->h, k) << h->pre) | (uint64_t)i, hash_mask);
			uint64_t kv = kh_val(g->h, k);
			uint64_t cb1 = f_val_cb1(kv);
			uint64_t cb2 = f_val_cb2(kv);

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
				if (f_val_pgnm_count1(v) > 0 && !f_val_filt(v)) { // dump the first k-mer
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



void write_tsv(const char *out_fn, pg_msht_t *h, const char *gnm_fn, int snp)
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
		pg_sht1_t *g = &h->h[i];
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

				if (!s_val_filt(v)) {
					fprintf(fp, "%s\t%u\n", seq, s_val_count1(v));
				}
				if (!s_val_filt2(v)) {
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