#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>
#include <unistd.h>

#include "bed.h"
#include "kseq.h"
#include "kthread.h"
#include "shtab.h"
#include "parser.h"
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
    int64_t n_ins;
    seq_t *a;
    seq_info_t *i;  // snp pass
    int n, m;
} buf_t;

typedef struct {
    const pg_opt_t *opt;
    kseq_t *ks; 			// file-specific sequence reader
	bedmap1_t *b;
	const char *fa_fn;
	const char *bed_fn;
    pg_msht_t *h;
	int filt;				// whether the intermediate k-mer filtering has been done for the first time, do not insert new k-mers after it is set to 1
	int cnt;				// whether it is first or second pass
	int snp;				// whether you count SNP-mers or k-mers
	int n_done;				// processed genomes
	int n_fns;
} pldat_t;

typedef struct { 			// data structure for each step in kt_pipeline()
    pldat_t *p;
	int *len;
	buf_t *buf;
    char **seq;
	int *name_idx;			// contig names' index
    int n, m, sum_len, nk;
} stepdat_t;


static inline int get_name_idx(cnames_t *nt, const char *name) {
    if (nt->n > 0 && strcmp(nt->names[nt->n - 1], name) == 0)
        return nt->n - 1;

    if (nt->n == nt->m) {
        nt->m = nt->m < 8 ? 8 : nt->m + (nt->m >> 1);
        REALLOC(nt->names, nt->m);
    }
    MALLOC(nt->names[nt->n], strlen(name) + 1);
    memcpy(nt->names[nt->n], name, strlen(name) + 1);
    return nt->n++;
}


static inline void insert_buf(buf_t *buf, pldat_t *p, uint64_t flanks, uint64_t center, uint32_t pos, int cname_idx) // insert a k-mer $y to a linear buffer
{	
	int pre = flanks & ((1<<p->opt->pre) - 1);
	buf_t *b = &buf[pre];

	if (p->cnt && p->opt->write_info) {
		if (b->n == b->m) {
			b->m = b->m < 8? 8 : b->m + (b->m>>1);
			REALLOC(b->a, b->m);
			REALLOC(b->i, b->m);
		}
		b->a[b->n].h_flanks = flanks;
		b->a[b->n].cb = center;
		b->i[b->n].pos = pos;
		b->i[b->n].allele = center;
		b->i[b->n].idx = cname_idx;
	} else {
		if (b->n == b->m) {
			b->m = b->m < 8? 8 : b->m + (b->m>>1);
			REALLOC(b->a, b->m);
		}
		b->a[b->n].h_flanks = flanks;
		b->a[b->n].cb = center;
	}
	b->n++;
}

static void count_seq_buf(buf_t *buf, pldat_t *p, int len, const char *seq, int cname_idx) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	uint64_t hash_mask = (1ULL<<((p->opt->k-1)*2)) - 1; // to hash only the flanks
	uint64_t x[2], mask = (1ULL<<p->opt->k*2) - 1, shift = (p->opt->k - 1) * 2;

	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  								// forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  								// reverse strand
			if (++l >= p->opt->k) { // we find a k-mer
				uint64_t y = x[0] < x[1] ? x[0] : x[1];
				uint64_t y_rev = x[0] < x[1] ? x[1] : x[0];
				uint64_t center = (y >> ((p->opt->k/2)*2)) & 3;           					// extract center from raw y (already contains the right allele)
				uint64_t flanks = (y & ((1ULL<<(p->opt->k/2)*2)-1))          				// right flank from raw y
								| ((y >> ((p->opt->k/2+1)*2)) << ((p->opt->k/2)*2)); 		// left flank from raw y
				uint64_t rev_flanks = (y_rev & ((1ULL<<(p->opt->k/2)*2)-1))          		// right flank from raw y
								| ((y_rev >> ((p->opt->k/2+1)*2)) << ((p->opt->k/2)*2)); 	// left flank from raw y

				if (flanks == rev_flanks) continue;											// palindromic
			
				insert_buf(buf, p, pg_hash64(flanks, hash_mask), center, (uint32_t)i-p->opt->k/2, cname_idx); // i-k/2 is the 0-based position of center
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}


static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
	stepdat_t *s = (stepdat_t*)data;
	buf_t *b = &s->buf[i];
	pg_msht_t *h = s->p->h;

	if (s->p->cnt && s->p->opt->write_info)
		pg_msht_count_list(h, b->n, b->a, b->i);
	else if (s->p->cnt && !s->p->opt->write_info)
		pg_msht_count_list(h, b->n, b->a, 0);
	else
		b->n_ins += pg_msht_insert_list(h, b->n, b->a, s->p->filt);
}

static void clear_for(void *data, long i, int tid) // callback for kt_for()
{
	pldat_t *p = (pldat_t*)data;
	if (p->cnt)
		pg_msht_clear2(p->h, i, p->opt->write_info);
	else
		pg_msht_clear1(p->h, i, p->opt->filt_type, p->opt->mko);
}

static void filter_for(void *data, long i, int tid) // callback for kt_for()
{
	pldat_t *p = (pldat_t*)data;
	int ff = p->n_done == p->n_fns; 
	int64_t n_del = pg_msht_filter(p->h, i, p->n_done, p->n_fns, ff, p->opt);

	pthread_mutex_lock(&p->h->mutex);
	p->h->n_del_tot += n_del;
	p->h->n_ins_tot -= n_del;
	pthread_mutex_unlock(&p->h->mutex);

	p->filt = 1;
}

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
	pldat_t *p = (pldat_t*)data;	

	if (step == 0) { // step 1: read a block of sequences
		int ret;
		stepdat_t *s;
		CALLOC(s, 1);
		s->p = p;
		while ((ret = kseq_read(p->ks)) >= 0) {
			int32_t l = p->ks->seq.l;
			const bed_ctg_t *c = 0;
			int64_t nk_ctg, len_ctg;
			if (l < p->opt->k) continue;
			
			// load bed file if passed
			if (p->b) {
				c = bed_get(p->b, p->ks->name.s);	// get contig bed entries
				if (c == 0) {
					// if (p->opt->verbose) fprintf(stderr, "[E::%s] contig %s not found in BED file %s\n", __func__, p->ks->name.s, p->bed_fn);
					continue;
				};
				nk_ctg = bed_nk(c, l, p->opt->k); // compute number of possible k-mers
				if (nk_ctg == 0) {
					// if (p->opt->verbose) fprintf(stderr, "[E::%s] contig %s has no regions longer than %d bases in the BED file %s\n", __func__, p->ks->name.s, p->opt->k, p->bed_fn);
					continue;
				};
			} else {
				nk_ctg  = l - p->opt->k + 1;
			}

			if (s->n == s->m) {
				s->m = s->m < 16? 16 : s->m + (s->n>>1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
				REALLOC(s->name_idx, s->m);
			}
			MALLOC(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);

			// mask the sequence outside the BED entries if a BED file is provided
			if (c) mask_fa(s->seq[s->n], l, c);

			s->name_idx[s->n] = p->cnt ? get_name_idx(&p->h->cnames, p->ks->name.s) : -1;
			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += nk_ctg;
			if (s->sum_len >= p->opt->chunk_size)
				break;
		}
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1 << p->opt->pre, m;
		CALLOC(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i) {
			s->buf[i].m = m;
			CALLOC(s->buf[i].a, m);
			if (p->cnt && p->opt->write_info) {
				CALLOC(s->buf[i].i, m);
			}
		}
		for (i = 0; i < s->n; ++i) {
			count_seq_buf(s->buf, p, s->len[i], s->seq[i], s->name_idx[i]); // CHANGED CONDITION
			free(s->seq[i]);
		}
		free(s->seq); free(s->len); free(s->name_idx);
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->opt->pre;
		int n_ins = 0;
		kt_for(p->opt->n_threads-2, worker_for, s, n);
		for (i = 0; i < n; ++i) {
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
			if (p->cnt && p->opt->write_info) {
				free(s->buf[i].i);
			}				
		}
		p->h->n_ins_tot += n_ins;

		free(s->buf); free(s);
	}
	return 0;
}

pg_msht_t *pg_detect(const char **fa_fns, const char **bed_fns, const int n_fns, const pg_opt_t *opt, const char *out_fn)
{	
	pldat_t pl;
	pl.n_fns = n_fns;
	pl.n_done = 0;
	pl.h = pg_msht_init(opt->k, opt->pre, opt->write_info);
	pl.opt = opt;
	pl.filt = 0;
	pl.cnt = 0; // first pass
	pl.snp = opt->snp;
	pl.h->n_del_tot = 0;
	const char *fa_fn;
	const char *bed_fn;
	for (int i = 0; i < n_fns; ++i) {
		fa_fn = fa_fns[i];
		bed_fn = bed_fns ? bed_fns[i] : 0;
		if (!bed_fn && opt->verbose)
			fprintf(stderr, "[M::%s] Processing '%s'\n", __func__, fa_fn);
		if (bed_fn && opt->verbose)
			fprintf(stderr, "[M::%s] Processing '%s' in the regions specified by %s\n", __func__, fa_fn, bed_fn);

		// open fasta file
		gzFile fp;
		fp = fa_fn == 0 || strcmp(fa_fn, "-") == 0? gzdopen(0, "r") : gzopen(fa_fn, "r");
		if (fp == 0) return 0;
		pl.ks = kseq_init(fp);
		pl.fa_fn = fa_fn;

		// open bed file if passed
		pl.bed_fn = bed_fn;
		pl.b = 0;
		if (pl.bed_fn) {
			pl.b = bed_read(pl.bed_fn);
		}

		kt_pipeline(3, worker_pipeline, &pl, 3);

		kseq_destroy(pl.ks);
		gzclose(fp);

		if (pl.b) {
			bed_destroy(pl.b); pl.b = 0;
		}

		// reset counters for the next round of counting (if any)
		kt_for(pl.opt->n_threads, clear_for, &pl, 1 << pl.opt->pre);

		// update number of processed files and filter k-mers if needed
		pl.n_done++;

		int check_fr = (int)round(pl.n_fns * (1.0 - opt->msf)) + 1;
		if (!(pl.n_done % check_fr) && pl.n_done < n_fns) {
			if (opt->verbose) {
				fprintf(stderr, "[M::%s] Filtering k-mers\n", __func__);
			}
			kt_for(pl.opt->n_threads, filter_for, &pl, 1 << pl.opt->pre);
			pg_msht_tighten(pl.h);
			if (opt->verbose) {
				fprintf(stderr, "[M::%s] Current number of SNP-mer entries in the hash table: %ld\n", __func__, pl.h->n_ins_tot);
			}
		}

		if (opt->verbose) {
			fprintf(stderr, "[M::%s] Processed %d genomes\n", __func__, pl.n_done);
		}
	}

	if (n_fns > 1) {
		fprintf(stderr, "[M::%s] Final filtering to get specific k-mers\n", __func__);
		kt_for(pl.opt->n_threads, filter_for, &pl, 1 << pl.opt->pre); // final filter
		pg_msht_tighten(pl.h);
		if (opt->verbose) {
			fprintf(stderr, "[M::%s] Current number of SNP-mer entries in the hash table: %ld\n", __func__, pl.h->n_ins_tot);
		}
	}

	pg_msht_tighten(pl.h);

	if (out_fn) {
		if (pl.snp)
			pg_dump_snpmers(out_fn, pl.h);
	}
	
    return pl.h;
}


static void rearrange_for(void *data, long i, int tid) // callback for kt_for()
{
	pg_msht_t *h = (pg_msht_t*)data;
	pg_msht_rearrange(h, i);
}

void pg_count(const char *fa_fn, const char *bed_fn, const pg_opt_t *opt, pg_msht_t *h, const char *out_fn)
{
	pldat_t pl;
	pl.n_done = 0;
	pl.h = h;
	pl.opt = opt;
	pl.filt = 0;
	pl.cnt = 1;
	pl.snp = opt->snp;

	if (!bed_fn && opt->verbose)
		fprintf(stderr, "[M::%s] Processing '%s'\n", __func__, fa_fn);
	if (bed_fn && opt->verbose)
		fprintf(stderr, "[M::%s] Processing '%s' in the regions specified by %s\n", __func__, fa_fn, bed_fn);

	// open fasta file
	gzFile fp;
	fp = fa_fn == 0 || strcmp(fa_fn, "-") == 0? gzdopen(0, "r") : gzopen(fa_fn, "r");
	if (fp == 0) return 0;
	pl.ks = kseq_init(fp);
	pl.fa_fn = fa_fn;

	// open bed file if passed
	pl.bed_fn = bed_fn;
	pl.b = 0;
	if (pl.bed_fn) {
		pl.b = bed_read(pl.bed_fn);
	}

	kt_pipeline(3, worker_pipeline, &pl, 3);

	kseq_destroy(pl.ks);
	gzclose(fp);

	if (pl.b) {
		bed_destroy(pl.b); pl.b = 0;
	}

	if (out_fn) {
		if (pl.snp)
			write_snpmer_tsv(out_fn, pl.h, pl.fa_fn, pl.opt->write_info, pl.opt->write_mko);
	}
	
	pg_msht_destroy(pl.h, opt->write_info);
}