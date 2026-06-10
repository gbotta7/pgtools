#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>

#include "kseq.h"
#include "kthread.h"
#include "htab.h"
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
	int64_t n_ins;
	ch_seq_t *a;
	int n, m;
} ch_buf_t;

typedef struct {
	pthread_mutex_t mutex;
	pg_id_map_t *id_maps;
	pg_csr_t *csr;
	const char **fns;
	int n_done;
	int n_fns;
} filedat_t;

typedef struct {
	filedat_t *f;
    const pg_opt_t *opt;
    kseq_t *ks; 			// file-specific sequence reader
	const char *fn;
    pg_mht_t *h;
	int filt;				// whether the intermediate k-mer filtering has been done for the first time
	int snp;				// whether it is the kmer pass or the snpmer pass
	int f_threads;
} pldat_t;

typedef struct { 			// data structure for each step in kt_pipeline()
    pldat_t *p;
	int *len;
	ch_buf_t *buf;
    char **seq;
    int n, m, sum_len, nk;
} stepdat_t;

static inline void ch_insert_buf(ch_buf_t *buf, int p, int k, uint64_t flanks, uint64_t center) // insert a k-mer $y to a linear buffer
{
	int pre = flanks & ((1<<p) - 1);
	ch_buf_t *b = &buf[pre];
	// fprintf(stderr, "[M::%s] b.m = %d, b.n = %d, flanks = %lu, center = %lu, pre = %d\n", __func__, b->m, b->n, (unsigned long)flanks, (unsigned long)center, pre);
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->a, b->m);
	}
	b->a[b->n].h_flanks = flanks;
    b->a[b->n].cb = center;
    b->n++;
}

static void count_seq_buf(ch_buf_t *buf, int k, int p, int len, const char *seq) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	uint64_t hash_mask = (1ULL<<((k-1)*2)) - 1; // to hash only the flanks
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				uint64_t center = (y >> ((k/2)*2)) & 3;           		// extract center from raw y
				uint64_t flanks = (y & ((1ULL<<(k/2)*2)-1))          	// right flank from raw y
								| ((y >> ((k/2+1)*2)) << ((k/2)*2)); 	// left flank from raw y
				ch_insert_buf(buf, p, k, pg_hash64(flanks, hash_mask), center);
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}


static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
	stepdat_t *s = (stepdat_t*)data;
	ch_buf_t *b = &s->buf[i];
	pg_mht_t *h = s->p->h;

	if (s->p->snp)
		pg_mht_count_list(h, b->n, b->a);
	else
		b->n_ins += pg_mht_insert_list(h, b->n, b->a, s->p->filt);
}

static void clear_for(void *data, long i, int tid) // callback for kt_for()
{
	pldat_t *p = (pldat_t*)data;
	if (p->snp)
		pg_mht_clear_s(p->h, i);
	else
		pg_mht_clear_k(p->h, i, p->opt->filt_type);
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
			int l = p->ks->seq.l;
			if (l < p->opt->k) continue;
			if (s->n == s->m) {
				s->m = s->m < 16? 16 : s->m + (s->n>>1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
			}
			MALLOC(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);
			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += l - p->opt->k + 1;
			if (s->sum_len >= p->opt->chunk_size)
				break;
		}
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->opt->pre, m;
		CALLOC(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i) {
			s->buf[i].m = m;
			MALLOC(s->buf[i].a, m);
		}
		for (i = 0; i < s->n; ++i) {
			count_seq_buf(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i]);
			free(s->seq[i]);
		}
		free(s->seq); free(s->len);
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->opt->pre;
		int n_ins = 0;
		int f_threads = p->snp ? p->f_threads : p->opt->n_threads-2;
		kt_for(f_threads, worker_for, s, n);
		for (i = 0; i < n; ++i) {
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
		}
		p->h->n_ins_tot += n_ins;

		free(s->buf); free(s);
	}
	return 0;
}

pg_mht_t *pg_count(const char **fns, const int n_fns, const pg_opt_t *opt)
{	
	pldat_t pl;
	pl.h = pg_mht_init(opt->k, opt->pre);
	pl.opt = opt;
	pl.filt = 0;
	pl.snp = 0; // kmer pass
	pl.h->n_del_tot = 0;
	const char *fn;

	for (int i = 0; i < n_fns; ++i) {
		fn = fns[i];
		if (opt->verbose)
			fprintf(stderr, "[M::%s] Processing '%s'\n", __func__, fn);

		gzFile fp;
		fp = fn == 0 || strcmp(fn, "-") == 0? gzdopen(0, "r") : gzopen(fn, "r");
		if (fp == 0) return 0;
		pl.ks = kseq_init(fp);
		pl.fn = fn;

		kt_pipeline(3, worker_pipeline, &pl, 3);

		kseq_destroy(pl.ks);
		gzclose(fp);

		// reset counters for the next round of counting (if any)
		kt_for(pl.opt->n_threads, clear_for, &pl, 1 << pl.opt->pre);

		// update number of processed files and filter k-mers if needed
		pl.h->n_done++;

		int check_fr = (int)round(n_fns * (1.0 - opt->min_freq)) + 1;
		if (!(pl.h->n_done % check_fr) && pl.h->n_done < n_fns) {
			if (opt->verbose) {
				fprintf(stderr, "[M::%s] Filtering k-mers\n", __func__);
			}
			int64_t n_del = pg_mht_filter(pl.h, pl.h->n_done, n_fns, pl.opt->min_freq, 0);
			pl.h->n_del_tot += n_del;
			pl.h->n_ins_tot -= n_del;
			pl.filt = 1;
			if (opt->verbose) {
				fprintf(stderr, "[M::%s] Filtered %ld k-mer entries\n", __func__, n_del);
			}
		}

		if (opt->verbose) {
			fprintf(stderr, "[M::%s] Processed %d genomes\n", __func__, pl.h->n_done);
			fprintf(stderr, "[M::%s] Current number of k-mer entries in the hash table: %ld\n", __func__, pl.h->n_ins_tot);
		}
	}

	if (n_fns > 1) {
		fprintf(stderr, "[M::%s] Final filtering to get SNP-mers\n", __func__);
		int64_t n_del = pg_mht_filter(pl.h, n_fns, n_fns, pl.opt->min_freq, 1); // filter to keep only SNP-mers
		pl.h->n_del_tot += n_del;
		pl.h->n_ins_tot -= n_del;
		if (opt->verbose) {
			fprintf(stderr, "[M::%s] Filtered %ld k-mer entries\n", __func__, n_del);
		}
	}

	pg_mht_tighten(pl.h);

    return pl.h;
}


static void rearrange_for(void *data, long i, int tid) // callback for kt_for()
{
	pg_mht_t *h = (pg_mht_t*)data;
	pg_mht_rearrange(h, i);
}


static void *worker_file(void *data)
{
    pldat_t *pl = (pldat_t*)data;
	filedat_t *fd = pl->f; // alias for convenience (not to use pl->f every time)

	while (1) { // cannot use n_fns because different threads are concurrently updating the number of processed genomes
        pthread_mutex_lock(&fd->mutex);
		// grab next genome index
        int i = fd->n_done++;
		if (i >= fd->n_fns) {
			pthread_mutex_unlock(&fd->mutex);;
			break;  // all genomes have been processed, exit thread
		}
		if (pl->opt->verbose) {
			fprintf(stderr, "[M::%s] Counted SNPs in %d genomes\n", __func__, fd->n_done);
		}
        pthread_mutex_unlock(&fd->mutex);

        pl->fn = fd->fns[i];
        gzFile fp = gzopen(pl->fn, "r");
        if (fp == 0) {
            fprintf(stderr, "[E::%s] failed to open '%s'\n", __func__, pl->fn);
            continue;
        }
        pl->ks = kseq_init(fp);
        kt_pipeline(3, worker_pipeline, pl, 3);
        kseq_destroy(pl->ks);
        gzclose(fp);

		// store counts to sparse matrix
		pthread_mutex_lock(&fd->mutex);
		pg_csr_insert(fd->csr, pl->h, fd->id_maps, i);
		pthread_mutex_unlock(&fd->mutex);

		kt_for(pl->f_threads, clear_for, pl, 1 << pl->opt->pre); // clear mht in this thread
    }
	
    pthread_exit(0);
}


pg_csr_t *pg_findsnp(const char **fns, const int n_fns, int n_snps, const pg_opt_t *opt, pg_mht_t *h)
{	
	// shift bits of the hash table values to count SNPs
	kt_for(opt->n_threads, rearrange_for, h, 1 << opt->pre);
	
	filedat_t fd;
	fd.id_maps = pg_mht_idx(h); // index hash table to later convert it to sparse matrix
	fd.csr = pg_csr_init(n_snps, n_fns, h, fd.id_maps); // init sparse matrix and store snpmers
	fd.fns = fns;
	fd.n_fns = n_fns;
	fd.n_done = 0;
	pthread_mutex_init(&fd.mutex, 0);

	// count SNPmers in each genome in parallel
	int *batch_threads = (int*)calloc(n_fns, sizeof(int));
	int n_batch = assign_threads(opt->n_threads, n_fns, batch_threads);

	pthread_t *tid = (pthread_t*)calloc(n_batch, sizeof(pthread_t));
	pldat_t *pl  = (pldat_t*)calloc(n_batch, sizeof(pldat_t));
	for (int i = 0; i < n_batch; ++i) {
		pl[i].f = &fd;
		pl[i].f_threads = batch_threads[i];
		pl[i].h = pg_mht_copy(h);
		pl[i].opt = opt;
		pl[i].snp = 1; // snpmer pass
		pthread_create(&tid[i], 0, worker_file, &pl[i]);
	}
	for (int i = 0; i < n_batch; ++i) {
		pthread_join(tid[i], 0);
		pg_mht_destroy(pl[i].h);
	}

    free(tid); free(pl); free(batch_threads);

	pthread_mutex_destroy(&fd.mutex);

	h = pl->h; // update for main

    return fd.csr;
}