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
	int n, m;
	int n_ins;
	ch_seq_t *a;
} ch_buf_t;

typedef struct { 			// global data shared across all file threads
	const pg_opt_t *opt; 	// terminal options
    const char **fns;
	pg_mht_t *h;
	pthread_mutex_t mutex;
	pthread_rwlock_t rwlock;
    int n_fns;
	int next;          		// index of next genome to process
	int done;				// number of genomes fully processed
	int filt;				// flag to indicate if filtering has been done
    int n_batch;            // number of parallel genome slots
    int *batch_threads;		// number of threads assigned to each file
} filedat_t;

typedef struct { 			// per-file data for each kt_pipeline() instance
    filedat_t *f; 			// points back to shared global data
    kseq_t *ks; 			// file-specific sequence reader
    const char *fn;
	int f_threads; 			// number of threads for processing fn
	int gnm_id;
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
	pg_mht_t *h = s->p->f->h;

	if (s->p->f->filt)
		pg_mht_insert_list(h, b->n, b->a, s->p->gnm_id);
	else
		b->n_ins += pg_mht_insert_list(h, b->n, b->a, 0);
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
			if (l < p->f->opt->k) continue;
			if (s->n == s->m) {
				s->m = s->m < 16? 16 : s->m + (s->n>>1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
			}
			MALLOC(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);
			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += l - p->f->opt->k + 1;
			if (s->sum_len >= p->f->opt->chunk_size)
				break;
		}
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->f->opt->pre, m;
		CALLOC(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i) {
			s->buf[i].m = m;
			MALLOC(s->buf[i].a, m);
		}
		for (i = 0; i < s->n; ++i) {
			count_seq_buf(s->buf, p->f->opt->k, p->f->opt->pre, s->len[i], s->seq[i]);
			free(s->seq[i]);
		}
		free(s->seq); free(s->len);
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->f->opt->pre;
		uint64_t n_ins = 0;
		kt_for(p->f_threads-2, worker_for, s, n);
		for (i = 0; i < n; ++i) {
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
		}
		__sync_fetch_and_add(&p->f->h->n_ins_tot, n_ins); // multiple I/O threads that update the number of k-mers

		free(s->buf); free(s);
	}
	return 0;
}


static void *file_worker(void *data)
{
    pldat_t *pl = (pldat_t*)data;
	filedat_t *fd = pl->f; // alias for convenience (not to use pl->f every time)

	while (1) { // cannot use n_fns because different threads are concurrently updating the number of processed genomes

		pthread_mutex_lock(&fd->mutex);
		// grab next genome index
        int i = fd->next++;
		if (i >= fd->n_fns) {
			pthread_mutex_unlock(&fd->mutex);
			break;  // all genomes have been processed, exit thread
		}
        pl->fn = fd->fns[i];
		pl->gnm_id = i + 1; // assign a genome ID and keep it fixed (next and done in filedat_t mutate)

		if (fd->opt->verbose)
			fprintf(stderr, "[M::%s] Processing '%s'\n", __func__, pl->fn);

		pthread_mutex_unlock(&fd->mutex);


		pthread_rwlock_rdlock(&fd->rwlock);
        gzFile fp = gzopen(pl->fn, "r");
        if (fp == 0) {
            fprintf(stderr, "[E::%s] failed to open '%s'\n", __func__, pl->fn);
            continue;
        }
        pl->ks = kseq_init(fp);
        kt_pipeline(3, worker_pipeline, pl, 3);
        kseq_destroy(pl->ks);
        gzclose(fp);
		pthread_rwlock_unlock(&fd->rwlock);

		pthread_mutex_lock(&fd->mutex);
		fd->done++;
		if (fd->opt->verbose)
			fprintf(stderr, "[M::%s] %d genomes done, %ld distinct k-mers in the hash table\n", __func__, fd->done, (long)fd->h->n_ins_tot);
		// int do_filt = 0;
		// if (fd->opt->min_freq < 1.0) { // if min_freq is 1, only filter at the end
		// 	int t = (1.0 - fd->opt->min_freq) * fd->n_fns;
		// 	do_filt = (fd->done > t && !fd->filt);
		// }
		pthread_mutex_unlock(&fd->mutex);

		// // filter low frequency k-mers every n files processed (based on min_freq)
		// if (do_filt) {
		// 	pthread_rwlock_wrlock(&fd->rwlock);
		// 	if (!fd->filt) {
		// 		fprintf(stderr, "[M::%s] Filtering k-mers with frequency < %.2f (after processing %d files)\n", __func__, fd->opt->min_freq, fd->done);
		// 		int n_del = pg_mht_filter(fd->h, fd->done, fd->n_fns, fd->opt->min_freq);
		// 		fprintf(stderr, "[M::%s] Filtered %d k-mer entries\n", __func__, n_del);
		// 		fd->filt = 1; // set the flag to indicate that filtering has been done
		// 	}
		// 	pthread_rwlock_unlock(&fd->rwlock);
		// }
	}
	
    pthread_exit(0);
}


pg_mht_t *pg_count(const char **fns, int n_fns, const pg_opt_t *opt, int filt, pg_mht_t *h_init)
{
    filedat_t fd;
    fd.fns = fns;
    fd.n_fns = n_fns;
    fd.opt = opt;
	if (h_init)
		fd.h = h_init;
	else
		fd.h = pg_mht_init(opt->k, opt->pre);
	fd.next = 0;
	fd.done = 0;
	fd.filt = filt;
	// split threads among input files as evenly as possible
	fd.batch_threads = (int*)calloc(n_fns, sizeof(int));
	fd.n_batch = assign_threads(opt->n_threads, n_fns, fd.batch_threads);

	pthread_mutex_init(&fd.mutex, 0);
	pthread_rwlock_init(&fd.rwlock, NULL);

	pthread_t *tid = (pthread_t*)calloc(fd.n_batch, sizeof(pthread_t));
	pldat_t   *pl  = (pldat_t*)calloc(fd.n_batch, sizeof(pldat_t));
	for (int i = 0; i < fd.n_batch; ++i) {
		pl[i].f = &fd;
		pl[i].f_threads = fd.batch_threads[i];
		pthread_create(&tid[i], 0, file_worker, &pl[i]);
	}
	for (int i = 0; i < fd.n_batch; ++i) pthread_join(tid[i], 0);

    free(tid); free(pl); free(fd.batch_threads);

	pthread_mutex_destroy(&fd.mutex);
	pthread_rwlock_destroy(&fd.rwlock);

    return fd.h;
}