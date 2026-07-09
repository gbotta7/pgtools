#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>
#include <unistd.h>

#include "kseq.h"
#include "kthread.h"
#include "htab.h"
#include "parser.h"
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
    int64_t n_ins;
    k_ch_seq_t *ak;  // first pass
    s_ch_seq_t *as;  // snp pass
    int n, m;
} ch_buf_t;

typedef struct {
	pthread_mutex_t mutex;
	const char *tmpdir;
	const char **fns;
	const char *ref;
	pg_mht_t *ref_h;		// master hash table for reference data
	int n_done;				// processed genomes
	int scan;				// idx of the scan when processing genomes (shifted from n_done due to ref processing)
	int n_fns;
	int n_snps;
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
	int *name_idx;			// contig names' index
    int n, m, sum_len, nk;
} stepdat_t;


static inline int ch_get_name_idx(cnames_t *nt, const char *name) {
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


static inline void ch_insert_buf(ch_buf_t *buf, pldat_t *p, uint64_t flanks, uint64_t center, uint32_t pos, int cname_idx, uint8_t strand) // insert a k-mer $y to a linear buffer
{	
	int pre = flanks & ((1<<p->opt->pre) - 1);
	ch_buf_t *b = &buf[pre];

	if (p->snp) {
		if (b->n == b->m) {
			b->m = b->m < 8? 8 : b->m + (b->m>>1);
			REALLOC(b->as, b->m);
		}
		b->as[b->n].h_flanks = flanks;
		b->as[b->n].cb = center;
		b->as[b->n].pos = pos;
		b->as[b->n].strand = strand;
		b->as[b->n].idx = cname_idx;
	} else {
		if (b->n == b->m) {
			b->m = b->m < 8? 8 : b->m + (b->m>>1);
			REALLOC(b->ak, b->m);
		}
		b->ak[b->n].h_flanks = flanks;
		b->ak[b->n].cb = center;
	}
	b->n++;
}

static void count_seq_buf(ch_buf_t *buf, pldat_t *p, int len, const char *seq, int cname_idx) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	uint64_t hash_mask = (1ULL<<((p->opt->k-1)*2)) - 1; // to hash only the flanks
	uint64_t x[2], mask = (1ULL<<p->opt->k*2) - 1, shift = (p->opt->k - 1) * 2;

	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= p->opt->k) { // we find a k-mer
				uint64_t y = x[0] < x[1] ? x[0] : x[1];
				uint64_t y_rev = x[0] < x[1] ? x[1] : x[0];
				uint64_t center = (y >> ((p->opt->k/2)*2)) & 3;           					// extract center from raw y
				// uint64_t cb = x[0] < x[1] ? center : 3 - center;
				uint64_t flanks = (y & ((1ULL<<(p->opt->k/2)*2)-1))          				// right flank from raw y
								| ((y >> ((p->opt->k/2+1)*2)) << ((p->opt->k/2)*2)); 		// left flank from raw y
				uint64_t rev_flanks = (y_rev & ((1ULL<<(p->opt->k/2)*2)-1))          		// right flank from raw y
								| ((y_rev >> ((p->opt->k/2+1)*2)) << ((p->opt->k/2)*2)); 	// left flank from raw y
				uint8_t strand = x[0] < x[1] ? 1 : 0;

				if (flanks == rev_flanks) continue;
			
				ch_insert_buf(buf, p, pg_hash64(flanks, hash_mask), center, (uint32_t)i-p->opt->k/2, cname_idx, strand); // i-k/2 is the 0-based position of center
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
		pg_mht_count_list(h, b->n, b->as);
	else
		b->n_ins += pg_mht_insert_list(h, b->n, b->ak, s->p->filt);
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
				REALLOC(s->name_idx, s->m);
				// REALLOC(s->name, s->m);
			}
			MALLOC(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);
			s->name_idx[s->n] = p->snp ? ch_get_name_idx(&p->h->cnames, p->ks->name.s) : -1;
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
		m = p->snp ? (int)(p->h->n_ins_tot * 1.2 / n) + 1 : (int)(s->nk * 1.2 / n) + 1; // pre-allocate much less memory in SNP pass
		for (i = 0; i < n; ++i) {
			s->buf[i].m = m;
			if (p->snp)
				CALLOC(s->buf[i].as, m);
			else
				CALLOC(s->buf[i].ak, m);
		}
		for (i = 0; i < s->n; ++i) {
			count_seq_buf(s->buf, p, s->len[i], s->seq[i],
              p->snp ? s->name_idx[i] : 0);
			free(s->seq[i]);
		}
		free(s->seq); free(s->len); free(s->name_idx);
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->opt->pre;
		int n_ins = 0;
		int f_threads = p->snp ? p->f_threads : p->opt->n_threads-2;
		kt_for(f_threads, worker_for, s, n);
		for (i = 0; i < n; ++i) {
			n_ins += s->buf[i].n_ins;
			if (p->snp)
				free(s->buf[i].as);
			else
				free(s->buf[i].ak);
		}
		p->h->n_ins_tot += n_ins;

		free(s->buf); free(s);
	}
	return 0;
}

pg_mht_t *pg_count_k(const char **fns, const int n_fns, const pg_opt_t *opt)
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
			pg_mht_tighten(pl.h);
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
		pg_mht_tighten(pl.h);
		if (opt->verbose) {
			fprintf(stderr, "[M::%s] Filtered %ld k-mer entries\n", __func__, n_del);
		}
	}

	pg_mht_tighten(pl.h);

	// if (out) {
	// 	pg_dump_snps(out, pl.h);
	// }

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
		int i = fd->scan++;
		if (i >= fd->n_fns) {
			pthread_mutex_unlock(&fd->mutex);
			break;
		}
		if (strcmp(fd->fns[i], fd->ref) == 0) {
			pthread_mutex_unlock(&fd->mutex);
			continue;
		}
		fd->n_done++;
		if (pl->opt->verbose) {
			fprintf(stderr, "[M::%s] Counting SNPs in genome number %d\n", __func__, fd->n_done);
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

		// store this genome's counts to its own file (no mutex: own file)
		char gnm_path[4096];
		snprintf(gnm_path, sizeof gnm_path, "%s/gnm.%d.vcf", fd->tmpdir, i);

		write_vcf(gnm_path, pl->h, fd->ref_h, pl->fn, pl->opt->write_info);

		kt_for(pl->f_threads, clear_for, pl, 1 << pl->opt->pre); // clear mht in this thread
    }
	
    pthread_exit(0);
}


void pg_count_snp(const char **fns, const int n_fns, int64_t n_snps, const pg_opt_t *opt, pg_mht_t *h, const char *ref_fn, const char *out_fn)
{	
	// shift bits of the hash table values to count SNPs
	kt_for(opt->n_threads, rearrange_for, h, 1 << opt->pre);

	filedat_t fd;
	fd.fns = fns;
	fd.ref = ref_fn;
	fd.n_snps = n_snps;
	fd.n_fns = n_fns;
	fd.n_done = 0;
	fd.scan = 0;
	pthread_mutex_init(&fd.mutex, 0);
	
	// create a temp directory for the per-genome count files
	char tmpdir[1024];
	snprintf(tmpdir, sizeof tmpdir, "%s.snptmp.XXXXXX", strcmp(out_fn, "-") ? out_fn : "pg");
	if (mkdtemp(tmpdir) == 0) {
		fprintf(stderr, "[E::%s] failed to create temp dir\n", __func__);
		return;
	}
	fd.tmpdir = tmpdir;

	// find the reference file
	int ref_idx = -1;
	if (fd.ref == NULL) {
		ref_idx = 0; // use the first input file as a reference
		fd.ref = fns[0];
	} else {
		for (int i = 0; i < fd.n_fns; ++i) {
			if (strcmp(fd.fns[i], fd.ref) == 0) {
				ref_idx = i;
				break;
			}
		}
	}
    if (ref_idx < 0) {
        fprintf(stderr, "[M::%s] Reference file %s not found in the list of filenames\n", __func__, fd.ref);
    }

	// count SNPmers in the reference first
	pldat_t pl_ref;
	pl_ref.f = &fd;
	pl_ref.f_threads = opt->n_threads;
	pl_ref.h = pg_mht_copy(h);
	pl_ref.fn = ref_fn;
	pl_ref.opt = opt;
	pl_ref.snp = 1; // snpmer pass
	// grab next genome index
	int i = fd.n_done++;
	if (pl_ref.opt->verbose) {
		fprintf(stderr, "[M::%s] Counting SNPs in genome number %d\n", __func__, fd.n_done);
	}

	gzFile fp = gzopen(ref_fn, "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open '%s'\n", __func__, ref_fn);
		return;
	}
	pl_ref.ks = kseq_init(fp);
	kt_pipeline(3, worker_pipeline, &pl_ref, 3);
	kseq_destroy(pl_ref.ks);
	gzclose(fp);

	// point fd to the reference mht
	fd.ref_h = pl_ref.h;

	// store this genome's counts to its own file (no mutex: own file)
	char gnm_path[4096];
	snprintf(gnm_path, sizeof gnm_path, "%s/gnm.%d.vcf", fd.tmpdir, ref_idx);
	write_vcf(gnm_path, pl_ref.h, fd.ref_h, pl_ref.fn, pl_ref.opt->write_info);
	

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

	// merge the per-genome files (opened in FD-safe batches) into the matrix
	if (opt->verbose)
		fprintf(stderr, "[M::%s] Merging %d per-genome files into '%s'\n", __func__, n_fns, out_fn);

	merge_vcfs(out_fn, tmpdir, n_fns, n_snps, ref_idx);

	// clean up per-genome VCFs
	for (int i = 0; i < n_fns; ++i) {
		char p[4096];
		snprintf(p, sizeof p, "%s/gnm.%d.vcf", tmpdir, i);
		remove(p);
	}
	remove(tmpdir);

	pg_mht_destroy(pl_ref.h);
	pg_mht_destroy(h);
}