#include <stdio.h>
#include <stdlib.h>

#include "htab.h"
#include "ketopt.h"
#include "sys.h"
#include "utils.h"

int main_count(int argc, char *argv[])
{   
    pg_mht_t *h;
	char *ref_fn = 0;
	char *fn_out = 0;
	int c, i;
	pg_opt_t opt;
	ketopt_t o = KETOPT_INIT;
	pg_opt_init(&opt);

	static ko_longopt_t long_opts[] = {
        { "kmer", ko_required_argument, 301 },
        { "min_freq", ko_required_argument, 302 },
        { "pre", ko_required_argument, 303 },
        { "filt_type", ko_required_argument, 304 },
        { "chunk_size", ko_required_argument, 305 },
        { "threads", ko_required_argument, 306 },
		{ "write_info", ko_no_argument, 307 },
		{ "ref", ko_required_argument, 308 },
		{ "output", ko_required_argument, 309 },
        { "verbose", ko_no_argument, 310 },
        { 0, 0, 0 }
    };

	while ((c = ketopt(&o, argc, argv, 1, "k:m:p:f:K:t:wr:o:v", long_opts)) >= 0) {
        if      (c == 'k' || c == 301) opt.k = atoi(o.arg);
        else if (c == 'm' || c == 302) opt.min_freq = atof(o.arg);
        else if (c == 'p' || c == 303) opt.pre = atoi(o.arg);
        else if (c == 'f' || c == 304) opt.filt_type = atoi(o.arg);
        else if (c == 'K' || c == 305) opt.chunk_size = mm_parse_num(o.arg);
        else if (c == 't' || c == 306) opt.n_threads = atoi(o.arg);
		else if (c == 'w' || c == 307) opt.write_info = 1;
		else if (c == 'r' || c == 308) ref_fn = o.arg;
		else if (c == 'o' || c == 309) fn_out = o.arg;
        else if (c == 'v' || c == 310) opt.verbose = 1;
    }

	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: pgtools count [options] <in1.fa> [in2.fa [...]]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  %-10s k-mer size [%d]\n",            "-k INT",  opt.k);
		fprintf(stderr, "  %-10s minimum frequency of k-mers across inputs to be kept [%g]\n",
				"-m FLOAT", opt.min_freq);
		fprintf(stderr, "  %-10s prefix length [%d]\n",         "-p INT",  opt.pre);
		fprintf(stderr, "  %-10s filter type [%d]\n",           "-f INT",  opt.filt_type);
		fprintf(stderr, "  %-10s number of worker threads [%d]\n",
				"-t INT",  opt.n_threads);
		fprintf(stderr, "  %-10s chunk size [1.9g]\n",          "-K INT");
		fprintf(stderr, "  %-10s writing all hits in the pangenome in the info field\n",          "-K INT");
		fprintf(stderr, "  %-10s path of the reference genome\n","-r FILE");
		fprintf(stderr, "  %-10s verbose output\n",             "-v");
		fprintf(stderr, "  %-10s output genome-specific SNPs in VCF format\n",
				"-o FILE");
		return 1;
	}
	if (opt.k >= 32 || !(opt.k % 2)) {
		fprintf(stderr, "ERROR: -k must be odd and <=31\n");
		return 1;
	}

	// first step: count k-mers in the input files argv and filter for SNP-mers
	h = pg_count_k(argv + o.ind, argc - o.ind, &opt);

	pg_mht_tighten(h);

	// third step: count SNPmers in each file
	if (ref_fn == NULL) {
		ref_fn = argv[o.ind];
		fprintf(stderr, "[M::%s] No reference passed, taking %s as a reference\n", __func__, argv[o.ind]);
	}

	// reduce chunk_size for second pass
	opt.chunk_size = opt.chunk_size / 10;
	if (opt.chunk_size < 1024*1024) opt.chunk_size = 1024*1024; // minimum 1MB

	if (fn_out == NULL) fn_out = "-"; // redirect output to stdout
	pg_count_snp(argv + o.ind, argc - o.ind, h->n_ins_tot, &opt, h, ref_fn, fn_out);

	fprintf(stderr, "[M::%s] Analyzed %d files\n", __func__, argc - o.ind);
	
    return 0;
}


int main(int argc, char *argv[])
{   
	int ret = 1;
    pg_reset_realtime();
    if (strcmp(argv[1], "count") == 0) ret = main_count(argc-1, argv+1);

    if (ret == 0) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, PG_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (int i = 0; i < argc; ++i) {
			fprintf(stderr, " %s", argv[i]);
		}
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, pg_realtime(), pg_cputime(), pg_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}

	return ret;
}