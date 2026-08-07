#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "htab.h"
#include "ketopt.h"
#include "sys.h"
#include "utils.h"

int main_count(int argc, char *argv[])
{   
    pg_mht_t *h;
	char *ref_fn = 0;
	char *fn_out = 0;
	const char *bed_list_fn = 0;
	const char **bed_fns = 0;
	int n_bed = 0;
	int n_fa; 
	int c;
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
		{ "bed", ko_required_argument, 311 },
        { 0, 0, 0 }
    };

	while ((c = ketopt(&o, argc, argv, 1, "k:m:p:f:K:t:wr:o:vb:", long_opts)) >= 0) {
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
		else if (c == 'b' || c == 311) bed_list_fn = o.arg;
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
		fprintf(stderr, "  %-10s chunk size [%ld]\n",           "-K INT", (long)opt.chunk_size);
		fprintf(stderr, "  %-10s write all pangenome hits in the INFO field\n", "-w");
		fprintf(stderr, "  %-10s text file listing one BED path per line, one line per input FASTA, in the same order\n",
				"-b FILE");
		fprintf(stderr, "  %-10s path of the reference genome\n","-r FILE");
		fprintf(stderr, "  %-10s verbose output\n",             "-v");
		fprintf(stderr, "  %-10s output genome-specific SNPs in VCF format\n",
				"-o FILE");

		free(bed_fns);
		return 1;
	}
	if (opt.k >= 32 || !(opt.k % 2)) {
		fprintf(stderr, "ERROR: -k must be odd and <=31\n");

		free(bed_fns);
		return 1;
	}

	n_fa = argc - o.ind;

	// load bed filenames
	if (bed_list_fn) {
		CALLOC(bed_fns, n_fa); // bed files have to be as many as fasta filesd
		char buf[1024];
		FILE *fp = fopen(bed_list_fn, "r");
		if (fp == 0) {
			fprintf(stderr, "[E::%s] failed to open %s\n", __func__, bed_list_fn);
			free(bed_fns);
			return 1;
		}
		while (n_bed < n_fa && fgets(buf, sizeof(buf), fp)) {
			buf[strcspn(buf, "\r\n")] = 0; // strip the \n or \r 
			if (buf[0] == 0) continue;	   // check for empty lines
			bed_fns[n_bed++] = strdup(buf);
		}
		fclose(fp);
	}

	if (n_bed && n_bed != n_fa) {
		fprintf(stderr, "[E::%s] %d BED files given for %d FASTA files; -b must be given once per input or not at all\n", __func__, n_bed, n_fa);

		free(bed_fns);
		return 1;
	}
	
	// check if the passed reference is in the fasta files
	if (ref_fn) {
		int na_ref = 1;
		for (int i = 0; i < n_fa; ++i) {
			if (strcmp(argv[o.ind + i], ref_fn) == 0) {
				na_ref = 0;
				break;
			}
		}
		if (na_ref) {
			fprintf(stderr, "[E::%s] reference file %s not found in the list of inputs\n", __func__, ref_fn);
			free(bed_fns);
			return 1;
		}
	}

	// first step: count k-mers in the input files argv and filter for SNP-mers
	// h = pg_count_k(argv + o.ind, argc - o.ind, &opt);
	h = pg_count_k(argv + o.ind, n_bed ? bed_fns : 0, n_fa, &opt);

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
	// pg_count_snp(argv + o.ind, argc - o.ind, h->n_ins_tot, &opt, h, ref_fn, fn_out);
	pg_count_snp(argv + o.ind, n_bed ? bed_fns : 0, n_fa, h->n_ins_tot, &opt, h, ref_fn, fn_out);

	free(bed_fns);
	fprintf(stderr, "[M::%s] Analyzed %d files\n", __func__, n_fa);
	
    return 0;
}


int main(int argc, char *argv[])
{   
	int ret = 1;
    pg_reset_realtime();

	if (argc < 2) {
		fprintf(stderr, "Usage: pgtools <command> [options]\n");
		fprintf(stderr, "Commands:\n");
		fprintf(stderr, "  count    count k-mers and call genome-specific SNPs\n");
		return 1;
	}

    if (strcmp(argv[1], "count") == 0) ret = main_count(argc-1, argv+1);
	else {
		fprintf(stderr, "[E::main] unknown command '%s'\n", argv[1]);
		return 1;
	}

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