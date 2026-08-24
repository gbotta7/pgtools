#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shtab.h"
#include "ketopt.h"
#include "sys.h"
#include "utils.h"

int main_detect(int argc, char *argv[])
{
	pg_msht_t *h;
	char *fn_out = 0;
	char *bed_list_fn = 0;
	char **bed_fns = 0;
	int n_bed = 0;
	int n_fa; 
	int c;
	pg_opt_t opt;
	ketopt_t o = KETOPT_INIT;
	pg_opt_init(&opt);

	static ko_longopt_t long_opts[] = {
        { "k_length", ko_required_argument, 301 },
        { "msf", ko_required_argument, 302 },
		{ "maf", ko_required_argument, 303 },
		{ "snp", ko_no_argument, 304 },
		{ "mko", ko_required_argument, 305 },
        { "pre", ko_required_argument, 306 },
        { "filt_type", ko_required_argument, 307 },
        { "chunk_size", ko_required_argument, 308 },
        { "threads", ko_required_argument, 309 },
		{ "output", ko_required_argument, 310 },
        { "verbose", ko_no_argument, 311 },
		{ "bed", ko_required_argument, 312 },
        { 0, 0, 0 }
    };

	while ((c = ketopt(&o, argc, argv, 1, "k:p:f:K:t:wo:vb:", long_opts)) >= 0) {
        if      (c == 'k' || c == 301) opt.k = atoi(o.arg);
        else if (c == 302) opt.msf = atof(o.arg);
		else if (c == 303) opt.maf = atof(o.arg);
		else if (c == 304) opt.snp = 1;
		else if (c == 305) opt.mko = atoi(o.arg);
        else if (c == 'p' || c == 306) opt.pre = atoi(o.arg);
        else if (c == 'f' || c == 307) opt.filt_type = atoi(o.arg);
        else if (c == 'K' || c == 308) opt.chunk_size = mm_parse_num(o.arg);
        else if (c == 't' || c == 309) opt.n_threads = atoi(o.arg);
		else if (c == 'o' || c == 310) fn_out = o.arg;
        else if (c == 'v' || c == 311) opt.verbose = 1;
		else if (c == 'b' || c == 312) bed_list_fn = o.arg;
    }

	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: pgtools detect [options] <in1.fa> [in2.fa [...]]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  %-10s k-mer size [%d]\n",            "-k INT",  opt.k);
		fprintf(stderr, "  %-10s minimum sample frequency of k-mers across inputs to be kept [%g]\n",
				"-t FLOAT", opt.msf);
		fprintf(stderr, "  %-10s minimum allelic frequency of SNP-mers across input samples [%g]\n",
				"-t FLOAT", opt.maf);
		fprintf(stderr, "  %-10s snp mode to count SNP-mers\n",         "--snp");
		fprintf(stderr, "  %-10s maximum k-mer occurrences [%d]\n",         "--mko INT",  opt.mko);
		fprintf(stderr, "  %-10s prefix length [%d]\n",         "-p INT",  opt.pre);
		fprintf(stderr, "  %-10s filter type [%d]\n",           "-f INT",  opt.filt_type);
		fprintf(stderr, "  %-10s number of worker threads [%d]\n",
				"-t INT",  opt.n_threads);
		fprintf(stderr, "  %-10s chunk size [%ld]\n",           "-K INT", (long)opt.chunk_size);
		fprintf(stderr, "  %-10s text file listing one BED path per line, one line per input FASTA, in the same order\n",
				"-b FILE");
		fprintf(stderr, "  %-10s verbose output\n",             "-v");
		fprintf(stderr, "  %-10s output k-mers or SNP-mers in a text file\n",
				"-o FILE");

		free(bed_fns);
		return 1;
	}

	// warnings
	if (opt.snp) {
		if (opt.pre < F_VAL_INFO_BITS) {
			fprintf(stderr, "[E::%s] -p/--pre must not be smaller than %d\n", __func__, F_VAL_INFO_BITS);
			return 0;
		}
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
	if (opt.mko < 1) {
		fprintf(stderr, "[E::%s] mko is %d must be >= 1\n", __func__, opt.mko);
		return 1;
	}
	if (opt.filt_type == 2) {
		fprintf(stderr, "[W::%s] -f/--filt_type = 2: parameter --mko has no effect in this mode\n", __func__);
	}
	if (!opt.snp) {
		fprintf(stderr, "[W::%s] --snp not activated: parameter --maf has no effect in this mode\n", __func__);
	}

	// start
	h = pg_detect(argv + o.ind, n_bed ? bed_fns : 0, n_fa, &opt, fn_out);
	fprintf(stderr, "[M::%s] Detected required k-mers in %d files\n", __func__, n_fa);

	return 0;
}


int main_count(int argc, char *argv[])
{   
    pg_msht_t *h;
	char *fn_out = 0;
	char *bed_list_fn = 0;
	char **bed_fns = 0;
	char *kmer_file = 0;
	int n_bed = 0;
	int n_fa; 
	int c;
	pg_opt_t opt;
	ketopt_t o = KETOPT_INIT;
	pg_opt_init(&opt);

	static ko_longopt_t long_opts[] = {
        { "k_length", ko_required_argument, 301 },
        { "msf", ko_required_argument, 302 },
		{ "maf", ko_required_argument, 303 },
		{ "snp", ko_no_argument, 304 },
		{ "mko", ko_required_argument, 305 },
        { "pre", ko_required_argument, 306 },
        { "filt_type", ko_required_argument, 307 },
        { "chunk_size", ko_required_argument, 308 },
        { "threads", ko_required_argument, 309 },
		{ "write_info", ko_no_argument, 310 },
		{ "write_mko", ko_required_argument, 311 },
		{ "output", ko_required_argument, 312 },
        { "verbose", ko_no_argument, 313 },
		{ "bed", ko_required_argument, 314 },
		{ "kmers", ko_required_argument, 315 },
        { 0, 0, 0 }
    };

	while ((c = ketopt(&o, argc, argv, 1, "k:p:f:K:t:wo:vb:", long_opts)) >= 0) {
        if      (c == 'k' || c == 301) opt.k = atoi(o.arg);
        else if (c == 302) opt.msf = atof(o.arg);
		else if (c == 303) opt.maf = atof(o.arg);
		else if (c == 304) opt.snp = 1;
		else if (c == 305) opt.mko = atoi(o.arg);
        else if (c == 'p' || c == 306) opt.pre = atoi(o.arg);
        else if (c == 'f' || c == 307) opt.filt_type = atoi(o.arg);
        else if (c == 'K' || c == 308) opt.chunk_size = mm_parse_num(o.arg);
        else if (c == 't' || c == 309) opt.n_threads = atoi(o.arg);
		else if (c == 'w' || c == 310) opt.write_info = 1;
		else if (c == 311) opt.write_mko = atoi(o.arg);
		else if (c == 'o' || c == 312) fn_out = o.arg;
        else if (c == 'v' || c == 313) opt.verbose = 1;
		else if (c == 'b' || c == 314) bed_list_fn = o.arg;
		else if (c == 315) kmer_file = o.arg;
    }

	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: pgtools count [options] <in1.fa> [in2.fa [...]]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  %-10s k-mer size [%d]\n",            "-k INT",  opt.k);
		fprintf(stderr, "  %-10s minimum sample frequency of k-mers across inputs to be kept [%g]\n",
				"-t FLOAT", opt.msf);
		fprintf(stderr, "  %-10s minimum allelic frequency of SNP-mers across input samples [%g]\n",
				"-t FLOAT", opt.maf);
		fprintf(stderr, "  %-10s snp mode to count SNP-mers\n",         "--snp");
		fprintf(stderr, "  %-10s maximum k-mer occurrences [%d]\n",         "--mko INT",  opt.mko);
		fprintf(stderr, "  %-10s prefix length [%d]\n",         "-p INT",  opt.pre);
		fprintf(stderr, "  %-10s filter type [%d]\n",           "-f INT",  opt.filt_type);
		fprintf(stderr, "  %-10s number of worker threads [%d]\n",
				"-t INT",  opt.n_threads);
		fprintf(stderr, "  %-10s chunk size [%ld]\n",           "-K INT", (long)opt.chunk_size);
		fprintf(stderr, "  %-10s write pangenome hits in the info file up to --write_mko occurrences\n", "-w");
		fprintf(stderr, "  %-10s write up to INT pangenome hits in the info file [%d]\n", "--write_mko INT", opt.write_mko);
		fprintf(stderr, "  %-10s text file listing one BED path per line, one line per input FASTA, in the same order\n",
				"-b FILE");
		fprintf(stderr, "  %-10s verbose output\n",             "-v");
		fprintf(stderr, "  %-10s output genome-specific k-mers or SNP-mers in TSV format\n",
				"-o FILE");
		fprintf(stderr, "  %-10s list of k-mers or SNP-mers as a text file\n",
				"--kmers FILE");

		free(bed_fns);
		return 1;
	}

	// warnings
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
	if (opt.mko < 1) {
		fprintf(stderr, "[E::%s] mko is %d must be >= 1\n", __func__, opt.mko);
	}
	if (!opt.snp) {
		fprintf(stderr, "[W::%s] --snp not activated: parameter --maf has no effect in this mode\n", __func__);
	}

	// start
	// first step: detect specific k-mers
	if (kmer_file) {
		fprintf(stderr, "[W::%s] -s/--kmers given: k-mers are loaded directly from '%s'\n-f/--filt_type, --msf, --maf, and --mko have no effect in this mode (they only apply to pgtools detect)\n", __func__, kmer_file);
		// repopulate hash table from SNP-mers file
		fprintf(stderr, "[M::%s] repopulating hash table from file '%s'\n", __func__, kmer_file);
		h = pg_msht_repopulate(kmer_file, &opt);
	}
	else {
		h = pg_detect(argv + o.ind, n_bed ? bed_fns : 0, n_fa, &opt, 0);
		fprintf(stderr, "[M::%s] Detected required k-mers in %d files\n", __func__, n_fa);
	}

	// second step: count k-mers in each file
	// reduce chunk_size for second pass
	opt.chunk_size = opt.chunk_size / 10;
	if (opt.chunk_size < 1024*1024) opt.chunk_size = 1024*1024; // minimum 1MB

	if (fn_out == NULL) fn_out = "-"; // redirect output to stdout
	pg_count(argv + o.ind, n_bed ? bed_fns : 0, n_fa, h->n_ins_tot, &opt, h, fn_out);

	free(bed_fns);
	fprintf(stderr, "[M::%s] Counted required k-mers in %d files\n", __func__, n_fa);
	
    return 0;
}

int main(int argc, char *argv[])
{   
	int ret = 1;
    pg_reset_realtime();

	if (argc < 2) {
		fprintf(stderr, "Usage: pgtools <command> [options]\n");
		fprintf(stderr, "Commands:\n");
		fprintf(stderr, "  detect    detect specific k-mers in the passed fasta files\n");
		fprintf(stderr, "  count    count specific k-mers in the passed fasta files\n");
		return 1;
	}

    if (strcmp(argv[1], "detect") == 0) ret = main_detect(argc-1, argv+1);
	else if (strcmp(argv[1], "count") == 0) ret = main_count(argc-1, argv+1);
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