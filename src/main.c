#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* for fork */
#include <sys/types.h> /* for pid_t */
#include <sys/wait.h> /* for wait */

#include "htab.h"
#include "ketopt.h"
#include "sys.h"
#include "utils.h"

double *main_count(int argc, char *argv[])
{   
    pg_mht_t *h;
	// pg_csr_t *csr;
	char *fn_snps = 0;
	char *fn_ropebwt3 = 0;
	char *fn_out = 0;
	char *index_fn = 0;
	int c, i;
	pg_opt_t opt;
	ketopt_t o = KETOPT_INIT;
	pg_opt_init(&opt);

	static ko_longopt_t long_opts[] = {
        { "kmer", ko_required_argument, 301 },
        { "min-freq", ko_required_argument, 302 },
        { "pre", ko_required_argument, 303 },
        { "filt-type", ko_required_argument, 304 },
        { "chunk-size", ko_required_argument, 305 },
        { "threads", ko_required_argument, 306 },
        { "index", ko_required_argument, 307 },
        { "output1", ko_required_argument, 308 },
		{ "output2", ko_required_argument, 309 },
		{ "output", ko_required_argument, 310 },
        { "verbose", ko_no_argument, 311 },
        { 0, 0, 0 }
    };

	while ((c = ketopt(&o, argc, argv, 1, "k:m:p:f:K:t:i:a:b:o:v", long_opts)) >= 0) {
        if      (c == 'k' || c == 301) opt.k = atoi(o.arg);
        else if (c == 'm' || c == 302) opt.min_freq = atof(o.arg);
        else if (c == 'p' || c == 303) opt.pre = atoi(o.arg);
        else if (c == 'f' || c == 304) opt.filt_type = atoi(o.arg);
        else if (c == 'K' || c == 305) opt.chunk_size = mm_parse_num(o.arg);
        else if (c == 't' || c == 306) opt.n_threads = atoi(o.arg);
        else if (c == 'i' || c == 307) index_fn = o.arg;
        else if (c == 'a' || c == 308) fn_snps = o.arg;
		else if (c == 'b' || c == 309) fn_ropebwt3 = o.arg;
		else if (c == 'o' || c == 310) fn_out = o.arg;
        else if (c == 'v' || c == 311) opt.verbose = 1;
    }

	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: pgtools count [options] <in1.fa> [in2.fa [...]]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -m FLOAT   minimum frequency of k-mers across inputs to be kept [%g]\n", opt.min_freq);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", opt.pre);
		fprintf(stderr, "  -f INT	  filter type [%d]\n", opt.filt_type);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_threads);
		fprintf(stderr, "  -K INT     chunk size [1.9g]\n");
		fprintf(stderr, "  -i FILE    index file to run ropebwt3\n");
		fprintf(stderr, "  -v         verbose output\n");
		fprintf(stderr, "  -a FILE   output pangenome SNP-mers\n");
		fprintf(stderr, "  -b FILE   output ropebwt3 paf file for SNP-mers\n");
		fprintf(stderr, "  -o FILE   output genome-specific SNPs in VCF format\n");
		return 1;
	}
	if (opt.k >= 32) {
		fprintf(stderr, "ERROR: -k must be <=31\n");
		return 1;
	}

	// first step: count k-mers in the input files argv and filter for SNP-mers
	h = pg_count(argv + o.ind, argc - o.ind, &opt, fn_snps);

    // second step: find SNPs positions with ropebwt3
	double *rb3_stats = calloc(3, sizeof(double));
	if (fn_ropebwt3) {
		if (!index_fn) {
			fprintf(stderr, "ERROR: -i index file is needed to run ropebwt3\n");
			return 1;
		}

		char *ropebwt3_path = find_cli_tool("ropebwt3");
		if (!ropebwt3_path) {
			fprintf(stderr, "ERROR: ropebwt3 not found in PATH (including ~/.bashrc)\n");
			return 1;
		}
		fprintf(stderr, "[M::main_count] Found ropebwt3 at: %s\n", ropebwt3_path);

		fprintf(stderr, "[M::main_count] Searching the SNP-mers against the FM-index using ropebwt3\n");

		char arg_t[32], arg_m[32], arg_p[32];
		snprintf(arg_t, sizeof(arg_t), "-t%d", opt.n_threads);
		snprintf(arg_m, sizeof(arg_m), "-m%d", opt.k);
		snprintf(arg_p, sizeof(arg_p), "-p%d", argc - o.ind);
		
		int pipe_fd[2];
		pipe(pipe_fd);

		pid_t pid = fork();
		if (pid == 0) { // child process
			FILE *out = fopen(fn_ropebwt3, "w");
			if (!out) { perror("fopen failed"); exit(127); }
			dup2(fileno(out), STDOUT_FILENO);
			fclose(out);
			dup2(pipe_fd[1], STDERR_FILENO);
			close(pipe_fd[0]);
			close(pipe_fd[1]);

			char *child_argv[] = {
				ropebwt3_path,
				"sw",
				"-N100",
				"-L",
				"-b",
				arg_t,
				arg_m,
				arg_p,
				index_fn,
				fn_snps,
				NULL
			};
			execv(ropebwt3_path, child_argv); // use exact path found by which
			perror("execv failed");
			exit(127);
		} else {
			close(pipe_fd[1]);

			FILE *pipe_read = fdopen(pipe_fd[0], "r");
			char line[1024];
			while (fgets(line, sizeof(line), pipe_read)) {
				double rt, cpu, rss;
				if (sscanf(line, "[M::main] Real time: %lf sec; CPU: %lf sec; Peak RSS: %lf GB",
						&rt, &cpu, &rss) == 3) {
					rb3_stats[0] = rt;
					rb3_stats[1] = cpu;
					rb3_stats[2] = rss;
				} 
				// else {
					// fputs(line, stderr);
				// }
			}
			fclose(pipe_read);
			waitpid(pid, 0, 0);
		}
		free(ropebwt3_path);
	}

	// third step: count SNPmers in each file
	pg_findsnp(argv + o.ind, argc - o.ind, h->n_ins_tot, &opt, h, fn_ropebwt3, fn_out);

	fprintf(stderr, "[M::%s] Analyzed %d files\n", __func__, argc - o.ind);
	
    return rb3_stats;
}


int main(int argc, char *argv[])
{   
	double *rb3_stats;
    pg_reset_realtime();
    if (strcmp(argv[1], "count") == 0) rb3_stats = main_count(argc-1, argv+1);

    if (rb3_stats != NULL) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, PG_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (int i = 0; i < argc; ++i) {
			fprintf(stderr, " %s", argv[i]);
		}
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, pg_realtime() + rb3_stats[0], pg_cputime() + rb3_stats[1], (pg_peakrss() / 1024.0 / 1024.0 / 1024.0 > rb3_stats[2]) ? (pg_peakrss() / 1024.0 / 1024.0 / 1024.0) : (rb3_stats[2]));
	}

    return 0;
}