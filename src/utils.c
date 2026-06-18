#include "utils.h"

unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

char nt4_seq_table[5] = {
    'A', // 0
    'C', // 1
    'G', // 2
    'T', // 3
    'N'  // 4 (invalid/ambiguous)
};


void pg_opt_init(pg_opt_t *o)
{
	memset(o, 0, sizeof(pg_opt_t));
	o->k = 31;
	o->min_freq = 0.95;
	o->pre = 10;
	o->filt_type = 2;
	o->n_threads = 4;
	o->chunk_size = mm_parse_num("1.9g");
	o->verbose = 0;
}

char *find_cli_tool(const char *tool) {
    char cmd[256];
    snprintf(cmd, sizeof(cmd), "bash -l -c 'which %s'", tool);

    FILE *fp = popen(cmd, "r");
    if (!fp) return NULL;

    char *path = malloc(1024);
    if (!fgets(path, 1024, fp)) {
        free(path);
        pclose(fp);
        return NULL;
    }
    pclose(fp);

    path[strcspn(path, "\n")] = '\0';
    return path;
}

int64_t mm_parse_num(const char *str)
{
    double x;
    char *p;
    x = strtod(str, &p);
    if      (*p == 'G' || *p == 'g') x *= 1e9;
    else if (*p == 'M' || *p == 'm') x *= 1e6;
    else if (*p == 'K' || *p == 'k') x *= 1e3;
    return (int64_t)(x + .499);
}