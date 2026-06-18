#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parser.h"
#include "kseq.h"


KSEQ_INIT(gzFile, gzread)


// FASTA/FASTQ parser
struct bseq_file_s {
	gzFile fp;
	kseq_t *ks;
};

bseq_file_t *bseq_open(const char *fn)
{
	bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (f == 0) return 0;
	fp = calloc(1, sizeof(bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

void bseq_close(bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

bseq1_t *bseq_read(bseq_file_t *fp, int64_t chunk_size, int keep_comment, int *n_)
{
	int m, n;
	int64_t size = 0;
	bseq1_t *seqs;
	kseq_t *ks = fp->ks;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		bseq1_t *s;
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
		s->name = strdup(ks->name.s);
		s->comment = (ks->comment.s && keep_comment)? strdup(ks->comment.s) : 0;
		s->seq = strdup(ks->seq.s);
		s->qual = ks->qual.l? strdup(ks->qual.s) : 0;
		s->l_seq = ks->seq.l;
		s->aux = 0;
		size += seqs[n++].l_seq;
		if (size >= chunk_size) break;
	}
	*n_ = n;
	return seqs;
}


// PAF parser
paf_rec_t *parse_paf(const char *fn, int n_snps, int32_t k)
{
    paf_rec_t *recs = malloc(n_snps * sizeof(paf_rec_t));

    FILE *fp = fopen(fn, "r");

    char line[LINE_BUF_SIZE];
    int i = 0;
    for (int i = 0; i < n_snps && fgets(line, sizeof(line), fp); i++) {
		char chrom[MAX_NAME_LEN];
		long pos;
		if (sscanf(line, "%*s\t%*s\t%*s\t%*s\t%*s\t%s\t%*s\t%ld", chrom, &pos) == 2) {
			strncpy(recs[i].chrom_name, chrom, MAX_NAME_LEN - 1);
			recs[i].chrom_name[MAX_NAME_LEN - 1] = '\0';
			recs[i].target_pos = pos + ((k >> 1) + 1);
		}
	}

    fclose(fp);
    return recs;
}