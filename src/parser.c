#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parser.h"
#include "kseq.h"
#include "utils.h"

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



// paf_idx_t *idx_paf(const char *fn, int32_t k)
// {
//     paf_idx_t *m = paf_idx_init();
//     FILE *fp = fopen(fn, "r");

//     char line[LINE_BUF_SIZE];
//     while (1) {
//         off_t off = ftello(fp);
//         if (!fgets(line, sizeof(line), fp)) break;

//         char qname[MAX_NAME_LEN], chrom[MAX_NAME_LEN];
//         if (sscanf(line, "%s\t%*s\t%*s\t%*s\t%*s\t%s", qname, chrom) != 2)
//             continue;

//         uint64_t x, mask = (1ULL << k*2) - 1;
//         int i, valid = 1;
//         for (i = 0, x = 0; i < k; ++i) {
//             int c = seq_nt4_table[(uint8_t)qname[i]];
//             if (c >= 4) { valid = 0; break; }
//             x = (x << 2 | c) & mask;
//         }
//         if (!valid) continue;

//         uint64_t y = pg_hash64(x, mask);

//         int absent;
//         khint_t b = paf_idx_put(m, y, &absent);
//         if (absent) {
//             kh_val(m, b).offset = off;
//             strncpy(kh_val(m, b).chrom_name, chrom, MAX_NAME_LEN - 1);
//             kh_val(m, b).chrom_name[MAX_NAME_LEN - 1] = '\0';
//         }
//     }
//     fclose(fp);
//     return m;
// }

// paf_rec_t fetch_paf_rec(paf_idx_t *m, FILE *fp, uint64_t hflanks, int32_t k)
// {
//     paf_rec_t rec = {0};
//     khint_t b = paf_idx_get(m, hflanks);
//     if (b == kh_end(m)) return rec;

//     fseeko(fp, kh_val(m, b).offset, SEEK_SET);
//     char line[LINE_BUF_SIZE];
//     fgets(line, sizeof(line), fp);

//     char chrom[MAX_NAME_LEN];
//     long pos;
//     if (sscanf(line, "%*s\t%*s\t%*s\t%*s\t%*s\t%s\t%*s\t%ld", chrom, &pos) == 2) {
//         strncpy(rec.chrom_name, chrom, MAX_NAME_LEN - 1);
//         rec.target_pos = pos + ((k >> 1) + 1);
//     }

//     return rec;
// }