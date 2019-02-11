#ifndef MM_BSEQ_H
#define MM_BSEQ_H

#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

struct mg_bseq_file_s;
typedef struct mg_bseq_file_s mg_bseq_file_t;

typedef struct {
	int l_seq, rid;
	char *name, *seq, *qual, *comment;
} mg_bseq1_t;

mg_bseq_file_t *mg_bseq_open(const char *fn);
void mg_bseq_close(mg_bseq_file_t *fp);
mg_bseq1_t *mg_bseq_read(mg_bseq_file_t *fp, int chunk_size, int with_qual, int with_comment, int frag_mode, int *n_);
mg_bseq1_t *mg_bseq_read_frag(int n_fp, mg_bseq_file_t **fp, int chunk_size, int with_qual, int with_comment, int *n_);
int mg_bseq_eof(mg_bseq_file_t *fp);

extern unsigned char seq_nt4_table[256];
extern unsigned char seq_comp_table[256];

static inline int mg_qname_len(const char *s)
{
	int l;
	l = strlen(s);
	return l >= 3 && s[l-1] >= '0' && s[l-1] <= '9' && s[l-2] == '/'? l - 2 : l;
}

static inline int mg_qname_same(const char *s1, const char *s2)
{
	int l1, l2;
	l1 = mg_qname_len(s1);
	l2 = mg_qname_len(s2);
	return (l1 == l2 && strncmp(s1, s2, l1) == 0);
}

static inline void mg_revcomp_bseq(mg_bseq1_t *s)
{
	int i, t, l = s->l_seq;
	for (i = 0; i < l>>1; ++i) {
		t = s->seq[l - i - 1];
		s->seq[l - i - 1] = seq_comp_table[(uint8_t)s->seq[i]];
		s->seq[i] = seq_comp_table[t];
	}
	if (l&1) s->seq[l>>1] = seq_comp_table[(uint8_t)s->seq[l>>1]];
	if (s->qual)
		for (i = 0; i < l>>1; ++i)
			t = s->qual[l - i - 1], s->qual[l - i - 1] = s->qual[i], s->qual[i] = t;
}

#ifdef __cplusplus
}
#endif

#endif
