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
	int32_t l_seq, rid;
	char *name, *seq, *qual, *comment;
} mg_bseq1_t;

mg_bseq_file_t *mg_bseq_open(const char *fn);
void mg_bseq_close(mg_bseq_file_t *fp);
mg_bseq1_t *mg_bseq_read(mg_bseq_file_t *fp, int64_t chunk_size, int with_qual, int with_comment, int frag_mode, int *n_);
mg_bseq1_t *mg_bseq_read_frag(int n_fp, mg_bseq_file_t **fp, int64_t chunk_size, int with_qual, int with_comment, int *n_);
int mg_bseq_eof(mg_bseq_file_t *fp);

extern unsigned char seq_nt4_table[256];
extern unsigned char gfa_comp_table[256];

static inline int32_t mg_qname_len(const char *s)
{
	int32_t l;
	l = strlen(s);
	return l >= 3 && s[l-1] >= '0' && s[l-1] <= '9' && s[l-2] == '/'? l - 2 : l;
}

static inline int32_t mg_qname_same(const char *s1, const char *s2)
{
	int32_t l1, l2;
	l1 = mg_qname_len(s1);
	l2 = mg_qname_len(s2);
	return (l1 == l2 && strncmp(s1, s2, l1) == 0);
}

static inline void mg_toupper(int32_t len, char *seq)
{
	int32_t j;
	for (j = 0; j < len; ++j)
		seq[j] = seq[j] < 'a' || seq[j] > 'z'? seq[j] : seq[j] - 32;
}

static inline void mg_revcomp_seq(int32_t len, char *seq)
{
	int32_t i;
	for (i = 0; i < len>>1; ++i) {
		int32_t t = seq[len - i - 1];
		seq[len - i - 1] = gfa_comp_table[(uint8_t)seq[i]];
		seq[i] = gfa_comp_table[t];
	}
	if (len&1) seq[len>>1] = gfa_comp_table[(uint8_t)seq[len>>1]];
}

static inline void mg_revcomp_bseq(mg_bseq1_t *s)
{
	int32_t i, t, l = s->l_seq;
	mg_revcomp_seq(s->l_seq, s->seq);
	if (s->qual)
		for (i = 0; i < l>>1; ++i)
			t = s->qual[l - i - 1], s->qual[l - i - 1] = s->qual[i], s->qual[i] = t;
}

#ifdef __cplusplus
}
#endif

#endif
