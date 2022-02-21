#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define __STDC_LIMIT_MACROS
#include "bseq.h"
#include "kvec-km.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define CHECK_PAIR_THRES 1000000

struct mg_bseq_file_s {
	gzFile fp;
	kseq_t *ks;
	mg_bseq1_t s;
};

mg_bseq_file_t *mg_bseq_open(const char *fn)
{
	mg_bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (f == 0) return 0;
	fp = (mg_bseq_file_t*)calloc(1, sizeof(mg_bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

void mg_bseq_close(mg_bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

static inline char *kstrdup(const kstring_t *s)
{
	char *t;
	t = (char*)malloc(s->l + 1);
	memcpy(t, s->s, s->l + 1);
	return t;
}

static inline void kseq2bseq(kseq_t *ks, mg_bseq1_t *s, int with_qual, int with_comment)
{
	int i;
	if (ks->name.l == 0)
		fprintf(stderr, "[WARNING]\033[1;31m empty sequence name in the input.\033[0m\n");
	s->name = kstrdup(&ks->name);
	s->seq = kstrdup(&ks->seq);
	for (i = 0; i < (int)ks->seq.l; ++i) // convert U to T
		if (s->seq[i] == 'u' || s->seq[i] == 'U')
			--s->seq[i];
	s->qual = with_qual && ks->qual.l? kstrdup(&ks->qual) : 0;
	s->comment = with_comment && ks->comment.l? kstrdup(&ks->comment) : 0;
	s->l_seq = ks->seq.l;
}

mg_bseq1_t *mg_bseq_read(mg_bseq_file_t *fp, int64_t chunk_size, int with_qual, int with_comment, int frag_mode, int *n_)
{
	int64_t size = 0;
	int ret;
	kvec_t(mg_bseq1_t) a = {0,0,0};
	kseq_t *ks = fp->ks;
	*n_ = 0;
	if (fp->s.seq) {
		kv_resize(mg_bseq1_t, 0, a, 256);
		kv_push(mg_bseq1_t, 0, a, fp->s);
		size = fp->s.l_seq;
		memset(&fp->s, 0, sizeof(mg_bseq1_t));
	}
	while ((ret = kseq_read(ks)) >= 0) {
		mg_bseq1_t *s;
		assert(ks->seq.l <= INT32_MAX);
		if (a.m == 0) kv_resize(mg_bseq1_t, 0, a, 256);
		kv_pushp(mg_bseq1_t, 0, a, &s);
		kseq2bseq(ks, s, with_qual, with_comment);
		size += s->l_seq;
		if (size >= chunk_size) {
			if (frag_mode && a.a[a.n-1].l_seq < CHECK_PAIR_THRES) {
				while (kseq_read(ks) >= 0) {
					kseq2bseq(ks, &fp->s, with_qual, with_comment);
					if (mg_qname_same(fp->s.name, a.a[a.n-1].name)) {
						kv_push(mg_bseq1_t, 0, a, fp->s);
						memset(&fp->s, 0, sizeof(mg_bseq1_t));
					} else break;
				}
			}
			break;
		}
	}
	if (ret < -1)
		fprintf(stderr, "[WARNING]\033[1;31m wrong FASTA/FASTQ record. Continue anyway.\033[0m\n");
	*n_ = a.n;
	return a.a;
}

mg_bseq1_t *mg_bseq_read_frag(int n_fp, mg_bseq_file_t **fp, int64_t chunk_size, int with_qual, int with_comment, int *n_)
{
	int i;
	int64_t size = 0;
	kvec_t(mg_bseq1_t) a = {0,0,0};
	*n_ = 0;
	if (n_fp < 1) return 0;
	while (1) {
		int n_read = 0;
		for (i = 0; i < n_fp; ++i)
			if (kseq_read(fp[i]->ks) >= 0)
				++n_read;
		if (n_read < n_fp) {
			if (n_read > 0)
				fprintf(stderr, "[W::%s]\033[1;31m query files have different number of records; extra records skipped.\033[0m\n", __func__);
			break; // some file reaches the end
		}
		if (a.m == 0) kv_resize(mg_bseq1_t, 0, a, 256);
		for (i = 0; i < n_fp; ++i) {
			mg_bseq1_t *s;
			kv_pushp(mg_bseq1_t, 0, a, &s);
			kseq2bseq(fp[i]->ks, s, with_qual, with_comment);
			size += s->l_seq;
		}
		if (size >= chunk_size) break;
	}
	*n_ = a.n;
	return a.a;
}

int mg_bseq_eof(mg_bseq_file_t *fp)
{
	return (ks_eof(fp->ks->f) && fp->s.seq == 0);
}
