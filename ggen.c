#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "kthread.h"
#include "kalloc.h"
#include "mgpriv.h"
#include "bseq.h"

typedef struct {
    int n_seq;
	const mg_mapopt_t *opt;
	mg_idx_t *gi;
	mg_bseq1_t *seq;
	mg_gchains_t **gcs;
	mg_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *s = (step_t*)_data;
	s->gcs[i] = mg_map(s->gi, s->seq[i].l_seq, s->seq[i].seq, s->buf[tid], s->opt, s->seq[i].name);
}

int mg_ggen(gfa_t *g, const char *fn, const mg_idxopt_t *ipt, const mg_mapopt_t *opt, const mg_ggopt_t *go, int n_threads)
{
	int i;
	mg_bseq_file_t *fp;
	step_t *s;

	// read sequences
	fp = mg_bseq_open(fn);
	if (fp == 0) return -1;

	KCALLOC(0, s, 1);
	s->gi = mg_index_gfa(g, ipt->k, ipt->w, ipt->bucket_bits, n_threads);
	s->opt = opt;
	s->seq = mg_bseq_read(fp, 1ULL<<62, 0, 0, 0, &s->n_seq);
	for (i = 0; i < s->n_seq; ++i) s->seq[i].rid = i;
	KCALLOC(0, s->buf, n_threads);
	for (i = 0; i < n_threads; ++i) s->buf[i] = mg_tbuf_init();
	KCALLOC(0, s->gcs, s->n_seq);

	// mapping and graph generation
	kt_for(n_threads, worker_for, s, s->n_seq);
	mg_idx_destroy(s->gi);
	for (i = 0; i < n_threads; ++i) mg_tbuf_destroy(s->buf[i]);
	free(s->buf);

#if 0 // for debugging
	kstring_t str = {0,0,0};
	for (i = 0; i < s->n_seq; ++i) {
		mg_bseq1_t *t = &s->seq[i];
		mg_write_paf(&str, g, s->gcs[i], t->l_seq, t->name, opt->flag, 0);
		mg_err_fputs(str.s, stdout);
	}
	free(str.s);
#endif

	mg_ggsimple(0, go, g, s->n_seq, s->seq, s->gcs);

	// free the rest
	for (i = 0; i < s->n_seq; ++i) {
		mg_gchain_free(s->gcs[i]);
		free(s->seq[i].seq); free(s->seq[i].name);
	}
	free(s->gcs); free(s->seq);
	if (mg_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), s->n_seq);
	free(s);
	mg_bseq_close(fp);
	return 0;
}
