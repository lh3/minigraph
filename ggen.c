#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "kthread.h"
#include "kalloc.h"
#include "mgpriv.h"
#include "bseq.h"

typedef struct {
	int n_threads;
	const mg_mapopt_t *opt;
	const mg_ggopt_t *go;
	mg_bseq_file_t *fp;
	const mg_idx_t *gi;
	kstring_t str;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int n_seq;
	mg_bseq1_t *seq;
	mg_gchains_t **gcs;
	mg_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *s = (step_t*)_data;
	s->gcs[i] = mg_map(s->p->gi, s->seq[i].l_seq, s->seq[i].seq, s->buf[tid], s->p->opt, s->seq[i].name);
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = mg_bseq_read(p->fp, 1ULL<<62, 0, 0, 0, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i) s->seq[i].rid = i;
			s->buf = (mg_tbuf_t**)calloc(p->n_threads, sizeof(mg_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i) s->buf[i] = mg_tbuf_init();
			s->gcs = KCALLOC(0, mg_gchains_t*, s->n_seq);
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_seq);
		return in;
    } else if (step == 2) { // step 2: output
		void *km = 0;
        step_t *s = (step_t*)in;
		for (i = 0; i < p->n_threads; ++i) mg_tbuf_destroy(s->buf[i]);
		free(s->buf);
		if (!(mg_dbg_flag & MG_DBG_NO_KALLOC)) km = km_init();
		//mg_ggsimple(km, p->go, p->gi->g, s->n_seq, s->seq, s->gcs);
		for (i = 0; i < s->n_seq; ++i) {
			mg_bseq1_t *t = &s->seq[i];
			mg_write_paf(&p->str, p->gi->g, s->gcs[i], t->l_seq, t->name, p->opt->flag, km);
			mg_err_fputs(p->str.s, stdout);
			mg_gchain_free(s->gcs[i]);
			free(s->seq[i].seq); free(s->seq[i].name);
		}
		free(s->gcs); free(s->seq);
		if (km) km_destroy(km);
		if (mg_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), s->n_seq);
		free(s);
	}
    return 0;
}

int mg_ggen(const mg_idx_t *idx, const char *fn, const mg_mapopt_t *opt, const mg_ggopt_t *go, int n_threads)
{
	pipeline_t pl;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.fp = mg_bseq_open(fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.go = go, pl.gi = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	kt_pipeline(3, worker_pipeline, &pl, 3);
	free(pl.str.s);
	mg_bseq_close(pl.fp);
	return 0;
}
