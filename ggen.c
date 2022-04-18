#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "kthread.h"
#include "kalloc.h"
#include "sys.h"
#include "bseq.h"
#include "ggen.h"
#include "mgpriv.h"
#include "gfa-priv.h"

typedef struct {
	int n_seq;
	mg_bseq1_t *seq;
	mg_gchains_t **gcs;
} maprst_t;

typedef struct {
	const mg_mapopt_t *opt;
	const mg_idx_t *gi;
	mg_tbuf_t **buf;
	maprst_t *r;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *s = (step_t*)_data;
	if (mg_dbg_flag & MG_DBG_QNAME)
		fprintf(stderr, "QR\t%s\t%d\t%d\n", s->r->seq[i].name, tid, s->r->seq[i].l_seq);
	if ((s->opt->flag & MG_M_SKIP_GCHECK) == 0 && mg_verbose >= 2) {
		if (gfa_sseq_get(s->gi->g, s->r->seq[i].name) >= 0)
			fprintf(stderr, "[W::%s] stable sequence \"%s\" already present in the graph. This will lead to inconsistent rGFA.\n",
					__func__, s->r->seq[i].name);
	}
	s->r->gcs[i] = mg_map(s->gi, s->r->seq[i].l_seq, s->r->seq[i].seq, s->buf[tid], s->opt, s->r->seq[i].name);
}

static maprst_t *ggen_map(const mg_idx_t *gi, const mg_mapopt_t *opt, const char *fn, int n_threads)
{
	mg_bseq_file_t *fp;
	maprst_t *r;
	step_t s;
	int i;

	fp = mg_bseq_open(fn);
	if (fp == 0) return 0;

	KCALLOC(0, r, 1);
	r->seq = mg_bseq_read(fp, 1ULL<<62, 0, 0, 0, &r->n_seq);
	mg_bseq_close(fp);
	if (mg_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] loaded file \"%s\"\n", __func__,
				realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), fn);
	for (i = 0; i < r->n_seq; ++i) {
		r->seq[i].rid = i;
		mg_toupper(r->seq[i].l_seq, r->seq[i].seq);
	}
	KCALLOC(0, r->gcs, r->n_seq);

	s.gi = gi, s.opt = opt, s.r = r;
	KCALLOC(0, s.buf, n_threads);
	for (i = 0; i < n_threads; ++i) s.buf[i] = mg_tbuf_init();
	kt_for(n_threads, worker_for, &s, r->n_seq);
	if (mg_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequence(s) to the graph\n", __func__,
				realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), r->n_seq);
	for (i = 0; i < n_threads; ++i) mg_tbuf_destroy(s.buf[i]);
	free(s.buf);
	return r;
}

static void mg_free_maprst(maprst_t *r)
{
	int i;
	for (i = 0; i < r->n_seq; ++i) {
		mg_gchain_free(r->gcs[i]);
		free(r->seq[i].seq); free(r->seq[i].name);
	}
	free(r->gcs); free(r->seq);
	free(r);
}

int mg_ggen_aug(gfa_t *g, int32_t n_fn, const char **fn, const mg_idxopt_t *ipt, const mg_mapopt_t *opt0, const mg_ggopt_t *go, int n_threads)
{
	int i;
	mg_mapopt_t opt = *opt0;
	if (g == 0) return -1;
	for (i = 0; i < n_fn; ++i) {
		mg_idx_t *gi;
		maprst_t *r;
		if ((gi = mg_index(g, ipt, n_threads, &opt)) == 0) return -1;
		r = ggen_map(gi, &opt, fn[i], n_threads);
		if (opt0->flag & MG_M_CIGAR)
			mg_ggsimple_cigar(0, go, g, r->n_seq, r->seq, r->gcs);
		else
			mg_ggsimple(0, go, g, r->n_seq, r->seq, r->gcs);
		mg_free_maprst(r);
		mg_idx_destroy(gi);
	}
	return 0;
}

int mg_ggen_cov(gfa_t *g, int32_t n_fn, const char **fn, const mg_idxopt_t *ipt, const mg_mapopt_t *opt0, const mg_ggopt_t *go, int n_threads)
{
	int32_t i;
	mg_mapopt_t opt = *opt0;
	mg_idx_t *gi;
	double *cov_seg, *cov_link;
	int64_t j;
	if ((gi = mg_index(g, ipt, n_threads, &opt)) == 0) return -1;
	KCALLOC(0, cov_seg,  g->n_seg);
	KCALLOC(0, cov_link, g->n_arc);
	for (i = 0; i < n_fn; ++i) {
		maprst_t *r;
		r = ggen_map(gi, &opt, fn[i], n_threads);
		mg_cov_asm(g, r->n_seq, r->gcs, go->min_mapq, go->min_map_len, cov_seg, cov_link);
		mg_free_maprst(r);
	}
	mg_idx_destroy(gi);
	for (j = 0; j < g->n_seg; ++j) cov_seg[j] /= n_fn;
	for (j = 0; j < g->n_arc; ++j) cov_link[j] /= n_fn;
	gfa_aux_update_cv(g, "cf", cov_seg, cov_link);
	free(cov_seg); free(cov_link);
	return 0;
}

int mg_ggen_call(gfa_t *g, const char *fn, const mg_idxopt_t *ipt, const mg_mapopt_t *opt0, const mg_ggopt_t *go, int n_threads)
{
	mg_mapopt_t opt = *opt0;
	mg_idx_t *gi;
	maprst_t *r;
	if ((gi = mg_index(g, ipt, n_threads, &opt)) == 0) return -1;
	r = ggen_map(gi, &opt, fn, n_threads);
	mg_call_asm(g, r->n_seq, r->seq, r->gcs, go->min_mapq, go->min_map_len);
	mg_free_maprst(r);
	mg_idx_destroy(gi);
	return 0;
}

int mg_ggen(gfa_t *g, int32_t n_fn, const char **fn, const mg_idxopt_t *ipt, const mg_mapopt_t *opt, const mg_ggopt_t *go, int n_threads)
{
	if (go->flag & MG_G_CALL) return mg_ggen_call(g, fn[0], ipt, opt, go, n_threads);
	else if (go->flag & MG_G_CAL_COV) return mg_ggen_cov(g, n_fn, fn, ipt, opt, go, n_threads);
	else return mg_ggen_aug(g, n_fn, fn, ipt, opt, go, n_threads);
}

int32_t mg_path2seq(void *km, const gfa_t *g, const mg_gchains_t *gcs, int32_t ls, int32_t le, int32_t voff[2], char **seq_, int32_t *cap_) // NB: [ls,le] is a CLOSED interval
{
	extern unsigned char gfa_comp_table[256];
	int32_t i, k, l = 0, cap = *cap_;
	char *seq = *seq_;
	assert(0 <= ls && ls <= le && le < gcs->n_lc);
	for (k = ls; k <= le; ++k) {
		uint32_t v = gcs->lc[k].v, len = g->seg[v>>1].len;
		int32_t st = 0, en = len, tmp;
		if (k == ls) st = voff[0];
		if (k == le) en = voff[1];
		assert(0 <= st && st <= en && en <= len);
		if (en - st + l + 1 > cap) {
			cap = en - st + l + 1;
			kroundup32(cap);
			KREALLOC(km, seq, cap);
		}
		if (v&1) {
			uint8_t *ss = (uint8_t*)g->seg[v>>1].seq;
			tmp = st, st = len - en, en = len - tmp;
			for (i = en - 1; i >= st; --i)
				seq[l++] = gfa_comp_table[ss[i]];
		} else {
			memcpy(&seq[l], &g->seg[v>>1].seq[st], en - st);
			l += en - st;
		}
	}
	if (l == 0 && cap == 0) {
		cap = 8;
		KREALLOC(km, seq, cap);
	}
	seq[l] = 0;
	*seq_ = seq, *cap_ = cap;
	return l;
}
