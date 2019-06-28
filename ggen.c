#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
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
	if (mg_dbg_flag & MG_DBG_QNAME)
		fprintf(stderr, "QR\t%s\t%d\t%d\n", s->seq[i].name, tid, s->seq[i].l_seq);
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
	if (mg_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] indexed the graph\n", __func__,
				realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));
	s->opt = opt;
	s->seq = mg_bseq_read(fp, 1ULL<<62, 0, 0, 0, &s->n_seq);
	if (mg_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] loaded file \"%s\"\n", __func__,
				realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), fn);
	for (i = 0; i < s->n_seq; ++i) s->seq[i].rid = i;
	KCALLOC(0, s->buf, n_threads);
	for (i = 0; i < n_threads; ++i) s->buf[i] = mg_tbuf_init();
	KCALLOC(0, s->gcs, s->n_seq);

	// mapping and graph generation
	kt_for(n_threads, worker_for, s, s->n_seq);
	if (mg_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequence(s) to the graph\n", __func__,
				realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), s->n_seq);
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
	free(s);
	mg_bseq_close(fp);
	return 0;
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
