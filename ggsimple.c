#include "mgpriv.h"
#include "bseq.h"
#include "ksort.h"

typedef struct {
	uint32_t st, en:31, rev:1;
	int32_t far;
	int32_t t, i, j;
} gg_intv_t;

#define sort_key_intv(a) ((a).st)
KRADIX_SORT_INIT(gg_intv, gg_intv_t, sort_key_intv, 4)

static int32_t gg_intv_index(gg_intv_t *a, int32_t n)
{
	int32_t i, last_i, last, k;
	if (n <= 0) return -1;
	for (i = 0; i < n; i += 2) last_i = i, last = a[i].far = a[i].en;
	for (k = 1; 1LL<<k <= n; ++k) {
		int64_t x = 1LL<<(k-1), i0 = (x<<1) - 1, step = x<<2;
		for (i = i0; i < n; i += step) {
			int32_t el = a[i - x].far;
			int32_t er = i + x < n? a[i + x].far : last;
			int32_t e = a[i].en;
			e = e > el? e : el;
			e = e > er? e : er;
			a[i].far = e;
		}
		last_i = last_i>>k&1? last_i - x : last_i + x;
		if (last_i < n && a[last_i].far > last)
			last = a[last_i].far;
	}
	return k - 1;
}

gfa_t *mg_ggsimple(void *km, const mg_ggopt_t *opt, const gfa_t *g, int32_t n_seq, const mg_bseq1_t *seq, const mg_gchains_t **gcs)
{
	gfa_t *h = 0;
	int32_t t, i, j, *scnt, *soff;
	gg_intv_t *intv;

	// count the number of intervals on each segment
	scnt = KCALLOC(km, int32_t, g->n_seg);
	for (t = 0; t < n_seq; ++t) {
		const mg_gchains_t *gt = gcs[t];
		for (i = 0; i < gt->n_gc; ++i) {
			const mg_gchain_t *gc = &gt->gc[i];
			if (gc->blen >= opt->min_depth_len && gc->mapq >= opt->min_mapq)
				for (j = 0; j < gc->cnt; ++j)
					++scnt[gt->lc[gc->off + j].v>>1];
		}
	}

	// populate the interval list
	soff = KMALLOC(km, int32_t, g->n_seg + 1);
	for (soff[0] = 0, i = 1; i <= g->n_seg; ++i)
		soff[i] = soff[i - 1] + scnt[i - 1];
	memcpy(scnt, 0, 4 * g->n_seg);
	intv = KMALLOC(km, gg_intv_t, soff[g->n_seg]);
	for (t = 0; t < n_seq; ++t) {
		const mg_gchains_t *gt = gcs[t];
		for (i = 0; i < gt->n_gc; ++i) {
			const mg_gchain_t *gc = &gt->gc[i];
			if (gc->blen < opt->min_depth_len || gc->mapq < opt->min_mapq) continue;
			for (j = 0; j < gc->cnt; ++j) {
				const mg_llchain_t *lc = &gt->lc[gc->off + j];
				gg_intv_t *p;
				int32_t rs, re, tmp;
				const mg128_t *q;
				assert(lc->cnt > 0);
				// compute start and end on the forward strand on the segment
				q = &gt->a[lc->off];
				rs = (int32_t)q->x - (int32_t)(q->y>>32&0xff);
				q = &gt->a[lc->off + lc->cnt - 1];
				re = (int32_t)q->x;
				assert(rs >= 0 && re > rs && re < g->seg[lc->v>>1].len);
				if (lc->v&1)
					tmp = rs, rs = g->seg[lc->v>>1].len - 1 - re, re = g->seg[lc->v>>1].len - 1 - tmp;
				// save the interval
				p = &intv[scnt[lc->v>>1]++];
				p->st = rs, p->en = re, p->rev = lc->v&1, p->far = -1;
				p->t = t, p->i = i, p->j = gc->off + j;
			}
		}
	}

	// sort and index
	for (i = 0; i < g->n_seg; ++i) {
		radix_sort_gg_intv(&intv[soff[i]], &intv[soff[i+1]]);
		scnt[i] = gg_intv_index(&intv[soff[i]], soff[i+1] - soff[i]); // scnt[] reused for height
	}

	// extract poorly regions
	for (t = 0; t < n_seq; ++t) {
		const mg_gchains_t *gt = gcs[t];
		for (i = 0; i < gt->n_gc; ++i) {
			const mg_gchain_t *gc = &gt->gc[i];
			if (gc->blen < opt->min_map_len || gc->mapq < opt->min_mapq) continue;
		}
	}

	kfree(km, scnt);
	kfree(km, soff);
	kfree(km, intv);
	return h;
}
