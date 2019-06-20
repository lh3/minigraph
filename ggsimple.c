#include "mgpriv.h"
#include "bseq.h"
#include "ksort.h"
#include "mss.h"

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

gfa_t *mg_ggsimple(void *km, const mg_ggopt_t *opt, const gfa_t *g, int32_t n_seq, const mg_bseq1_t *seq, mg_gchains_t *const* gcs)
{
	gfa_t *h = 0;
	int32_t t, i, j, *scnt, *soff, max_acnt;
	int64_t sum_acnt, sum_alen, *sc;
	gg_intv_t *intv;
	double a_dens;

	// count the number of intervals on each segment
	KCALLOC(km, scnt, g->n_seg);
	for (t = 0, max_acnt = 0; t < n_seq; ++t) {
		const mg_gchains_t *gt = gcs[t];
		for (i = 0; i < gt->n_gc; ++i) {
			const mg_gchain_t *gc = &gt->gc[i];
			if (gc->blen < opt->min_depth_len || gc->mapq < opt->min_mapq) continue;
			if (gc->n_anchor > max_acnt) max_acnt = gc->n_anchor;
			for (j = 0; j < gc->cnt; ++j)
				++scnt[gt->lc[gc->off + j].v>>1];
		}
	}
	if (max_acnt == 0) { // no gchain
		kfree(km, scnt);
		return 0;
	}

	// populate the interval list
	KMALLOC(km, soff, g->n_seg + 1);
	for (soff[0] = 0, i = 1; i <= g->n_seg; ++i)
		soff[i] = soff[i - 1] + scnt[i - 1];
	memset(scnt, 0, 4 * g->n_seg);
	KMALLOC(km, intv, soff[g->n_seg]);
	sum_acnt = sum_alen = 0;
	for (t = 0; t < n_seq; ++t) {
		const mg_gchains_t *gt = gcs[t];
		for (i = 0; i < gt->n_gc; ++i) {
			const mg_gchain_t *gc = &gt->gc[i];
			if (gc->blen < opt->min_depth_len || gc->mapq < opt->min_mapq) continue;
			for (j = 0; j < gc->cnt; ++j) {
				const mg_llchain_t *lc = &gt->lc[gc->off + j];
				gg_intv_t *p;
				int32_t rs, re, tmp;
				if (lc->cnt > 0) { // compute start and end on the forward strand on the segment
					const mg128_t *q = &gt->a[lc->off];
					rs = (int32_t)q->x - (int32_t)(q->y>>32&0xff);
					q = &gt->a[lc->off + lc->cnt - 1];
					re = (int32_t)q->x;
					assert(rs >= 0 && re > rs && re < g->seg[lc->v>>1].len);
					sum_alen += re - rs, sum_acnt += (q->x>>32) - (gt->a[lc->off].x>>32) + 1;
					if (lc->v&1)
						tmp = rs, rs = g->seg[lc->v>>1].len - 1 - re, re = g->seg[lc->v>>1].len - 1 - tmp;
				} else rs = 0, re = g->seg[lc->v>>1].len;
				// save the interval
				p = &intv[scnt[lc->v>>1]++];
				p->st = rs, p->en = re, p->rev = lc->v&1, p->far = -1;
				p->t = t, p->i = i, p->j = gc->off + j;
			}
		}
	}
	a_dens = (double)sum_acnt / sum_alen;

	// sort and index intervals
	for (i = 0; i < g->n_seg; ++i) {
		radix_sort_gg_intv(&intv[soff[i]], &intv[soff[i+1]]);
		scnt[i] = gg_intv_index(&intv[soff[i]], soff[i+1] - soff[i]); // scnt[] reused for height
	}

	// extract poorly regions
	KMALLOC(km, sc, max_acnt);
	for (t = 0; t < n_seq; ++t) {
		const mg_gchains_t *gt = gcs[t];
		for (i = 0; i < gt->n_gc; ++i) {
			const mg_gchain_t *gc = &gt->gc[i];
			int32_t off_a, off_l, n_ss;
			msseg_t *ss;
			if (gc->blen < opt->min_map_len || gc->mapq < opt->min_mapq) continue;
			assert(gc->cnt > 0);

			// fill sc[]. This part achieves a similar goal to the one in mg_gchain_extra(). It makes more assumptions, but is logically simpler.
			off_l = gc->off;
			off_a = gt->lc[off_l].off + 1;
			for (j = 1; j < gc->n_anchor; ++j, ++off_a) {
				const mg128_t *q = &gt->a[off_a - 1], *p = &gt->a[off_a];
				const mg_llchain_t *lc = &gt->lc[off_l];
				int32_t pd, qd = (int32_t)p->y - (int32_t)q->y, c = (int32_t)(p->x>>32) - (int32_t)(q->x>>32) - 1;
				if (off_a == lc->off + lc->cnt) { // we are at the end of the current lchain
					pd = g->seg[lc->v>>1].len - (int32_t)q->x - 1;
					for (++off_l; off_l < gc->off + gc->cnt && gt->lc[off_l].cnt == 0; ++off_l)
						pd += g->seg[gt->lc[off_l].v>>1].len;
					assert(off_l < gc->off + gc->cnt);
					pd += (int32_t)p->x + 1;
				} else pd = (int32_t)p->x - (int32_t)q->x;
				sc[j - 1] = pd == qd && c == 0? -opt->match_pen : pd > qd? (int64_t)(c + (pd - qd) * a_dens + .499) : c;
			}

			ss = mss_find_all(0, gc->n_anchor - 1, sc, 10, 0, &n_ss);
			for (j = 0; j < n_ss; ++j) {
			}
			kfree(0, ss);
		}
	}
	kfree(km, sc);

	kfree(km, scnt);
	kfree(km, soff);
	kfree(km, intv);
	return h;
}
