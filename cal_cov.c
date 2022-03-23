#include <assert.h>
#include <string.h>
#include "mgpriv.h"
#include "gfa-priv.h"
#include "algo.h"
#include "kalloc.h"

void mg_cov_map(const gfa_t *g, const mg_gchains_t *gt, int32_t min_mapq, int32_t min_blen, double *c_seg, double *c_link, const char *qname)
{
	int32_t i, j;
	if (c_seg == 0 && c_link == 0) return;
	if (gt == 0 || gt->n_gc == 0) return;
	for (i = 0; i < gt->n_gc; ++i) {
		const mg_gchain_t *gc = &gt->gc[i];
		const mg128_t *last_an;
		assert(gc->cnt > 0 && gc->n_anchor > 0);
		if (gc->mapq < min_mapq || gc->blen < min_blen) continue;
		// count segment coverage
		for (j = 0; j < gc->cnt; ++j) {
			const mg_llchain_t *lc = &gt->lc[gc->off + j];
			int32_t s, e;
			s = 0, e = g->seg[lc->v>>1].len;
			if (j == 0) s = (int32_t)gt->a[lc->off].x + 1 - (int32_t)(gt->a[lc->off].y>>32&0xff);
			if (j == gc->cnt - 1) e = (int32_t)gt->a[lc->off + lc->cnt - 1].x + 1;
			if (c_seg) c_seg[lc->v>>1] += (double)(e - s) / g->seg[lc->v>>1].len;
		}
		// count link
		assert(gt->lc[gc->off].cnt > 0);
		last_an = &gt->a[gt->lc[gc->off].off + gt->lc[gc->off].cnt - 1];
		for (j = 1; j < gc->cnt; ++j) {
			const mg_llchain_t *lc0 = &gt->lc[gc->off + j - 1];
			const mg_llchain_t *lc1 = &gt->lc[gc->off + j];
			int64_t a01, a10;
			if (lc1->cnt > 0) {
				const mg128_t *curr_an = &gt->a[lc1->off];
				int32_t is_skip = (mg_seg_id(*curr_an) != mg_seg_id(*last_an));
				last_an = &gt->a[lc1->off + lc1->cnt - 1];
				if (is_skip) continue;
			}
			a01 = gfa_find_arc(g, lc0->v, lc1->v);
			a10 = gfa_find_arc(g, lc1->v^1, lc0->v^1);
			if (a01 < 0 || a10 < 0) {
				if (mg_verbose >= 2)
					fprintf(stderr, "[W] Multi/disconnected link: %c%s[%d] -> %c%s[%d] (%s, %ld, %ld). Continue anyway!\n",
							"><"[lc0->v&1], g->seg[lc0->v>>1].name, lc0->v,
							"><"[lc1->v&1], g->seg[lc1->v>>1].name, lc1->v, qname, (long)a01, (long)a10);
				continue;
			}
			assert((g->arc[a01].comp ^ g->arc[a10].comp) == 1);
			if (c_link) c_link[a01] += 1.0, c_link[a10] += 1.0;
		}
	}
}

void mg_cov_asm(const gfa_t *g, int32_t n_seq, mg_gchains_t *const *gcs, int32_t min_mapq, int32_t min_blen, double *cov_seg, double *cov_link)
{
	int32_t i, j, t, *soff, *scnt, *cnt_link;
	int64_t k;
	mg_intv_t *sintv = 0;
	void *km = 0;

	// precalculate the size of sintv[] for each segment
	KCALLOC(km, scnt, g->n_seg);
	for (t = 0; t < n_seq; ++t) {
		const mg_gchains_t *gt = gcs[t];
		if (gt == 0 || gt->n_gc == 0) continue;
		for (i = 0; i < gt->n_gc; ++i) {
			const mg_gchain_t *gc = &gt->gc[i];
			assert(gc->cnt > 0 && gc->n_anchor > 0);
			if (gc->mapq < min_mapq || gc->blen < min_blen) continue;
			for (j = 0; j < gc->cnt; ++j) {
				const mg_llchain_t *lc = &gt->lc[gc->off + j];
				++scnt[lc->v>>1];
			}
		}
	}
	KMALLOC(km, soff, g->n_seg + 1);
	for (soff[0] = 0, i = 1; i <= g->n_seg; ++i)
		soff[i] = soff[i - 1] + scnt[i - 1];
	memset(scnt, 0, 4 * g->n_seg);
	KMALLOC(km, sintv, soff[g->n_seg]);

	// fill sintv[]
	KCALLOC(km, cnt_link, g->n_arc);
	for (t = 0; t < n_seq; ++t) {
		const mg_gchains_t *gt = gcs[t];
		if (gt == 0 || gt->n_gc == 0) continue;
		for (i = 0; i < gt->n_gc;) {
			const mg_gchain_t *gc = &gt->gc[i];
			if (gc->mapq < min_mapq || gc->blen < min_blen) continue;
			// count segment coverage
			for (j = 0; j < gc->cnt; ++j) {
				const mg_llchain_t *lc = &gt->lc[gc->off + j];
				int32_t s, e, tmp;
				mg_intv_t *p;
				s = 0, e = g->seg[lc->v>>1].len;
				if (j == 0) s = (int32_t)gt->a[lc->off].x + 1 - (int32_t)(gt->a[lc->off].y>>32&0xff);
				if (j == gc->cnt - 1) e = (int32_t)gt->a[lc->off + lc->cnt - 1].x + 1;
				if (lc->v&1) // convert to the forward strand of segment lc->v>>1
					tmp = g->seg[lc->v>>1].len - s, s = g->seg[lc->v>>1].len - e, e = tmp;
				p = &sintv[soff[lc->v>>1] + scnt[lc->v>>1]];
				++scnt[lc->v>>1];
				p->st = s, p->en = e, p->rev = lc->v&1, p->far = -1, p->i = -1;
			}
			// count link
			for (j = 1; j < gc->cnt; ++j) {
				const mg_llchain_t *lc0 = &gt->lc[gc->off + j - 1];
				const mg_llchain_t *lc1 = &gt->lc[gc->off + j];
				int64_t a01, a10;
				a01 = gfa_find_arc(g, lc0->v, lc1->v);
				a10 = gfa_find_arc(g, lc1->v^1, lc0->v^1);
				assert(a01 >= 0 && a10 >= 0);
				assert((g->arc[a01].comp ^ g->arc[a10].comp) == 1);
				++cnt_link[a01];
				++cnt_link[a10];
			}
		}
	}

	// update cov_link[] and cov_seg[]
	for (k = 0; k < g->n_arc; ++k)
		if (cnt_link[k] > 0) cov_link[k] += 1.0;
	for (i = 0; i < g->n_seg; ++i) {
		int32_t st = 0, en = 0, cov = 0;
		assert(scnt[i] == soff[i+1] - soff[i]);
		radix_sort_mg_intv(&sintv[soff[i]], &sintv[soff[i+1]]);
		for (j = soff[i]; j < soff[i+1]; ++j) {
			if (sintv[j].st > en)
				cov += en - st, st = sintv[j].st, en = sintv[j].en;
			else en = sintv[j].en > en? sintv[j].en : en;
		}
		cov += en - st;
		cov_seg[i] += (double)cov / g->seg[i].len;
	}

	// free
	kfree(km, cnt_link);
	kfree(km, sintv); kfree(km, soff); kfree(km, scnt);
}
