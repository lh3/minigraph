#include <assert.h>
#include "mgpriv.h"
#include "gfa-priv.h"

void mg_count_cov_simple(const gfa_t *g, const mg_gchains_t *gt, int32_t min_mapq, int32_t min_blen, int64_t *c_seg, int32_t *c_link, const char *qname)
{
	int32_t i, j, a_off = 0;
	if (c_seg == 0 && c_link == 0) return;
	if (gt == 0 || gt->n_gc == 0) return;
	for (i = 0; i < gt->n_gc; a_off += gt->gc[i++].n_anchor) {
		const mg_gchain_t *gc = &gt->gc[i];
		assert(gc->cnt > 0 && gc->n_anchor > 0);
		if (gc->mapq < min_mapq || gc->blen < min_blen) continue;
		// count segment coverage
		for (j = 0; j < gc->cnt; ++j) {
			const mg_llchain_t *lc = &gt->lc[gc->off + j];
			int32_t s, e;
			s = 0, e = g->seg[lc->v>>1].len;
			if (j == 0) s = (int32_t)gt->a[lc->off].x + 1 - (int32_t)(gt->a[lc->off].y>>32&0xff);
			if (j == gc->cnt - 1) e = (int32_t)gt->a[lc->off + lc->cnt - 1].x + 1;
			if (c_seg) c_seg[lc->v>>1] += e - s;
		}
		// count link
		for (j = 1; j < gc->cnt; ++j) {
			const mg_llchain_t *lc0 = &gt->lc[gc->off + j - 1];
			const mg_llchain_t *lc1 = &gt->lc[gc->off + j];
			int64_t a01, a10;
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
			if (c_link) ++c_link[g->arc[a01].comp == 0? a01 : a10];
		}
	}
}
