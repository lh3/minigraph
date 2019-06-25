#include "mgpriv.h"
#include "bseq.h"
#include "ksort.h"
#include "mss.h"

/*********************************
 * Overlap query (from cgranges) *
 *********************************/

typedef struct {
	uint32_t st, en:31, rev:1;
	int32_t far;
	int32_t t, i, j;
} gg_intv_t;

#define sort_key_intv(a) ((a).st)
KRADIX_SORT_INIT(gg_intv, gg_intv_t, sort_key_intv, 4)

static int32_t gg_intv_index(int32_t n, gg_intv_t *a)
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

typedef struct {
	int64_t x;
	int32_t k, w;
} istack_t;

static int32_t gg_intv_overlap(void *km, int32_t n_a, const gg_intv_t *a, int32_t st, int32_t en, int32_t **b_, int32_t *m_b_)
{
	int32_t t = 0, h, *b = *b_, m_b = *m_b_, n = 0;
	istack_t stack[64], *p;

	for (h = 0; 1<<h <= n_a; ++h);
	--h;
	p = &stack[t++];
	p->k = h, p->x = (1LL<<p->k) - 1, p->w = 0; // push the root into the stack
	while (t) { // stack is not empyt
		istack_t z = stack[--t];
		if (z.k <= 3) { // the subtree is no larger than (1<<(z.k+1))-1; do a linear scan
			int32_t i, i0 = z.x >> z.k << z.k, i1 = i0 + (1LL<<(z.k+1)) - 1;
			if (i1 >= n_a) i1 = n_a;
			for (i = i0; i < i1 && a[i].st < en; ++i)
				if (st < a[i].en) {
					if (n == m_b) KEXPAND(km, b, m_b);
					b[n++] = i;
				}
		} else if (z.w == 0) { // if left child not processed
			int32_t y = z.x - (1LL<<(z.k-1));
			p = &stack[t++];
			p->k = z.k, p->x = z.x, p->w = 1;
			if (y >= n_a || a[y].far > st) {
				p = &stack[t++];
				p->k = z.k - 1, p->x = y, p->w = 0; // push the left child to the stack
			}
		} else if (z.x < n_a && a[z.x].st < en) {
			if (st < a[z.x].en) { // then z.x overlaps the query; write to the output array
				if (n == m_b) KEXPAND(km, b, m_b);
				b[n++] = z.x;
			}
			p = &stack[t++];
			p->k = z.k - 1, p->x = z.x + (1LL<<(z.k-1)), p->w = 0; // push the right child
		}
	}
	*b_ = b, *m_b_ = m_b;
	return n;
}

/**********************
 * Graph augmentation *
 **********************/

void mg_ggsimple(void *km, const mg_ggopt_t *opt, gfa_t *g, int32_t n_seq, const mg_bseq1_t *seq, mg_gchains_t *const* gcs)
{
	int32_t t, i, j, *scnt, *soff, max_acnt, *sc, m_ovlp = 0, *ovlp = 0, n_ins, m_ins;
	int64_t sum_acnt, sum_alen;
	uint64_t *meta;
	gg_intv_t *intv;
	double a_dens;
	gfa_ins_t *ins;

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
		return;
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
					rs = (int32_t)q->x + 1 - (int32_t)(q->y>>32&0xff);
					q = &gt->a[lc->off + lc->cnt - 1];
					re = (int32_t)q->x;
					assert(rs >= 0 && re > rs && re < g->seg[lc->v>>1].len);
					sum_alen += re - rs, sum_acnt += (q->x>>32) - (gt->a[lc->off].x>>32) + 1;
					if (lc->v&1)
						tmp = rs, rs = g->seg[lc->v>>1].len - re, re = g->seg[lc->v>>1].len - tmp;
				} else rs = 0, re = g->seg[lc->v>>1].len;
				// save the interval
				p = &intv[soff[lc->v>>1] + scnt[lc->v>>1]];
				++scnt[lc->v>>1];
				p->st = rs, p->en = re, p->rev = lc->v&1, p->far = -1;
				p->t = t, p->i = i, p->j = gc->off + j;
			}
		}
	}
	a_dens = (double)sum_acnt / sum_alen;

	// sort and index intervals
	for (i = 0; i < g->n_seg; ++i) {
		assert(soff[i+1] - soff[i] == scnt[i]);
		radix_sort_gg_intv(&intv[soff[i]], &intv[soff[i+1]]);
		gg_intv_index(soff[i+1] - soff[i], &intv[soff[i]]);
	}
	kfree(km, scnt);

	// extract poorly regions
	m_ins = n_ins = 0, ins = 0;
	KMALLOC(km, sc, max_acnt);
	KMALLOC(km, meta, max_acnt);
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
				int32_t s, off_l0 = off_l, pd, qd = (int32_t)p->y - (int32_t)q->y, c = (int32_t)(p->x>>32) - (int32_t)(q->x>>32) - 1;
				if (off_a == lc->off + lc->cnt) { // we are at the end of the current lchain
					pd = g->seg[lc->v>>1].len - (int32_t)q->x - 1;
					for (++off_l; off_l < gc->off + gc->cnt && gt->lc[off_l].cnt == 0; ++off_l)
						pd += g->seg[gt->lc[off_l].v>>1].len;
					assert(off_l < gc->off + gc->cnt);
					pd += (int32_t)p->x + 1;
				} else pd = (int32_t)p->x - (int32_t)q->x;
				if (pd == qd && c == 0) s = -opt->match_pen;
				else if (pd > qd) s = (int32_t)(c + (pd - qd) * a_dens + .499);
				else s = c; // TODO: check if this is an underestimate when there are overlaps on the query; perhaps the line above has addressed this.
				sc[j - 1] = s;
				meta[j-1] = (uint64_t)pd << 32 | off_l0;
			}

			// get regions to insert
			ss = mss_find_all(0, gc->n_anchor - 1, sc, 10, 0, &n_ss);
			off_a = gt->lc[gc->off].off;
			for (j = 0; j < n_ss; ++j) {
				const mg128_t *p, *q;
				int32_t st, en, ls, le, span, pd, k;
				gfa_ins_t I;

				// find the initial positions
				if (ss[j].st <= 1 || ss[j].en >= gt->n_a) continue; // not at the ends
				st = ss[j].st - 1, en = ss[j].en;
				q = &gt->a[off_a + st];
				p = &gt->a[off_a + en];
				span = p->y>>32&0xff;
				I.ctg = t;
				ls = (int32_t)meta[st], le = (int32_t)meta[en]; // first and last lchain
				I.v[0] = gt->lc[ls].v;
				I.v[1] = gt->lc[le].v;
				I.voff[0] = (int32_t)q->x + 1;
				I.voff[1] = (int32_t)p->x + 1 - span;
				I.coff[0] = (int32_t)q->y + 1;
				I.coff[1] = (int32_t)p->y + 1 - span;
				for (k = st, pd = 0; k < en; ++k) pd += meta[k] >> 32;
				pd -= span;

				// adjust for overlapping poistions
				if (I.coff[0] > I.coff[1]) {
					I.voff[1] += I.coff[0] - I.coff[1];
					pd += I.coff[0] - I.coff[1];
					I.coff[1] = I.coff[0];
				}
				if (I.v[0] == I.v[1] && I.voff[0] > I.voff[1]) {
					I.coff[1] += I.voff[0] - I.voff[1];
					pd += I.voff[0] - I.voff[1];
					I.voff[1] = I.voff[0];
				}
				pd -= gfa_ins_adj(g, 3, &I, seq[t].seq); // "3" is not used for now

				// filtering
				if (I.coff[1] - I.coff[0] < opt->min_var_len && pd < opt->min_var_len)
					continue;
				for (k = I.coff[0]; k < I.coff[1]; ++k) { // test ambiguous bases
					int c = seq[t].seq[k];
					if (c == 'n' || c == 'N') break;
				}
				if (k != I.coff[1]) continue; // no ambiguous bases on the insert
				for (k = ls; k <= le; ++k) { // find other mappings overlapping with the insert on the graph
					uint32_t v = gt->lc[k].v, len = g->seg[v>>1].len;
					int32_t s = 0, e = len, tmp, n_ovlp;
					if (k == ls && k == le) {
						s = (int32_t)gt->a[off_a+st].x + 1 - (int32_t)(gt->a[off_a+st].y>>32&0xff);
						e = (int32_t)gt->a[off_a+en].x + 1;
					} else if (k == ls) {
						s = (int32_t)gt->a[off_a+st].x + 1 - (int32_t)(gt->a[off_a+st].y>>32&0xff);
					} else if (k == le) {
						e = (int32_t)gt->a[off_a+en].x + 1;
					}
					if (v&1) tmp = s, s = len - e, e = len - tmp;
					n_ovlp = gg_intv_overlap(km, soff[(v>>1)+1] - soff[v>>1], &intv[soff[v>>1]], s, e, &ovlp, &m_ovlp);
					assert(n_ovlp > 0);
					if (n_ovlp > 1) continue;
				}
				//fprintf(stderr, "IN\t[%u:%d,%u:%d|%d] <=> [%d,%d|%d]\n", I.v[0], I.voff[0], I.v[1], I.voff[1], pd, I.coff[0], I.coff[1], I.coff[1] - I.coff[0]);
				if (n_ins == m_ins) KEXPAND(km, ins, m_ins);
				ins[n_ins++] = I;
			}
			kfree(0, ss);
		}
	}
	kfree(km, ovlp);
	kfree(km, sc);
	kfree(km, meta);
	kfree(km, soff);
	kfree(km, intv);

	if (n_ins > 0) {
		char **names, **seqs;
		KMALLOC(km, names, n_seq);
		KMALLOC(km, seqs, n_seq);
		for (i = 0; i < n_seq; ++i)
			names[i] = seq[i].name, seqs[i] = seq[i].seq;
		gfa_augment(g, n_ins, ins, n_seq, (const char*const*)names, (const char*const*)seqs);
		kfree(km, ins);
		kfree(km, names);
		kfree(km, seqs);
	}
}
