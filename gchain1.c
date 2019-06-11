#include <math.h>
#include <string.h>
#include "mgpriv.h"
#include "ksort.h" // for radix sort
#include "khash.h" // for __ac_Wang_hash()

typedef struct {
	uint32_t srt;
	int32_t i;
} gc_frag_t;

#define gc_frag_key(p) ((p).srt)
KRADIX_SORT_INIT(gc, gc_frag_t, gc_frag_key, 4)

static int32_t find_max(int32_t n, const gc_frag_t *gf, int32_t x)
{
	int32_t s = 0, e = n;
	if (n == 0) return -1;
	if (gf[n-1].srt < x) return n - 1;
	if (gf[0].srt >= x) return -1;
	while (e > s) { // TODO: finish this block
		int32_t m = s + (e - s) / 2;
		if (gf[m].srt >= x) e = m;
		else s = m + 1;
	}
	assert(s == e);
	return s;
}

static int32_t mg_target_dist(const gfa_t *g, const mg_lchain_t *l0, const mg_lchain_t *l1)
{
	return (l1->qs - l0->qe) - (g->seg[l0->v>>1].len - l0->re) + (g->seg[l1->v>>1].len - l1->rs);
}

int32_t mg_gchain1_dp(void *km, const gfa_t *g, int32_t n_lc, mg_lchain_t *lc, int32_t qlen, int32_t max_dist_g, int32_t max_dist_q, int32_t bw, uint64_t **u_)
{
	int32_t i, j, k, m_dst, n_dst, n_ext, j_st, n_u, n_v;
	int32_t *f, *p, *v, *t;
	uint64_t *u;
	gfa_path_dst_t *dst;
	gc_frag_t *a;
	mg_lchain_t *swap;

	*u_ = 0;
	if (n_lc == 0) return 0;

	a = KMALLOC(km, gc_frag_t, n_lc);
	for (i = n_ext = 0; i < n_lc; ++i) { // a[] is a view of frag[]; for sorting
		mg_lchain_t *r = &lc[i];
		gc_frag_t *ai = &a[i];
		int32_t is_isolated = 0;
		if (r->rs > max_dist_g && g->seg[r->v>>1].len - r->re > max_dist_g)
			is_isolated = 1;
		ai->srt = (uint32_t)is_isolated<<31 | r->qe;
		ai->i = i;
		if (!is_isolated) ++n_ext;
	}
	if (n_ext < 2) { // no graph chaining needed; early return
		kfree(km, a);
		u = KMALLOC(km, uint64_t, n_lc);
		for (i = 0; i < n_lc; ++i)
			u[i] = (uint64_t)lc[i].score<<32 | 1;
		*u_ = u;
		return n_lc;
	}
	radix_sort_gc(a, a + n_lc);

	v = KMALLOC(km, int32_t, n_ext);
	f = KMALLOC(km, int32_t, n_ext);
	p = KMALLOC(km, int32_t, n_ext);
	t = KCALLOC(km, int32_t, n_ext);

	m_dst = n_dst = 0, dst = 0;
	for (i = 0, j_st = 0; i < n_ext; ++i) { // core loop
		gc_frag_t *ai = &a[i];
		mg_lchain_t *li = &lc[ai->i];
		int32_t max_f = li->score;
		int32_t max_j = -1;
		int32_t x = li->qs + bw;
		while (j_st < i && a[j_st].srt + max_dist_q < li->qs) ++j_st;
		if (x > qlen) x = qlen;
		x = j_st + find_max(i - j_st, a + j_st, x);
		n_dst = 0;
		for (j = x; j >= j_st; --j) { // collect potential destination vertices
			gc_frag_t *aj = &a[j];
			mg_lchain_t *lj = &lc[aj->i];
			gfa_path_dst_t *q;
			int32_t min_dist = li->rs + (g->seg[lj->v>>1].len - lj->re);
			if (min_dist > max_dist_g) continue;
			if (min_dist - bw > li->qs - lj->qe) continue; // TODO: double check this line
			if (n_dst == m_dst) KEXPAND(km, dst, m_dst); // TODO: watch out the quadratic behavior!
			q = &dst[n_dst++];
			q->v = lj->v^1;
			q->meta = j;
			q->target_dist = mg_target_dist(g, lj, li);
			if (q->target_dist < 0) q->target_dist = 0;
		}
		fprintf(stderr, "[src:%d] q_end=%d, src=%c%s[%d], n_dst=%d, max_dist=%d\n", i, ai->srt, "><"[(li->v&1)^1], g->seg[li->v>>1].name, li->v^1, n_dst, max_dist_g + (g->seg[li->v>>1].len - li->rs));
		gfa_shortest_k(km, g, li->v^1, n_dst, dst, max_dist_g + (g->seg[li->v>>1].len - li->rs), GFA_MAX_SHORT_K, 0);
		for (j = 0; j < n_dst; ++j) {
			gfa_path_dst_t *dj = &dst[j];
			int32_t gap, log_gap, sc;
			if (dj->n_path == 0) continue;
			gap = dj->dist - dj->target_dist;
			if (gap < 0) gap = -gap;
			if (gap > bw) continue;
			log_gap = gap? mg_ilog2_32(gap) : 0;
			sc = li->score;
			sc -= (int32_t)(gap * 0.2) + (log_gap >> 1);
			sc += f[dj->meta];
			fprintf(stderr, "  [dst:%d] dst=%c%s[%d], n_path=%d, target=%d, opt_dist=%d, score=%d\n", j, "><"[dj->v&1], g->seg[dj->v>>1].name, dj->v, dj->n_path, dj->target_dist, dj->dist, sc);
			if (sc > max_f) max_f = sc, max_j = dj->meta;
		}
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f;
	}
	kfree(km, dst);

	u = mg_chain_backtrack(km, n_ext, f, p, v, t, 0, 0, n_lc - n_ext, &n_u, &n_v);
	kfree(km, f); kfree(km, p); kfree(km, t);
	assert(n_v == n_ext && n_u > 0);

	for (i = 0; i < n_lc - n_ext; ++i) {
		u[n_u++] = (uint64_t)lc[a[n_ext + i].i].score << 32 | 1;
		v[n_v++] = n_ext + i;
	}

	swap = KMALLOC(km, mg_lchain_t, n_lc);
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			swap[k++] = lc[a[v[k0 + (ni - j - 1)]].i];
	}
	assert(k == n_lc);
	memcpy(lc, swap, n_lc * sizeof(mg_lchain_t));
	*u_ = u;

	kfree(km, a);
	kfree(km, swap);
	kfree(km, v);
	return n_u;
}

void mg_gchain_extra(const gfa_t *g, mg_gchains_t *gs)
{
	int32_t i, j, k;
	for (i = 0; i < gs->n_gc; ++i) { // iterate over gchains
		mg_gchain_t *p = &gs->gc[i];
		const mg_llchain_t *q;
		const mg128_t *last_a;
		int32_t q_span, rest_pl, tmp;

		p->qs = p->qe = p->ps = p->pe = -1, p->plen = p->blen = p->mlen = 0, p->div = -1.0f;
		if (p->cnt == 0) continue;

		assert(gs->lc[p->off].cnt > 0 && gs->lc[p->off + p->cnt - 1].cnt > 0); // first and last lchains can't be empty
		q = &gs->lc[p->off];
		q_span = (int32_t)(gs->a[q->off].y>>32&0xff);
		p->qs = (int32_t)gs->a[q->off].y + 1 - q_span;
		p->ps = (int32_t)gs->a[q->off].x + 1 - q_span;
		tmp = (int32_t)(gs->a[q->off].x>>32);
		assert(p->qs >= 0 && p->ps >= 0);
		q = &gs->lc[p->off + p->cnt - 1];
		p->qe = (int32_t)gs->a[q->off + q->cnt - 1].y + 1;
		p->pe = g->seg[q->v>>1].len - (int32_t)gs->a[q->off + q->cnt - 1].x - 1; // this is temporary
		tmp = (int32_t)(gs->a[q->off + q->cnt - 1].x>>32) - tmp + 1;
		assert(p->n_anchor > 0 && tmp >= p->n_anchor);
		p->div = log((double)tmp / p->n_anchor) / q_span;

		rest_pl = 0; // this value is never used if the first lchain is not empty (which should always be true)
		last_a = &gs->a[gs->lc[p->off].off];
		for (j = 0; j < p->cnt; ++j) { // iterate over lchains
			const mg_llchain_t *q = &gs->lc[p->off + j];
			int32_t vlen = g->seg[q->v>>1].len;
			p->plen += vlen;
			for (k = 0; k < q->cnt; ++k) { // iterate over anchors
				const mg128_t *r = &gs->a[q->off + k];
				int32_t pl, ql = (int32_t)r->y - (int32_t)last_a->y;
				int32_t span = (int32_t)(r->y>>32&0xff);
				if (j == 0 && k == 0) { // the first anchor on the first lchain
					pl = ql = span;
				} else if (j > 0 && k == 0) { // the first anchor but not on the first lchain
					pl = (int32_t)r->x + 1 + rest_pl;
				} else {
					pl = (int32_t)r->x - (int32_t)last_a->x;
				}
				p->blen += pl > ql? pl : ql;
				p->mlen += pl > span && ql > span? span : pl < ql? pl : ql;
				last_a = r;
			}
			if (q->cnt == 0) rest_pl += vlen;
			else rest_pl = vlen - (int32_t)gs->a[q->off + q->cnt - 1].x - 1;
		}
		p->pe = p->plen - p->pe;
		assert(p->pe >= p->ps);
	}
}

void mg_gchain_sort_by_score(void *km, mg_gchains_t *gcs)
{
	mg128_t *z;
	mg_gchain_t *gc;
	int32_t i;
	z = KMALLOC(km, mg128_t, gcs->n_gc);
	gc = KMALLOC(km, mg_gchain_t, gcs->n_gc);
	for (i = 0; i < gcs->n_gc; ++i)
		z[i].x = (uint64_t)gcs->gc[i].score << 32 | gcs->gc[i].hash, z[i].y = i;
	radix_sort_128x(z, z + gcs->n_gc);
	for (i = gcs->n_gc - 1; i >= 0; --i)
		gc[gcs->n_gc - 1 - i] = gcs->gc[z[i].y];
	memcpy(gcs->gc, gc, gcs->n_gc * sizeof(mg_gchain_t));
	kfree(km, z); kfree(km, gc);
}

static inline void copy_lchain(mg_llchain_t *q, const mg_lchain_t *p, int32_t *n_a, mg128_t *a_new, const mg128_t *a_old)
{
	q->cnt = p->cnt, q->v = p->v, q->score = p->score;
	memcpy(&a_new[*n_a], &a_old[p->off], p->cnt * sizeof(mg128_t));
	q->off = *n_a;
	(*n_a) += p->cnt;
}

// TODO: if frequent malloc() is a concern, filter first and then generate gchains; or generate gchains in thread-local pool and then move to global malloc()
mg_gchains_t *mg_gchain_gen(void *km_dst, void *km, const gfa_t *g, int32_t n_u, const uint64_t *u, const mg_lchain_t *lc, const mg128_t *a, uint32_t hash, int32_t min_gc_cnt, int32_t min_gc_score)
{
	mg_gchains_t *gc = 0;
	mg_llchain_t *tmp;
	int32_t i, j, k, st, n_g, n_a, s_tmp, n_tmp, m_tmp;

	// count the number of gchains and remaining anchors
	for (i = 0, st = 0, n_g = n_a = 0; i < n_u; ++i) {
		int32_t m = 0, nui = (int32_t)u[i];
		for (j = 0; j < nui; ++j) m += lc[st + j].cnt; // m is the number of anchors in this gchain
		if (m >= min_gc_cnt && u[i]>>32 >= min_gc_score)
			++n_g, n_a += m;
		st += nui;
	}
	if (n_g == 0) return 0;

	// preallocate
	gc = KCALLOC(km_dst, mg_gchains_t, 1);
	gc->km = km_dst;
	gc->n_gc = n_g, gc->n_a = n_a;
	gc->gc = KCALLOC(km_dst, mg_gchain_t, n_g);
	gc->a = KMALLOC(km_dst, mg128_t, n_a);

	// core loop
	tmp = 0; s_tmp = n_tmp = m_tmp = 0;
	for (i = k = 0, st = 0, n_a = 0; i < n_u; ++i) {
		int32_t m = 0, nui = (int32_t)u[i];
		for (j = 0; j < nui; ++j) m += lc[st + j].cnt;
		if (m >= min_gc_cnt && u[i]>>32 >= min_gc_score) {
			mg_llchain_t *q;
			uint32_t h = hash;

			gc->gc[k].n_anchor = m;
			gc->gc[k].score = u[i]>>32;
			gc->gc[k].off = s_tmp;

			for (j = 0; j < nui; ++j) {
				const mg_lchain_t *p = &lc[st + j];
				h += __ac_Wang_hash(p->qs) + __ac_Wang_hash(p->re) + __ac_Wang_hash(p->v);
			}
			gc->gc[k].hash = __ac_Wang_hash(h);

			if (n_tmp == m_tmp) KEXPAND(km, tmp, m_tmp);
			copy_lchain(&tmp[n_tmp++], &lc[st], &n_a, gc->a, a); // copy the first lchain

			for (j = 1; j < nui; ++j) {
				const mg_lchain_t *l0 = &lc[st + j - 1], *l1 = &lc[st + j];
				gfa_path_dst_t dst;
				int32_t s, n_pathv;
				gfa_pathv_t *p;
				dst.v = l0->v ^ 1;
				dst.target_dist = mg_target_dist(g, l0, l1);
				p = gfa_shortest_k(km, g, l1->v^1, 1, &dst, dst.target_dist, GFA_MAX_SHORT_K, &n_pathv);
				fprintf(stderr, "%c%s[%d] -> %c%s[%d], dist=%d, target=%d\n", "><"[(l1->v^1)&1], g->seg[l1->v>>1].name, l1->v^1, "><"[(l0->v^1)&1], g->seg[l0->v>>1].name, l0->v^1, dst.dist, dst.target_dist);
				assert(n_pathv > 0);
				for (s = 1; s < n_pathv - 1; ++s) {
					if (n_tmp == m_tmp) KEXPAND(km, tmp, m_tmp);
					q = &tmp[n_tmp++];
					q->off = q->cnt = q->score = 0;
					q->v = p[s].v;
				}
				if (n_tmp == m_tmp) KEXPAND(km, tmp, m_tmp);
				copy_lchain(&tmp[n_tmp++], l1, &n_a, gc->a, a);
				kfree(km, p);
			}
			gc->gc[k].cnt = n_tmp - s_tmp;
			++k, s_tmp = n_tmp;
		}
		st += nui;
	}
	assert(n_a == gc->n_a);

	gc->n_lc = n_tmp;
	gc->lc = KMALLOC(km_dst, mg_llchain_t, n_tmp);
	memcpy(gc->lc, tmp, n_tmp * sizeof(mg_llchain_t));
	kfree(km, tmp);

	mg_gchain_extra(g, gc);
	mg_gchain_sort_by_score(km, gc);
	return gc;
}

void mg_gchain_free(mg_gchains_t *gs)
{
	void *km;
	if (gs == 0) return;
	km = gs->km;
	kfree(km, gs->gc); kfree(km, gs->a); kfree(km, gs->lc); kfree(km, gs);
}
