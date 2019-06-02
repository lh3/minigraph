#include "mgpriv.h"
#include "ksort.h"

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

int32_t mg_gchain1(void *km, const gfa_t *g, int32_t n_frag, mg_gfrag_t *frag, int32_t qlen, int32_t max_dist_g, int32_t max_dist_q, int32_t bw, uint64_t **u_)
{
	int32_t i, j, k, m_dst, n_dst, n_ext, j_st, n_u, n_v;
	int32_t *f, *p, *v, *t;
	uint64_t *u;
	gfa_path_dst_t *dst;
	gc_frag_t *a;
	mg_gfrag_t *swap;

	*u_ = 0;
	if (n_frag == 0) return 0;

	a = KMALLOC(km, gc_frag_t, n_frag);
	for (i = n_ext = 0; i < n_frag; ++i) { // a[] is a view of frag[]; for sorting
		mg_gfrag_t *r = &frag[i];
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
		u = KMALLOC(km, uint64_t, n_frag);
		for (i = 0; i < n_frag; ++i)
			u[i] = (uint64_t)frag[i].sc_chain<<32 | 1;
		*u_ = u;
		return n_frag;
	}
	radix_sort_gc(a, a + n_frag);

	v = KMALLOC(km, int32_t, n_ext);
	f = KMALLOC(km, int32_t, n_ext);
	p = KMALLOC(km, int32_t, n_ext);
	t = KCALLOC(km, int32_t, n_ext);

	m_dst = n_dst = 0, dst = 0;
	for (i = 0, j_st = 0; i < n_ext; ++i) { // core loop
		gc_frag_t *ai = &a[i];
		mg_gfrag_t *fi = &frag[ai->i];
		int32_t max_f = fi->sc_chain;
		int32_t max_j = -1;
		int32_t x = fi->qs + bw;
		while (j_st < i && a[j_st].srt + max_dist_q < fi->qs) ++j_st;
		if (x > qlen) x = qlen;
		x = j_st + find_max(i - j_st, a + j_st, x);
		n_dst = 0;
		for (j = x; j >= j_st; --j) { // collect potential destination vertices
			gc_frag_t *aj = &a[j];
			mg_gfrag_t *fj = &frag[aj->i];
			gfa_path_dst_t *q;
			int32_t min_dist = fi->rs + (g->seg[fj->v>>1].len - fj->re);
			if (min_dist > max_dist_g) continue;
			if (min_dist - bw > fi->qs - fj->qe) continue; // TODO: double check this line
			if (n_dst == m_dst) KEXPAND(km, dst, m_dst); // TODO: watch out the quadratic behavior!
			q = &dst[n_dst++];
			q->v = fj->v^1;
			q->meta = j;
			q->target_dist = (fi->qs - fj->qe) - min_dist + g->seg[fi->v>>1].len;
			if (q->target_dist < 0) q->target_dist = 0;
		}
		fprintf(stderr, "[%d] q_end=%d, src=%c%s, n_dst=%d, max_dist=%d\n", i, ai->srt, "><"[(fi->v&1)^1], g->seg[fi->v>>1].name, n_dst, max_dist_g + (g->seg[fi->v>>1].len - fi->rs));
		gfa_shortest_k(km, g, fi->v^1, n_dst, dst, max_dist_g + (g->seg[fi->v>>1].len - fi->rs), GFA_MAX_SHORT_K, 0);
		for (j = 0; j < n_dst; ++j) {
			gfa_path_dst_t *dj = &dst[j];
			int32_t gap, log_gap, sc;
			fprintf(stderr, "  [%d] dst=%c%s, n_path=%d, target=%d, opt_dist=%d\n", j, "><"[dj->v&1], g->seg[dj->v>>1].name, dj->n_path, dj->target_dist, dj->dist);
			if (dj->n_path == 0) continue;
			gap = dj->dist - dj->target_dist;
			if (gap < 0) gap = -gap;
			if (gap > bw) continue;
			log_gap = gap? mg_ilog2_32(gap) : 0;
			sc = fi->sc_chain;
			sc -= (int32_t)(gap * 0.2) + (log_gap >> 1);
			sc += f[dj->meta];
			if (sc > max_f) max_f = sc, max_j = dj->meta;
		}
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f;
	}

	u = mg_chain_backtrack(km, n_ext, f, p, v, t, 0, 0, n_frag - n_ext, &n_u, &n_v);
	kfree(km, f); kfree(km, p); kfree(km, t);
	assert(n_v == n_ext && n_u > 0);

	for (i = 0; i < n_frag - n_ext; ++i) {
		u[n_u++] = (uint64_t)frag[a[n_ext + i].i].sc_chain << 32 | 1;
		v[n_v++] = n_ext + i;
	}

	swap = KMALLOC(km, mg_gfrag_t, n_frag);
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			swap[k++] = frag[a[v[k0 + (ni - j - 1)]].i];
	}
	assert(k == n_frag);
	memcpy(frag, swap, n_frag * sizeof(mg_gfrag_t));
	*u_ = u;

	kfree(km, a);
	kfree(km, swap);
	kfree(km, v);
	return n_u;
}
