#include <math.h>
#include <string.h>
#include "mgpriv.h"
#include "ksort.h" // for radix sort
#include "khashl.h" // for kh_hash_uint32()
#include "gfa-priv.h"

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
	// below equals (l1->qs - l0->qe) - min_dist + g->seg[l1->v>>1].len; see mg_gchain1_dp() for the calculation of min_dist
	return (l1->qs - l0->qe) - (g->seg[l0->v>>1].len - l0->re) + (g->seg[l1->v>>1].len - l1->rs);
	// when l0->v == l1->v, the above becomes (l1->qs - l0->qe) - (l1->rs - l0->re), which is what we want
}

static inline int32_t cal_sc(const mg_path_dst_t *dj, const mg_lchain_t *li, const mg_lchain_t *lc, const mg128_t *an, const gc_frag_t *a, const int32_t *f,
							 int bw, int ref_bonus, float chn_pen_gap)
{
	const mg_lchain_t *lj;
	int32_t gap, sc, segi, segj;
	float lin_pen, log_pen;
	if (dj->n_path == 0) return INT32_MIN;
	segi = (an[li->off].y & MG_SEED_SEG_MASK) >> MG_SEED_SEG_SHIFT;
	gap = dj->dist - dj->target_dist;
	lj = &lc[a[dj->meta].i];
	segj = (an[lj->off + lj->cnt - 1].y & MG_SEED_SEG_MASK) >> MG_SEED_SEG_SHIFT;
	if (gap < 0) gap = -gap;
	if (segi == segj && gap > bw) return INT32_MIN;
	if (lj->qe <= li->qs) sc = li->score;
	else sc = (int32_t)((double)(li->qe - lj->qe) / (li->qe - li->qs) * li->score + .499); // dealing with overlap on query
	//sc += dj->mlen; // TODO: is this line the right thing to do?
	if (dj->is_0) sc += ref_bonus;
	lin_pen = chn_pen_gap * (float)gap;
	log_pen = gap >= 2? mg_log2(gap) : 0.0f;
	sc -= (int32_t)(lin_pen + log_pen);
	sc += f[dj->meta];
	return sc;
}

int32_t mg_gchain1_dp(void *km, const gfa_t *g, int32_t *n_lc_, mg_lchain_t *lc, int32_t qlen, int32_t max_dist_g, int32_t max_dist_q, int32_t bw, int32_t max_skip,
					  int32_t ref_bonus, float chn_pen_gap, float chn_pen_skip, float mask_level, const mg128_t *an, uint64_t **u_)
{
	int32_t i, j, k, m_dst, n_dst, n_ext, n_u, n_v, n_lc = *n_lc_;
	int32_t *f, *v, *t;
	int64_t *p;
	uint64_t *u;
	mg_path_dst_t *dst;
	gc_frag_t *a;
	mg_lchain_t *swap;
	char *qs;

	*u_ = 0;
	if (n_lc == 0) return 0;

	KMALLOC(km, a, n_lc);
	for (i = n_ext = 0; i < n_lc; ++i) { // a[] is a view of frag[]; for sorting
		mg_lchain_t *r = &lc[i];
		gc_frag_t *ai = &a[i];
		int32_t is_isolated = 0, min_end_dist_g;
		r->dist_pre = -1;
		min_end_dist_g = g->seg[r->v>>1].len - r->re;
		if (r->rs < min_end_dist_g) min_end_dist_g = r->rs;
		if (min_end_dist_g > max_dist_g) is_isolated = 1; // if too far from segment ends
		else if (min_end_dist_g>>3 > r->score) is_isolated = 1; // if the lchain too small relative to distance to the segment ends
		ai->srt = (uint32_t)is_isolated<<31 | r->qe;
		ai->i = i;
		if (!is_isolated) ++n_ext;
	}
	if (n_ext < 2) { // no graph chaining needed; early return
		kfree(km, a);
		KMALLOC(km, u, n_lc);
		for (i = 0; i < n_lc; ++i)
			u[i] = (uint64_t)lc[i].score<<32 | 1;
		*u_ = u;
		return n_lc;
	}
	radix_sort_gc(a, a + n_lc);

	KMALLOC(km, v, n_lc);
	KMALLOC(km, f, n_ext);
	KMALLOC(km, p, n_ext);
	KCALLOC(km, t, n_ext);

	KMALLOC(km, qs, max_dist_q + 1);
	m_dst = n_dst = 0, dst = 0;
	for (i = 0; i < n_ext; ++i) { // core loop
		gc_frag_t *ai = &a[i];
		mg_lchain_t *li = &lc[ai->i];
		int32_t segi = (an[li->off].y & MG_SEED_SEG_MASK) >> MG_SEED_SEG_SHIFT;
		{ // collect end points potentially reachable from _i_
			int32_t x = li->qs + bw, n_skip = 0;
			if (x > qlen) x = qlen;
			x = find_max(i, a, x);
			n_dst = 0;
			for (j = x; j >= 0; --j) { // collect potential destination vertices
				gc_frag_t *aj = &a[j];
				mg_lchain_t *lj = &lc[aj->i];
				mg_path_dst_t *q;
				int32_t target_dist, segj, dq;
				if (lj->qs >= li->qs) continue; // lj is contained in li on the query coordinate
				if (lj->qe > li->qs) { // test overlap on the query
					int o = lj->qe - li->qs;
					if (o > (lj->qe - lj->qs) * mask_level || o > (li->qe - li->qs) * mask_level)
						continue;
				}
				dq = li->qs - lj->qe;
				segj = (an[lj->off + lj->cnt - 1].y & MG_SEED_SEG_MASK) >> MG_SEED_SEG_SHIFT;
				if (segi == segj) {
					if (dq > max_dist_q) break; // if query gap too large, stop
				} else {
					if (dq > max_dist_g && dq > max_dist_q) break;
				}
				if (li->v != lj->v) { // the two linear chains are on two different segments
					int32_t min_dist = li->rs + (g->seg[lj->v>>1].len - lj->re); // minimal graph gap
					if (min_dist > max_dist_g) continue; // graph gap too large
					if (segi == segj && min_dist - bw > li->qs - lj->qe) continue; // when li->qs < lj->qe, the condition turns to min_dist + (lj->qe - li->qs) > bw, which is desired
					target_dist = mg_target_dist(g, lj, li);
					if (target_dist < 0) continue; // this may happen if the query overlap is far too large
				} else if (lj->rs >= li->rs || lj->re >= li->re) { // not colinear
					continue;
				} else {
					int32_t dr = li->rs - lj->re, w = dr > dq? dr - dq : dq - dr;
					if (segi == segj && w > bw) continue; // test bandwidth
					if (dr > max_dist_g || dr < -max_dist_g) continue;
					if (lj->re > li->rs) { // test overlap on the graph segment
						int o = lj->re - li->rs;
						if (o > (lj->re - lj->rs) * mask_level || o > (li->re - li->rs) * mask_level)
							continue;
					}
					target_dist = mg_target_dist(g, lj, li);
				}
				if (n_dst == m_dst) KEXPAND(km, dst, m_dst); // TODO: watch out the quadratic behavior!
				q = &dst[n_dst++];
				memset(q, 0, sizeof(mg_path_dst_t));
				q->inner = (li->v == lj->v);
				q->v = lj->v^1;
				q->meta = j;
				q->qlen = li->qs - lj->qe;
				q->target_dist = target_dist;
				q->target_hash = 0;
				q->check_hash = 0;
				if (t[j] == i) {
					if (++n_skip > max_skip)
						break;
				}
				if (p[j] >= 0) t[p[j]] = i;
			}
		}
		{ // confirm reach-ability
			int32_t k;
			// test reach-ability without sequences
			mg_shortest_k(km, g, li->v^1, n_dst, dst, max_dist_g + (g->seg[li->v>>1].len - li->rs), MG_MAX_SHORT_K, 0);
			// remove unreachable destinations
			for (j = k = 0; j < n_dst; ++j) {
				mg_path_dst_t *dj = &dst[j];
				int32_t sc;
				if (dj->n_path == 0) continue; // not reachable
				sc = cal_sc(dj, li, lc, an, a, f, bw, ref_bonus, chn_pen_gap);
				if (sc == INT32_MIN) continue; // out of band
				if (sc + li->score < 0) continue; // negative score and too low
				dst[k++] = dst[j];
			}
			n_dst = k;
		}
		{ // DP
			int32_t max_f = li->score, max_j = -1, max_d = -1, max_inner = 0;
			uint32_t max_hash = 0;
			for (j = 0; j < n_dst; ++j) {
				mg_path_dst_t *dj = &dst[j];
				int32_t sc;
				sc = cal_sc(dj, li, lc, an, a, f, bw, ref_bonus, chn_pen_gap);
				if (sc == INT32_MIN) continue;
				if (mg_dbg_flag & MG_DBG_GC1) {
					mg_lchain_t *lj = &lc[a[dj->meta].i];
					fprintf(stderr, "  [dst:%d] dst=%c%s[%d], n_path=%d, target=%d, opt_dist=%d, score=%d, q_intv=[%d,%d), g_intv=[%d,%d)\n", dj->meta, "><"[dj->v&1], g->seg[dj->v>>1].name, dj->v, dj->n_path, dj->target_dist - g->seg[li->v>>1].len, dj->dist - g->seg[li->v>>1].len, sc, lj->qs, lj->qe, lj->rs, lj->re);
				}
				if (sc > max_f) max_f = sc, max_j = dj->meta, max_d = dj->dist, max_hash = dj->hash, max_inner = dj->inner;
			}
			f[i] = max_f, p[i] = max_j;
			li->dist_pre = max_d;
			li->hash_pre = max_hash;
			li->inner_pre = max_inner;
			v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f;
			if (mg_dbg_flag & MG_DBG_GC1) fprintf(stderr, " [opt:%d] opt=%d, max_f=%d\n", ai->i, max_j, max_f);
		}
	}
	kfree(km, dst);
	kfree(km, qs);
	if (mg_dbg_flag & MG_DBG_GC1) {
		int32_t mmax_f = 0, mmax_i = -1;
		for (i = 0; i < n_ext; ++i) if (f[i] > mmax_f) mmax_f = f[i], mmax_i = i;
		i = mmax_i; while (i >= 0) { fprintf(stderr, "[best] i=%d, seg=%s, max_f=%d, chn_pen_gap=%f\n", a[i].i, g->seg[lc[a[i].i].v>>1].name, f[i], chn_pen_gap); i = p[i]; }
	}

	u = mg_chain_backtrack(km, n_ext, f, p, v, t, 0, 0, INT32_MAX, n_lc - n_ext, &n_u, &n_v);
	kfree(km, f); kfree(km, p); kfree(km, t);

	for (i = 0; i < n_lc - n_ext; ++i) {
		u[n_u++] = (uint64_t)lc[a[n_ext + i].i].score << 32 | 1;
		v[n_v++] = n_ext + i;
	}

	KMALLOC(km, swap, n_v);
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			swap[k++] = lc[a[v[k0 + (ni - j - 1)]].i];
	}
	assert(k == n_v);
	memcpy(lc, swap, n_v * sizeof(mg_lchain_t));
	*n_lc_ = n_v;
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
		int32_t q_span, rest_pl, tmp, n_mini;

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
		n_mini = (int32_t)(gs->a[q->off + q->cnt - 1].x>>32) - tmp + 1;
		assert(p->n_anchor > 0);

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
				if (ql < 0) ql = -ql, n_mini += (int32_t)(last_a->x>>32) - (int32_t)(r->x>>32); // dealing with overlapping query at junctions
				p->blen += pl > ql? pl : ql;
				p->mlen += pl > span && ql > span? span : pl < ql? pl : ql;
				last_a = r;
			}
			if (q->cnt == 0) rest_pl += vlen;
			else rest_pl = vlen - (int32_t)gs->a[q->off + q->cnt - 1].x - 1;
		}
		p->pe = p->plen - p->pe;
		assert(p->pe >= p->ps);
		// here n_mini >= p->n_anchor should stand almost all the time
		p->div = n_mini >= p->n_anchor? log((double)n_mini / p->n_anchor) / q_span : log((double)p->n_anchor / n_mini) / q_span;
	}
}

/*
 * Generate graph chains
 */
typedef struct {
	void *km;
	const gfa_t *g;
	const gfa_edseq_t *es;
	const char *qseq;
	int32_t n_seg, n_llc, m_llc, n_a;
	mg_llchain_t *llc;
} bridge_aux_t;

static inline void copy_lchain(mg_llchain_t *q, const mg_lchain_t *p, int32_t *n_a, mg128_t *a_new, const mg128_t *a_old, int32_t ed)
{
	q->cnt = p->cnt, q->v = p->v, q->score = p->score, q->ed = ed;
	memcpy(&a_new[*n_a], &a_old[p->off], q->cnt * sizeof(mg128_t));
	q->off = *n_a;
	(*n_a) += q->cnt;
}

static int32_t bridge_shortk(bridge_aux_t *aux, const mg_lchain_t *l0, const mg_lchain_t *l1)
{
	int32_t s, n_pathv;
	mg_path_dst_t dst;
	mg_pathv_t *p;
	memset(&dst, 0, sizeof(mg_path_dst_t));
	dst.v = l0->v ^ 1;
	assert(l1->dist_pre >= 0);
	dst.target_dist = l1->dist_pre;
	dst.target_hash = l1->hash_pre;
	dst.check_hash = 1;
	p = mg_shortest_k(aux->km, aux->g, l1->v^1, 1, &dst, dst.target_dist, MG_MAX_SHORT_K, &n_pathv);
	if (n_pathv == 0 || dst.target_hash != dst.hash) {
		fprintf(stderr, "[W::%s] %c%s[%d] -> %c%s[%d], dist=%d, target_dist=%d; chain skiped.\n", __func__, "><"[(l1->v^1)&1], aux->g->seg[l1->v>>1].name, l1->v^1, "><"[(l0->v^1)&1],
				aux->g->seg[l0->v>>1].name, l0->v^1, dst.dist, dst.target_dist);
		kfree(aux->km, p);
		return -1;
	}
	for (s = n_pathv - 2; s >= 1; --s) { // path found in a backward way, so we need to reverse it
		mg_llchain_t *q;
		if (aux->n_llc == aux->m_llc) KEXPAND(aux->km, aux->llc, aux->m_llc);
		q = &aux->llc[aux->n_llc++];
		q->off = q->cnt = q->score = 0;
		q->v = p[s].v^1; // when reversing a path, we also need to flip the orientation
		q->ed = -1;
	}
	kfree(aux->km, p);
	return 0;
}

static int32_t bridge_gwfa(bridge_aux_t *aux, int32_t kmer_size, int32_t gdp_max_ed, const mg_lchain_t *l0, const mg_lchain_t *l1, int32_t *ed)
{
	uint32_t v0 = l0->v, v1 = l1->v;
	int32_t qs = l0->qe - kmer_size, qe = l1->qs + kmer_size, end0, end1, j;
	void *z;
	gfa_edopt_t opt;
	gfa_edrst_t r;

	*ed = -1;
	end0 = l0->re - kmer_size;
	end1 = l1->rs + kmer_size - 1;

	gfa_edopt_init(&opt);
	opt.traceback = 1, opt.max_chk = 1000, opt.bw_dyn = 1000, opt.max_lag = gdp_max_ed/2;
	opt.i_term = 500000000LL;
	z = gfa_ed_init(aux->km, &opt, aux->g, aux->es, qe - qs, &aux->qseq[qs], v0, end0);
	gfa_ed_step(z, v1, end1, gdp_max_ed, &r);
	gfa_ed_destroy(z);
	//fprintf(stdout, "qs=%d,qe=%d,v0=%c%s:%d:%d,v1=%c%s:%d,s=%d,nv=%d\n", qs, qe, "><"[v0&1], aux->g->seg[v0>>1].name, end0, aux->g->seg[v0>>1].len - end0 - 1, "><"[v1&1], aux->g->seg[v1>>1].name, end1, r.s, r.nv);
	if (r.s < 0) return 0;

	for (j = 1; j < r.nv - 1; ++j) {
		mg_llchain_t *q;
		if (aux->n_llc == aux->m_llc) KEXPAND(aux->km, aux->llc, aux->m_llc);
		q = &aux->llc[aux->n_llc++];
		q->off = q->cnt = q->score = 0;
		q->v = r.v[j];
		q->ed = -1;
	}
	kfree(aux->km, r.v);
	*ed = r.s;
	return 1;
}

static int32_t bridge_lchains(mg_gchains_t *gc, bridge_aux_t *aux, int32_t kmer_size, int32_t gdp_max_ed, const mg_lchain_t *l0, const mg_lchain_t *l1, const mg128_t *a)
{
	if (l1->v != l0->v) { // bridging two segments
		int32_t ed = -1, ret = 0;
		if (aux->n_seg > 1 || !bridge_gwfa(aux, kmer_size, gdp_max_ed, l0, l1, &ed))
			ret = bridge_shortk(aux, l0, l1);
		if (ret < 0) return -1;
		if (aux->n_llc == aux->m_llc) KEXPAND(aux->km, aux->llc, aux->m_llc);
		copy_lchain(&aux->llc[aux->n_llc++], l1, &aux->n_a, gc->a, a, ed);
	} else { // on one segment
		int32_t k;
		mg_llchain_t *t = &aux->llc[aux->n_llc - 1];
		for (k = 0; k < l1->cnt; ++k) { // FIXME: this part is made redundant by resolve_overlap()
			const mg128_t *ak = &a[l1->off + k];
			if ((int32_t)ak->x > l0->re && (int32_t)ak->y > l0->qe)
				break;
		}
		if (k < l1->cnt) { // l1 contained. TODO: check what is happening...
			t->cnt += l1->cnt - k, t->score += l1->score;
			memcpy(&gc->a[aux->n_a], &a[l1->off + k], (l1->cnt - k) * sizeof(mg128_t));
			aux->n_a += l1->cnt - k;
		}
	}
	return 0;
}

static void resolve_overlap(mg_lchain_t *l0, mg_lchain_t *l1, const mg128_t *a)
{
	int32_t j, x, y, shift0, shift1;
	// check the end of l0
	x = (int32_t)a[l1->off].x;
	y = (int32_t)a[l1->off].y;
	for (j = l0->cnt - 1; j >= 0; --j)
		if ((int32_t)a[l0->off + j].y <= y && (l0->v != l1->v || (int32_t)a[l0->off + j].x <= x))
			break;
	shift0 = l0->cnt - 1 - j;
	// check the start of l1
	x = (int32_t)a[l0->off + l0->cnt - 1].x;
	y = (int32_t)a[l0->off + l0->cnt - 1].y;
	for (j = 0; j < l1->cnt; ++j)
		if ((int32_t)a[l1->off + j].y >= y && (l0->v != l1->v || (int32_t)a[l1->off + j].x >= x))
			break;
	shift1 = j;
	assert(shift1 < l1->cnt); // this should never happen, or it is a bug
	// update
	if (shift0 > 0) {
		l0->cnt -= shift0;
		if (l0->cnt) { // l0->cnt may be 0 as the start of l0 may be changed and go into l1
			l0->qe = (int32_t)a[l0->off + l0->cnt - 1].y + 1;
			l0->re = (int32_t)a[l0->off + l0->cnt - 1].x + 1;
		}
	}
	if (shift1 > 0) {
		l1->off += shift1, l1->cnt -= shift1;
		l1->qs = (int32_t)a[l1->off].y + 1 - (int32_t)(a[l1->off].y>>32&0xff);
		l1->rs = (int32_t)a[l1->off].x + 1 - (int32_t)(a[l1->off].y>>32&0xff);
	}
	if (l0->cnt == 0) l0->qs = l0->qe = l1->qs, l0->rs = l0->re = l1->rs; // this line should have no effect
}

mg_gchains_t *mg_gchain_gen(void *km_dst, void *km, const gfa_t *g, const gfa_edseq_t *es, int32_t n_u, const uint64_t *u,
							mg_lchain_t *lc, const mg128_t *a, uint32_t hash, int32_t min_gc_cnt, int32_t min_gc_score,
							int32_t gdp_max_ed, int32_t n_seg, const char *qseq)
{
	mg_gchains_t *gc;
	int32_t i, j, k, st, kmer_size;
	bridge_aux_t aux;

	// preallocate gc->gc and gc->a
	KCALLOC(km_dst, gc, 1);
	for (i = 0, st = 0; i < n_u; ++i) {
		int32_t m = 0, nui = (int32_t)u[i];
		for (j = 0; j < nui; ++j) m += lc[st + j].cnt; // m is the number of anchors in this gchain
		if (m >= min_gc_cnt && u[i]>>32 >= min_gc_score)
			gc->n_gc++, gc->n_a += m;
		st += nui;
	}
	if (gc->n_gc == 0) return gc;
	gc->km = km_dst;
	KCALLOC(km_dst, gc->gc, gc->n_gc);
	KMALLOC(km_dst, gc->a, gc->n_a);

	// core loop
	memset(&aux, 0, sizeof(aux));
	aux.km = km, aux.g = g, aux.es = es, aux.n_seg = n_seg, aux.qseq = qseq;
	kmer_size = a[0].y>>32&0xff;
	for (i = k = 0, st = 0, aux.n_a = 0; i < n_u; ++i) {
		int32_t n_a0 = aux.n_a, n_llc0 = aux.n_llc, m = 0, nui = (int32_t)u[i];
		for (j = 0; j < nui; ++j) m += lc[st + j].cnt;
		if (m >= min_gc_cnt && u[i]>>32 >= min_gc_score) {
			uint32_t h = hash;
			int32_t j0;
			gc->gc[k].score = u[i]>>32;
			gc->gc[k].off = n_llc0;
			for (j = 0; j < nui; ++j) {
				const mg_lchain_t *p = &lc[st + j];
				h += kh_hash_uint32(p->qs) + kh_hash_uint32(p->re) + kh_hash_uint32(p->v);
			}
			gc->gc[k].hash = kh_hash_uint32(h);

			for (j = 1; j < nui; ++j)
				resolve_overlap(&lc[st + j - 1], &lc[st + j], a);

			if (aux.n_llc == aux.m_llc) KEXPAND(aux.km, aux.llc, aux.m_llc);
			copy_lchain(&aux.llc[aux.n_llc++], &lc[st], &aux.n_a, gc->a, a, -1); // copy the first lchain
			for (j0 = 0, j = 1; j < nui; ++j) {
				const mg_lchain_t *l0 = &lc[st + j0], *l1 = &lc[st + j];
				if (l1->cnt > 0) {
					int32_t ret, t;
					ret = bridge_lchains(gc, &aux, kmer_size, gdp_max_ed, l0, l1, a);
					if (ret < 0) {
						for (t = j0; t < j; ++t) {
							ret = bridge_lchains(gc, &aux, kmer_size, gdp_max_ed, &lc[st + t], &lc[st + t + 1], a);
							assert(ret >= 0);
						}
					}
					j0 = j;
				}
			}

			gc->gc[k].cnt = aux.n_llc - n_llc0;
			gc->gc[k].n_anchor = aux.n_a - n_a0;
			++k;
		}
		st += nui;
	}
	assert(aux.n_a <= gc->n_a);

	gc->n_a = aux.n_a;
	gc->n_lc = aux.n_llc;
	KMALLOC(km_dst, gc->lc, aux.n_llc);
	memcpy(gc->lc, aux.llc, aux.n_llc * sizeof(mg_llchain_t));
	kfree(km, aux.llc);

	mg_gchain_extra(g, gc);
	mg_gchain_sort_by_score(km, gc);
	return gc;
}

void mg_gchain_free(mg_gchains_t *gs)
{
	void *km;
	int32_t i;
	if (gs == 0) return;
	km = gs->km;
	for (i = 0; i < gs->n_gc; ++i) {
		if (gs->gc[i].p) kfree(km, gs->gc[i].p);
		if (gs->gc[i].ds.ds) kfree(km, gs->gc[i].ds.ds);
		if (gs->gc[i].ds.off) kfree(km, gs->gc[i].ds.off);
	}
	kfree(km, gs->gc); kfree(km, gs->a); kfree(km, gs->lc);
	kfree(km, gs);
}
