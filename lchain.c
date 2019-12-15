#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "mgpriv.h"
#include "kalloc.h"
#include "kavl.h"

uint64_t *mg_chain_backtrack(void *km, int64_t n, const int32_t *f, const int32_t *p, int32_t *v, int32_t *t, int32_t min_cnt, int32_t min_sc, int32_t extra_u, int32_t *n_u_, int32_t *n_v_)
{
	mg128_t *z;
	uint64_t *u;
	int64_t i, k, n_z, n_v;
	int32_t n_u;

	*n_u_ = *n_v_ = 0;
	for (i = 0, n_z = 0; i < n; ++i) // precompute n_z
		if (f[i] >= min_sc) ++n_z;
	if (n_z == 0) return 0;
	KMALLOC(km, z, n_z);
	for (i = 0, k = 0; i < n; ++i) // populate z[]
		if (f[i] >= min_sc) z[k].x = f[i], z[k++].y = i;
	radix_sort_128x(z, z + n_z);

	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // precompute n_u
		int64_t n_v0 = n_v;
		int32_t sc;
		for (i = z[k].y; i >= 0 && t[i] == 0; i = p[i])
			++n_v, t[i] = 1;
		sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
		if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
			++n_u;
		else n_v = n_v0;
	}
	KMALLOC(km, u, n_u + extra_u);
	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // populate u[]
		int64_t n_v0 = n_v;
		int32_t sc;
		for (i = z[k].y; i >= 0 && t[i] == 0; i = p[i])
			v[n_v++] = i, t[i] = 1;
		sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
		if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
			u[n_u++] = (uint64_t)sc << 32 | (n_v - n_v0);
		else n_v = n_v0;
	}
	kfree(km, z);
	assert(n_v < INT32_MAX);
	*n_u_ = n_u, *n_v_ = n_v;
	return u;
}

static inline int32_t comput_sc(const mg128_t *ai, const mg128_t *aj, int32_t max_dist_x, int32_t max_dist_y, int32_t bw, float chn_pen_gap, float chn_pen_skip, int is_cdna, int n_segs)
{
	int32_t dq = (int32_t)ai->y - (int32_t)aj->y, dr, dd, dg, same_seg, q_span, sc;
	float lin_pen, log_pen;
	if (dq <= 0 || dq > max_dist_x) return INT32_MIN;
	same_seg = ((ai->y & MG_SEED_SEG_MASK) == (aj->y & MG_SEED_SEG_MASK));
	dr = (int32_t)(ai->x - aj->x);
	if (same_seg && (dq > max_dist_y || dr == 0)) return INT32_MIN; // don't skip if an anchor is used by multiple segments
	if (n_segs > 1 && !is_cdna && same_seg && dr > max_dist_y) return INT32_MIN;
	dd = dr > dq? dr - dq : dq - dr;
	if (same_seg && dd > bw) return INT32_MIN;
	dg = dr < dq? dr : dq;
	q_span = aj->y>>32&0xff;
	sc = q_span < dg? q_span : dg;
	if (aj->y>>MG_SEED_WT_SHIFT < 255) {
		int tmp = (int)(0.00392156862745098 * (aj->y>>MG_SEED_WT_SHIFT) * sc); // 0.00392... = 1/255
		sc = tmp > 1? tmp : 1;
	}
	lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
	log_pen = dd >= 2? mg_log2(dd) : 0.0f; // mg_log2() only works for dd>=2
	if (is_cdna || !same_seg) {
		if (!same_seg && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
		else if (dr > dq || !same_seg) sc -= (int)(lin_pen < log_pen? lin_pen : log_pen); // deletion or jump between paired ends
		else sc -= (int)(lin_pen + log_pen);
	} else sc -= (int)(lin_pen + log_pen);
	return sc;
}

/* Input:
 *   a[].x: tid<<33 | rev<<32 | tpos
 *   a[].y: flags<<40 | q_span<<32 | q_pos
 * Output:
 *   n_u: #chains
 *   u[]: score<<32 | #anchors (sum of lower 32 bits of u[] is the returned length of a[])
 * input a[] is deallocated on return
 */
mg128_t *mg_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					  int is_cdna, int n_segs, int64_t n, mg128_t *a, int *n_u_, uint64_t **_u, void *km)
{ // TODO: make sure this works when n has more than 32 bits
	int32_t k, *f, *p, *t, *v, n_u, n_v;
	int64_t i, j, max_ii, st = 0;
	uint64_t *u, *u2;
	mg128_t *b, *w;

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) return 0;
	f = (int32_t*)kmalloc(km, n * 4);
	p = (int32_t*)kmalloc(km, n * 4);
	t = (int32_t*)kmalloc(km, n * 4);
	v = (int32_t*)kmalloc(km, n * 4);
	memset(t, 0, n * 4);

	// fill the score and backtrack arrays
	int64_t n_iter = 0; int32_t mmax_f = 0;
	for (i = 0, max_ii = -1; i < n; ++i) {
		int64_t max_j = -1, end_j;
		int32_t max_f = a[i].y>>32&0xff, n_skip = 0;
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist_x)) ++st;
		if (i - st > max_iter) st = i - max_iter;
		for (j = i - 1; j >= st; --j) {
			int32_t sc;
			sc = comput_sc(&a[i], &a[j], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_segs);
			++n_iter;
			if (sc == INT32_MIN) continue;
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == i) {
				if (++n_skip > max_skip)
					break;
			}
			if (p[j] >= 0) t[p[j]] = i;
		}
		end_j = j;
		if (max_ii < 0 || a[i].x - a[max_ii].x > (int64_t)max_dist_x) {
			int32_t max = INT32_MIN;
			max_ii = -1;
			for (j = i - 1; j >= st; --j)
				if (max < f[j]) max = f[j], max_ii = j;
		}
		if (max_ii >= 0 && max_ii < end_j) {
			int32_t tmp;
			tmp = comput_sc(&a[i], &a[max_ii], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_segs);
			if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		if (max_ii < 0 || (a[i].x - a[max_ii].x <= (int64_t)max_dist_x && f[max_ii] < f[i]))
			max_ii = i;
		if (mmax_f < max_f) mmax_f = max_f;
	}
	fprintf(stderr, "Z\tn_iter=%ld\tmmax_f=%d\n", (long)n_iter, mmax_f);

	u = mg_chain_backtrack(km, n, f, p, v, t, min_cnt, min_sc, 0, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, f); kfree(km, p); kfree(km, t);
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}

	// write the result to b[]
	KMALLOC(km, b, n_v);
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			b[k++] = a[v[k0 + (ni - j - 1)]];
	}
	kfree(km, v);

	// sort u[] and a[] by the target position, such that adjacent chains may be joined
	KMALLOC(km, w, n_u);
	for (i = k = 0; i < n_u; ++i) {
		w[i].x = b[k].x, w[i].y = (uint64_t)k<<32|i;
		k += (int32_t)u[i];
	}
	radix_sort_128x(w, w + n_u);
	KMALLOC(km, u2, n_u);
	for (i = k = 0; i < n_u; ++i) {
		int32_t j = (int32_t)w[i].y, n = (int32_t)u[j];
		u2[i] = u[j];
		memcpy(&a[k], &b[w[i].y>>32], n * sizeof(mg128_t));
		k += n;
	}
	memcpy(u, u2, n_u * 8);
	memcpy(b, a, k * sizeof(mg128_t)); // write _a_ to _b_ and deallocate _a_ because _a_ is oversized, sometimes a lot
	kfree(km, a); kfree(km, w); kfree(km, u2);
	return b;
}

typedef struct lc_elem_s {
	int32_t y;
	int64_t i;
	KAVL_HEAD(struct lc_elem_s) head;
} lc_elem_t;

#define lc_elem_cmp(a, b) ((a)->y < (b)->y? -1 : (a)->y > (b)->y? 1 : ((a)->i > (b)->i) - ((a)->i < (b)->i))
KAVL_INIT(lc_elem, lc_elem_t, head, lc_elem_cmp)

mg128_t *mg_lchain_alt(int max_dist, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip, int64_t n, mg128_t *a, int *n_u_, uint64_t **_u, void *km)
{
	int32_t k, *f, *p, *t, *v, n_u, n_v, n_del = 0, m_del = 0;
	int64_t i, j, st = 0;
	uint64_t *u, *u2;
	mg128_t *b, *w;
	lc_elem_t *root = 0, **del = 0;
	void *mem;

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) return 0;
	f = (int32_t*)kmalloc(km, n * 4);
	p = (int32_t*)kmalloc(km, n * 4);
	t = (int32_t*)kmalloc(km, n * 4);
	v = (int32_t*)kmalloc(km, n * 4);
	mem = km_init2(km, 0x10000);
	memset(t, 0, n * 4);

	// fill the score and backtrack arrays
	int64_t n_iter = 0; int32_t mmax_f = 0, max_size = 0, max_del = 0;
	for (i = 0; i < n; ++i) {
		int64_t max_j = -1;
		int32_t q_span = a[i].y>>32&0xff, max_f = q_span;
		lc_elem_t t, *q, *lower, *upper;
		const lc_elem_t *r;
		kavl_itr_t(lc_elem) itr;
		if (max_size < kavl_size(head, root)) max_size = kavl_size(head, root);
		// get rid of active chains out of range
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist)) {
			t.y = (int32_t)a[st].y, t.i = st;
			q = kavl_find(lc_elem, root, &t, 0);
			if (q) {
				q = kavl_erase(lc_elem, &root, q, 0);
				kfree(mem, q);
			}
			++st;
		}
		// traverse the neighbors
		t.i = 0, t.y = (int32_t)a[i].y;
		if (t.y < 0) t.y = 0;
		kavl_interval(lc_elem, root, &t, &lower, &upper);
		if (lower == 0) goto skip_tree;
		kavl_itr_find(lc_elem, root, lower, &itr);
		fprintf(stderr, "Y2\t%ld\tcut=%d\tsize=%d\tn_iter=%ld\n", (long)i, t.y, kavl_size(head, root), (long)n_iter);
		while ((r = kavl_at(&itr)) != 0) {
			int64_t j = r->i;
			int32_t sc, dq, dr, dd, dg;
			float lin_pen, log_pen;
			dq = (int32_t)a[i].y - (int32_t)a[j].y;
			if (dq > max_dist) break;
			++n_iter;
			fprintf(stderr, "X1\t(%d,%d) -> (%d,%d)\n", (int32_t)a[j].x, (int32_t)a[j].y, (int32_t)a[i].x, (int32_t)a[i].y);
			dr = (int32_t)(a[i].x - a[j].x);
			dd = dr > dq? dr - dq : dq - dr;
			dg = dr < dq? dr : dq;
			sc = q_span < dg? q_span : dg;
			if (a[j].y>>MG_SEED_WT_SHIFT < 255) {
				int tmp = (int)(0.00392156862745098 * (a[j].y>>MG_SEED_WT_SHIFT) * sc); // 0.00392... = 1/255
				sc = tmp > 1? tmp : 1;
			}
			lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
			log_pen = dd >= 2? mg_log2(dd) : 0.0f; // mg_log2() only works for dd>=2
			sc -= (int32_t)(lin_pen + log_pen);
			sc += f[j];
			if (sc > max_f) max_f = sc, max_j = j;
			if (!kavl_itr_prev(lc_elem, &itr)) break;
			//if (!kavl_itr_next(lc_elem, &itr)) break;
		}
		// update the tree
		if (upper == 0) goto skip_tree;
		kavl_itr_find(lc_elem, root, upper, &itr);
		n_del = 0;
		/*
		while ((r = kavl_at(&itr)) != 0) {
			int64_t j = r->i;
			int32_t dq, dr, dd;
			float thres;
			dq = (int32_t)a[j].y - (int32_t)a[i].y;
			if (dq > max_dist) break;
			++n_iter;
			dr = (int32_t)(a[i].x - a[j].x);
			dd = dq > dr? dq - dr : dr - dq;
			thres = max_f - chn_pen_gap * dd;
			if (thres < 0.0f) thres = 0.0f;
			if (f[j] - chn_pen_skip * dq < thres) {
				if (n_del == m_del) KEXPAND(km, del, m_del);
				del[n_del++] = (lc_elem_t*)r;
			}
			//fprintf(stderr, "X2\t(%d,%d) -> (%d,%d)\tmax_f=%d\tf[j]=%d\n", (int32_t)a[j].x, (int32_t)a[j].y, (int32_t)a[i].x, (int32_t)a[i].y, max_f, f[j]);
			if (!kavl_itr_next(lc_elem, &itr)) break;
		}
		*/
		for (j = 0; j < n_del; ++j) {
			q = kavl_erase(lc_elem, &root, del[j], 0);
			kfree(mem, q);
		}
		if (max_del < n_del) max_del = n_del;
skip_tree:
		KMALLOC(mem, q, 1);
		q->y = (int32_t)a[i].y, q->i = i;
		q = kavl_insert(lc_elem, &root, q, 0);
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		if (max_f > mmax_f) mmax_f = max_f;
		//fprintf(stderr, "Y3\ti=%ld\t%d\n", (long)i, kavl_size(head, root));
	}
	fprintf(stderr, "Z\tn=%ld\tn_iter=%ld\tmmax_f=%d\tmax_size=%d\tmax_del=%d\n", (long)n, (long)n_iter, mmax_f, max_size, max_del);
	kfree(km, del);
	km_destroy(mem);

	u = mg_chain_backtrack(km, n, f, p, v, t, min_cnt, min_sc, 0, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, f); kfree(km, p); kfree(km, t);
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}

	// write the result to b[]
	KMALLOC(km, b, n_v);
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			b[k++] = a[v[k0 + (ni - j - 1)]];
	}
	kfree(km, v);

	// sort u[] and a[] by the target position, such that adjacent chains may be joined
	KMALLOC(km, w, n_u);
	for (i = k = 0; i < n_u; ++i) {
		w[i].x = b[k].x, w[i].y = (uint64_t)k<<32|i;
		k += (int32_t)u[i];
	}
	radix_sort_128x(w, w + n_u);
	KMALLOC(km, u2, n_u);
	for (i = k = 0; i < n_u; ++i) {
		int32_t j = (int32_t)w[i].y, n = (int32_t)u[j];
		u2[i] = u[j];
		memcpy(&a[k], &b[w[i].y>>32], n * sizeof(mg128_t));
		k += n;
	}
	memcpy(u, u2, n_u * 8);
	memcpy(b, a, k * sizeof(mg128_t)); // write _a_ to _b_ and deallocate _a_ because _a_ is oversized, sometimes a lot
	kfree(km, a); kfree(km, w); kfree(km, u2);
	return b;
}

mg_lchain_t *mg_lchain_gen(void *km, uint32_t hash, int qlen, int n_u, uint64_t *u, mg128_t *a)
{
	mg128_t *z;
	mg_lchain_t *r;
	int i, k;

	if (n_u == 0) return 0;
	KCALLOC(km, r, n_u);

	// sort by query position
	KMALLOC(km, z, n_u);
	for (i = k = 0; i < n_u; ++i) {
		int32_t qs = (int32_t)a[k].y + 1 - (a[k].y>>32 & 0xff);
		z[i].x = (uint64_t)qs << 32 | u[i] >> 32;
		z[i].y = (uint64_t)k << 32 | (int32_t)u[i];
		k += (int32_t)u[i];
	}
	radix_sort_128x(z, z + n_u);

	// populate r[]
	for (i = 0; i < n_u; ++i) {
		mg_lchain_t *ri = &r[i];
		int32_t k = z[i].y >> 32, q_span = a[k].y >> 32 & 0xff;
		ri->off = k;
		ri->cnt = (int32_t)z[i].y;
		ri->score = (uint32_t)z[i].x;
		ri->v = a[k].x >> 32;
		ri->rs = (int32_t)a[k].x + 1 > q_span? (int32_t)a[k].x + 1 - q_span : 0; // for HPC k-mer
		ri->qs = z[i].x >> 32;
		ri->re = (int32_t)a[k + ri->cnt - 1].x + 1;
		ri->qe = (int32_t)a[k + ri->cnt - 1].y + 1;
	}
	kfree(km, z);
	return r;
}

static int32_t get_mini_idx(const mg128_t *a, int32_t n, const int32_t *mini_pos)
{
	int32_t x, L = 0, R = n - 1;
	x = (int32_t)a->y;
	while (L <= R) { // binary search
		int32_t m = ((uint64_t)L + R) >> 1;
		int32_t y = mini_pos[m];
		if (y < x) L = m + 1;
		else if (y > x) R = m - 1;
		else return m;
	}
	return -1;
}

/* Before:
 *   a[].x: tid<<33 | rev<<32 | tpos
 *   a[].y: flags<<40 | q_span<<32 | q_pos
 * After:
 *   a[].x: mini_pos<<32 | tpos
 *   a[].y: same
 */
void mg_update_anchors(int32_t n_a, mg128_t *a, int32_t n, const int32_t *mini_pos)
{
	int32_t st, j, k;
	if (n_a <= 0) return;
	st = get_mini_idx(&a[0], n, mini_pos);
	assert(st >= 0);
	for (k = 0, j = st; j < n && k < n_a; ++j)
		if ((int32_t)a[k].y == mini_pos[j])
			a[k].x = (uint64_t)j << 32 | (a[k].x & 0xffffffffU), ++k;
	assert(k == n_a);
}
