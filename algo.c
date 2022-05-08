#include <stdlib.h>
#include <string.h>
#include "kalloc.h"
#define __STDC_LIMIT_MACROS
#include "algo.h"
#include "miniwfa.h"

/************************
 * Max-scoring segments *
 ************************/

#include "kvec-km.h"

#define MSS_NEG_INF INT32_MIN

typedef struct {
    int32_t st, en;
    MG_MSS_TYPE L, R;
    int32_t pre;
} msseg_aux_t;

typedef kvec_t(mg_msseg_t) msseg_v;
typedef kvec_t(msseg_aux_t) msseg_aux_v;

static void move_segs(void *km, msseg_v *ret, msseg_aux_v *seg, MG_MSS_TYPE min_sc)
{
    int32_t i;
    for (i = 0; i < seg->n; ++i) {
        msseg_aux_t *p = &seg->a[i];
        if (p->R - p->L >= min_sc) {
            mg_msseg_t *q;
            kv_pushp(mg_msseg_t, km, *ret, &q);
            q->st = p->st, q->en = p->en, q->sc = p->R - p->L;
        }
    }
    seg->n = 0;
}

// Reference: Ruzzo and Tompa (1999) A linear time algorithm for finding all maximal scoring subsequencs
mg_msseg_t *mg_mss_all(void *km, int32_t n, const MG_MSS_TYPE *S, MG_MSS_TYPE min_sc, MG_MSS_TYPE xdrop, int32_t *n_seg)
{
    int32_t i, j;
    MG_MSS_TYPE L, max;
    msseg_v ret = {0,0,0};
    msseg_aux_v seg = {0,0,0};
    msseg_aux_t t;

	kv_resize(mg_msseg_t, km, ret, 16);
	kv_resize(msseg_aux_t, km, seg, 16);
    for (i = 0, L = 0, max = MSS_NEG_INF; i < n;) {
        if (S[i] > 0) {
            int32_t k;
            MG_MSS_TYPE R = L + S[i];
            for (k = i + 1; k < n && S[k] > 0; ++k)
                R += S[k];
			if (R > max) max = R;
            t.st = i, t.en = k, t.L = L, t.R = R;
            while (1) {
                msseg_aux_t *p;
                for (j = seg.n - 1; j >= 0;) {
                    p = &seg.a[j];
                    if (p->L < t.L) break;
                    j = p->pre >= 0? p->pre : j - 1;
                }
                if (j >= 0 && seg.a[j].R < t.R) {
                    p = &seg.a[j];
                    t.st = p->st, t.L = p->L, t.pre = p->pre;
                    seg.n = j;
                } else {
                    if (j < 0) {
						move_segs(km, &ret, &seg, min_sc);
						max = R;
					}
                    t.pre = j;
                    kv_push(msseg_aux_t, km, seg, t);
                    break;
                }
            }
            L = R, i = k;
        } else {
			if (xdrop > 0 && L + S[i] + xdrop < max) { // reset
				move_segs(km, &ret, &seg, min_sc);
				L = 0, max = MSS_NEG_INF;
			}
			L += S[i++];
		}
    }
    move_segs(km, &ret, &seg, min_sc);
    kfree(km, seg.a);
	KREALLOC(km, ret.a, ret.n);
    *n_seg = ret.n;
    return ret.a;
}

/**************************
 * Interval overlap query *
 **************************/

#include <assert.h>
#include "ksort.h"

#define sort_key_intv(a) ((a).st)
KRADIX_SORT_INIT(mg_intv, mg_intv_t, sort_key_intv, 4)

int32_t mg_intv_index(int32_t n, mg_intv_t *a)
{
	int32_t i, last_i, last, k;
	if (n <= 0) return -1;
	radix_sort_mg_intv(a, a + n);
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

int32_t mg_intv_overlap(void *km, int32_t n_a, const mg_intv_t *a, int32_t st, int32_t en, int32_t **b_, int32_t *m_b_)
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

/********************
 * Global alignment *
 ********************/

int32_t mg_wfa_cmp(void *km, int32_t l1, const char *s1, int32_t l2, const char *s2, int32_t max_pen, int32_t *mlen, int32_t *blen)
{
	mwf_opt_t opt;
	mwf_rst_t r;
	int32_t i;
	mwf_opt_init(&opt);
	opt.max_s = max_pen;
	opt.flag |= MWF_F_CIGAR;
	mwf_wfa_exact(km, &opt, l1, s1, l2, s2, &r);
	*mlen = *blen = 0;
	for (i = 0; i < r.n_cigar; ++i) {
		int32_t op = r.cigar[i]&0xf, len = r.cigar[i]>>4;
		*blen += len;
		if (op == 7) *mlen += len;
	}
	kfree(km, r.cigar);
	return r.s < 0? -(l1 + l2) : (l1 + l2) / 2 - r.s;
}
