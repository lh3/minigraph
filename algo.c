#include <stdlib.h>
#include <string.h>
#include "kalloc.h"
#define __STDC_LIMIT_MACROS
#include "algo.h"

/**********************************
 * Longest increasing subsequence *
 **********************************/

int32_t ks_lis_64(void *km, int32_t n, const LIS_TYPE *a, int32_t *b)
{
	int32_t i, k, L = 0, *M, *P = b;
	KMALLOC(km, M, n+1);
	for (i = 0; i < n; ++i) {
		int32_t lo = 1, hi = L, newL;
		while (lo <= hi) {
			int32_t mid = (lo + hi + 1) >> 1;
			if (a[M[mid]] < a[i]) lo = mid + 1;
			else hi = mid - 1;
		}
		newL = lo, P[i] = M[newL - 1], M[newL] = i;
		if (newL > L) L = newL;
	}
	k = M[L];
	memcpy(M, P, n * sizeof(int32_t));
	for (i = L - 1; i >= 0; --i) b[i] = k, k = M[k];
	kfree(km, M);
	return L;
}

/************************
 * Max-scoring segments *
 ************************/

#include "kvec.h"

#define MSS_NEG_INF INT32_MIN

typedef struct {
    int32_t st, en;
    MSS_TYPE L, R;
    int32_t pre;
} msseg_aux_t;

typedef kvec_t(msseg_t) msseg_v;
typedef kvec_t(msseg_aux_t) msseg_aux_v;

static void move_segs(void *km, msseg_v *ret, msseg_aux_v *seg, MSS_TYPE min_sc)
{
    int32_t i;
    for (i = 0; i < seg->n; ++i) {
        msseg_aux_t *p = &seg->a[i];
        if (p->R - p->L >= min_sc) {
            msseg_t *q;
            kv_pushp(msseg_t, km, *ret, &q);
            q->st = p->st, q->en = p->en, q->sc = p->R - p->L;
        }
    }
    seg->n = 0;
}

// Reference: Ruzzo and Tompa (1999) A linear time algorithm for finding all maximal scoring subsequencs
msseg_t *mss_find_all(void *km, int32_t n, const MSS_TYPE *S, MSS_TYPE min_sc, MSS_TYPE xdrop, int32_t *n_seg)
{
    int32_t i, j;
    MSS_TYPE L, max;
    msseg_v ret = {0,0,0};
    msseg_aux_v seg = {0,0,0};
    msseg_aux_t t;

	kv_resize(msseg_t, km, ret, 16);
	kv_resize(msseg_aux_t, km, seg, 16);
    for (i = 0, L = 0, max = MSS_NEG_INF; i < n;) {
        if (S[i] > 0) {
            int32_t k;
            MSS_TYPE R = L + S[i];
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
