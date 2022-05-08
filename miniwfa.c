#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "miniwfa.h"
#include "kalloc.h"

/*
 * Default setting
 */
void mwf_opt_init(mwf_opt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->x  = 4; // corresponding SW score: m=1, x=3, o1=4, e1=3/2, o2=15, e2=1/2
	opt->o1 = 4, opt->e1 = 2;
	opt->o2 = 15, opt->e2 = 1;
	opt->kmer = 13, opt->max_occ = 2, opt->min_len = 30;
}

/*
 * Structs and simple functions for traceback
 */
typedef struct {
	int32_t lo, hi;
	uint8_t *x;
} wf_tb1_t;

typedef struct {
	int32_t m, n;
	wf_tb1_t *a;
} wf_tb_t;

static wf_tb1_t *wf_tb_add(void *km, wf_tb_t *tb, int32_t lo, int32_t hi)
{
	wf_tb1_t *p;
	if (tb->n == tb->m) {
		tb->m += (tb->m>>1) + 4;
		tb->a = Krealloc(km, wf_tb1_t, tb->a, tb->m);
	}
	p = &tb->a[tb->n++];
	p->lo = lo, p->hi = hi;
	p->x = Kcalloc(km, uint8_t, hi - lo + 1);
	return p;
}

typedef struct {
	int32_t m, n;
	uint32_t *cigar;
} wf_cigar_t;

static void wf_cigar_push1(void *km, wf_cigar_t *c, int32_t op, int32_t len)
{
	if (c->n && op == (c->cigar[c->n-1]&0xf)) {
		c->cigar[c->n-1] += len<<4;
	} else {
		if (c->n == c->m) {
			c->m = c->m + (c->m>>1) + 4;
			c->cigar = Krealloc(km, uint32_t, c->cigar, c->m);
		}
		c->cigar[c->n++] = len<<4 | op;
	}
}

/*
 * The stripe data structure
 */
#define WF_NEG_INF (-0x40000000)

typedef struct {
	int32_t lo, hi;
	int32_t *mem, *H, *E1, *E2, *F1, *F2;
} wf_slice_t;

typedef struct {
	int32_t s, top, n, max_pen, lo, hi;
	wf_slice_t *a;
} wf_stripe_t;

void wf_stripe_add(void *km, wf_stripe_t *wf, int32_t lo, int32_t hi)
{
	int32_t i, n, m1 = wf->max_pen + 1, m2 = m1 * 2;
	wf_slice_t *f;
	++wf->s;
	++wf->top;
	if (wf->top == wf->n) wf->top = 0;
	f = &wf->a[wf->top];
	f->lo = lo, f->hi = hi;
	n = hi - lo + 1;
	kfree(km, f->mem);
	f->mem = Kmalloc(km, int32_t, 5 * (n + m2));
	f->H = f->mem + m1;
	f->E1 = f->H  + n + m2;
	f->F1 = f->E1 + n + m2;
	f->E2 = f->F1 + n + m2;
	f->F2 = f->E2 + n + m2;
	for (i = -m1; i < 0; ++i)
		f->H[i] = f->E1[i] = f->E2[i] = f->F1[i] = f->F2[i] = WF_NEG_INF;
	for (i = n; i < n + m1; ++i)
		f->H[i] = f->E1[i] = f->E2[i] = f->F1[i] = f->F2[i] = WF_NEG_INF;
	f->H -= lo, f->E1 -= lo, f->E2 -= lo, f->F1 -= lo, f->F2 -= lo; // such that f->H[lo] points to 0
}

static wf_stripe_t *wf_stripe_init(void *km, int32_t max_pen)
{
	int32_t i;
	wf_stripe_t *wf;
	wf = Kcalloc(km, wf_stripe_t, 1);
	wf->max_pen = max_pen;
	wf->n = max_pen + 1;
	wf->a = Kcalloc(km, wf_slice_t, wf->n);
	wf->lo = wf->hi = 0;
	for (i = 0; i < wf->n; ++i) {
		wf_slice_t *f;
		wf_stripe_add(km, wf, 0, 0);
		f = &wf->a[wf->top];
		f->H[0] = f->E1[0] = f->E2[0] = f->F1[0] = f->F2[0] = WF_NEG_INF;
	}
	wf->s = 0;
	wf->a[wf->top].H[0] = -1;
	return wf;
}

static void wf_stripe_destroy(void *km, wf_stripe_t *wf)
{
	int32_t i;
	for (i = 0; i < wf->n; ++i)
		kfree(km, wf->a[i].mem);
	kfree(km, wf->a);
	kfree(km, wf);
}

static inline wf_slice_t *wf_stripe_get(const wf_stripe_t *wf, int32_t x)
{
	int32_t y = wf->top - x;
	if (y < 0) y += wf->n;
	return &wf->a[y];
}

static inline int good_diag(int32_t d, int32_t k, int32_t tl, int32_t ql) // check if (d,k) falls within the DP matrix
{
	return ((k >= -1 && k < tl) && (d + k >= -1 && d + k < ql));
}

static void wf_stripe_shrink(wf_stripe_t *wf, int32_t tl, int32_t ql)
{
	int32_t j, d;
	for (d = wf->lo; d <= wf->hi; ++d) {
		for (j = 0; j < wf->n; ++j) {
			wf_slice_t *p = &wf->a[(wf->top + 1 + j) % wf->n];
			if (d < p->lo || d > p->hi) continue;
			if (good_diag(d, p->H[d], tl, ql)) break;
			if (good_diag(d, p->E1[d], tl, ql) || good_diag(d, p->F1[d], tl, ql)) break;
			if (good_diag(d, p->E2[d], tl, ql) || good_diag(d, p->F2[d], tl, ql)) break;
		}
		if (j < wf->n) break; // stop when we see a "good diagonal" in the stripe
	}
	assert(d <= wf->hi); // should never happen
	wf->lo = d;
	for (d = wf->hi; d >= wf->lo; --d) {
		for (j = 0; j < wf->n; ++j) {
			wf_slice_t *p = &wf->a[(wf->top + 1 + j) % wf->n];
			if (d < p->lo || d > p->hi) continue;
			if (good_diag(d, p->H[d], tl, ql)) break;
			if (good_diag(d, p->E1[d], tl, ql) || good_diag(d, p->F1[d], tl, ql)) break;
			if (good_diag(d, p->E2[d], tl, ql) || good_diag(d, p->F2[d], tl, ql)) break;
		}
		if (j < wf->n) break;
	}
	assert(d >= wf->lo);
	wf->hi = d;
}

typedef struct {
	int32_t s, d;
} wf_chkpt_t;

/*
 * Extend a diagonal along exact matches
 */

// pad strings with distinct characters
static void wf_pad_str(void *km, int32_t tl, const char *ts, int32_t ql, const char *qs, char **pts, char **pqs)
{
	uint8_t t[256];
	int32_t i, c1 = -1, c2 = -1;
	char *s1, *s2;
	*pts = *pqs = 0;
	// collect all used characters
	memset(t, 0, 256);
	for (i = 0; i < tl; ++i)
		if (t[(uint8_t)ts[i]] == 0)
			t[(uint8_t)ts[i]] = 1;
	for (i = 0; i < ql; ++i)
		if (t[(uint8_t)qs[i]] == 0)
			t[(uint8_t)qs[i]] = 1;
	for (i = 0; i < 256; ++i)
		if (t[i] == 0) {
			if (c1 < 0) c1 = i;
			else if (c2 < 0) c2 = i;
		}
	if (c1 < 0 || c2 < 0) return; // The two strings use >=255 characters. Unlikely for bio strings.
	s1 = Kmalloc(km, char, tl + ql + 16); // the two strings are allocated together
	s2 = s1 + tl + 8;
	memcpy(s1, ts, tl);
	for (i = tl; i < tl + 8; ++i) s1[i] = c1; // pad with c1
	memcpy(s2, qs, ql);
	for (i = ql; i < ql + 8; ++i) s2[i] = c2; // pad with c2
	*pts = s1, *pqs = s2;
}

// Extend a diagonal along exact matches.
static inline int32_t wf_extend1_padded(const char *ts, const char *qs, int32_t k, int32_t d)
{
	uint64_t cmp = 0;
	const char *ts_ = ts + 1;
	const char *qs_ = qs + d + 1;
	while (1) {
		uint64_t x = *(uint64_t*)(ts_ + k); // warning: unaligned memory access
		uint64_t y = *(uint64_t*)(qs_ + k);
		cmp = x ^ y;
		if (cmp == 0) k += 8;
		else break;
	}
	k += __builtin_ctzl(cmp) >> 3;
	return k;
}

/*
 * Core wf_next() routines
 */

// Force loop vectorization. Learned from WFA.
#if defined(__clang__)
  #define PRAGMA_LOOP_VECTORIZE _Pragma("clang loop vectorize(enable)")
#elif defined(__GNUC__)
  #define PRAGMA_LOOP_VECTORIZE _Pragma("GCC ivdep")
#else
  #define PRAGMA_LOOP_VECTORIZE _Pragma("ivdep")
#endif

#define wf_max(a, b) ((a) >= (b)? (a) : (b))

static void wf_next_prep(void *km, const mwf_opt_t *opt, wf_stripe_t *wf, int32_t lo, int32_t hi,
						 int32_t **H, int32_t **E1, int32_t **F1, int32_t **E2, int32_t **F2,
						 const int32_t **pHx, const int32_t **pHo1, const int32_t **pHo2,
						 const int32_t **pE1, const int32_t **pF1, const int32_t **pE2, const int32_t **pF2)
{
	const wf_slice_t *fx, *fo1, *fo2, *fe1, *fe2;
	wf_slice_t *ft;
	wf_stripe_add(km, wf, lo, hi);
	ft  = &wf->a[wf->top];
	fx  = wf_stripe_get(wf, opt->x);
	fo1 = wf_stripe_get(wf, opt->o1 + opt->e1);
	fo2 = wf_stripe_get(wf, opt->o2 + opt->e2);
	fe1 = wf_stripe_get(wf, opt->e1);
	fe2 = wf_stripe_get(wf, opt->e2);
	*pHx = fx->H, *pHo1 = fo1->H, *pHo2 = fo2->H, *pE1 = fe1->E1, *pE2 = fe2->E2, *pF1 = fe1->F1, *pF2 = fe2->F2;
	*H = ft->H, *E1 = ft->E1, *E2 = ft->E2, *F1 = ft->F1, *F2 = ft->F2;
}

static void wf_next_score(int32_t lo, int32_t hi, int32_t *H, int32_t *E1, int32_t *F1, int32_t *E2, int32_t *F2,
						  const int32_t *pHx, const int32_t *pHo1, const int32_t *pHo2,
						  const int32_t *pE1, const int32_t *pF1, const int32_t *pE2, const int32_t *pF2)
{
	int32_t d;
	PRAGMA_LOOP_VECTORIZE
	for (d = lo; d <= hi; ++d) {
		int32_t h, f, e;
		E1[d] = wf_max(pHo1[d-1], pE1[d-1]);
		E2[d] = wf_max(pHo2[d-1], pE2[d-1]);
		e = wf_max(E1[d], E2[d]);
		F1[d] = wf_max(pHo1[d+1], pF1[d+1]) + 1;
		F2[d] = wf_max(pHo2[d+1], pF2[d+1]) + 1;
		f = wf_max(F1[d], F2[d]);
		h = wf_max(e, f);
		H[d] = wf_max(pHx[d] + 1, h);
		// if (H[d] >= -1) fprintf(stderr, "s=%d, d=%d, k=%d, (%d,%d)\n", wf->s, d, H[d], E1[d], F1[d]);
	}
}

static void wf_next_tb(int32_t lo, int32_t hi, int32_t *H, int32_t *E1, int32_t *F1, int32_t *E2, int32_t *F2, uint8_t *ax,
					   const int32_t *pHx, const int32_t *pHo1, const int32_t *pHo2,
					   const int32_t *pE1, const int32_t *pF1, const int32_t *pE2, const int32_t *pF2)
{
	int32_t d;
	PRAGMA_LOOP_VECTORIZE
	for (d = lo; d <= hi; ++d) {
		int32_t h, f, e;
		uint8_t x = 0, ze, zf, z;
		x |= pHo1[d-1] >= pE1[d-1]? 0 : 0x08;
		E1[d] = wf_max(pHo1[d-1], pE1[d-1]);
		x |= pHo2[d-1] >= pE2[d-1]? 0 : 0x20;
		E2[d] = wf_max(pHo2[d-1], pE2[d-1]);
		ze = E1[d] >= E2[d]? 1 : 3;
		e = wf_max(E1[d], E2[d]);
		x |= pHo1[d+1] >= pF1[d+1]? 0 : 0x10;
		F1[d] = wf_max(pHo1[d+1], pF1[d+1]) + 1;
		x |= pHo2[d+1] >= pF2[d+1]? 0 : 0x40;
		F2[d] = wf_max(pHo2[d+1], pF2[d+1]) + 1;
		zf = F1[d] >= F2[d]? 2 : 4;
		f = wf_max(F1[d], F2[d]);
		z = e >= f? ze : zf;
		h = wf_max(e, f);
		z = pHx[d] + 1 >= h? 0 : z;
		H[d] = wf_max(pHx[d] + 1, h);
		ax[d] = x | z;
	}
}

/*
 * Core algorithm
 */
static void wf_next_basic(void *km, void *km_tb, const mwf_opt_t *opt, wf_stripe_t *wf, wf_tb_t *tb, int32_t lo, int32_t hi)
{
	int32_t *H, *E1, *E2, *F1, *F2;
	const int32_t *pHx, *pHo1, *pHo2, *pE1, *pE2, *pF1, *pF2;
	wf_next_prep(km, opt, wf, lo, hi, &H, &E1, &F1, &E2, &F2, &pHx, &pHo1, &pHo2, &pE1, &pF1, &pE2, &pF2);
	if (tb) {
		uint8_t *ax;
		ax = wf_tb_add(km_tb, tb, lo, hi)->x - lo;
		wf_next_tb(lo, hi, H, E1, F1, E2, F2, ax, pHx, pHo1, pHo2, pE1, pF1, pE2, pF2);
	} else {
		wf_next_score(lo, hi, H, E1, F1, E2, F2, pHx, pHo1, pHo2, pE1, pF1, pE2, pF2);
	}
	if (H[lo] >= -1 || E1[lo] >= -1 || F1[lo] >= -1 || E2[lo] >= -1 || F2[lo] >= -1) wf->lo = lo;
	if (H[hi] >= -1 || E1[hi] >= -1 || F1[hi] >= -1 || E2[hi] >= -1 || F2[hi] >= -1) wf->hi = hi;
}

static uint32_t *wf_traceback(void *km, const mwf_opt_t *opt, wf_tb_t *tb, int32_t t_end, const char *ts, int32_t q_end, const char *qs, int32_t last, int32_t *n_cigar)
{
	wf_cigar_t cigar = {0,0,0};
	int32_t i = q_end, k = t_end, s = tb->n - 1;
	while (i >= 0 && k >= 0) {
		int32_t k0 = k, j, x, state, ext;
		if (last == 0) { // if the previous state is 0, check exact matches
			while (i >= 0 && k >= 0 && qs[i] == ts[k])
				--i, --k;
			if (k0 - k > 0)
				wf_cigar_push1(km, &cigar, 7, k0 - k);
			if (i < 0 || k < 0) break;
		}
		assert(s >= 0);
		j = i - k - tb->a[s].lo;
		assert(j <= tb->a[s].hi - tb->a[s].lo);
		x = tb->a[s].x[j];
		state = last == 0? x&7 : last;
		ext = state > 0? x>>(state+2)&1 : 0; // whether an extension
		//fprintf(stderr, "s=%d, %d->%d, ext=%d%d%d%d, i=%d, k=%d\n", s, last, state, x>>3&1, x>>4&1, x>>5&1, x>>6&1, i, k);
		if (state == 0) {
			wf_cigar_push1(km, &cigar, 8, 1);
			--i, --k, s -= opt->x;
		} else if (state == 1) {
			wf_cigar_push1(km, &cigar, 1, 1);
			--i, s -= ext? opt->e1 : opt->o1 + opt->e1;
		} else if (state == 3) {
			wf_cigar_push1(km, &cigar, 1, 1);
			--i, s -= ext? opt->e2 : opt->o2 + opt->e2;
		} else if (state == 2) {
			wf_cigar_push1(km, &cigar, 2, 1);
			--k, s -= ext? opt->e1 : opt->o1 + opt->e1;
		} else if (state == 4) {
			wf_cigar_push1(km, &cigar, 2, 1);
			--k, s -= ext? opt->e2 : opt->o2 + opt->e2;
		} else abort();
		last = state > 0 && ext? state : 0;
	}
	if (opt->flag&MWF_F_DEBUG) fprintf(stderr, "s0=%d, s=%d, i=%d, k=%d\n", tb->n-1, s, i, k);
	if (i >= 0) wf_cigar_push1(km, &cigar, 1, i + 1);
	else if (k >= 0) wf_cigar_push1(km, &cigar, 2, k + 1);
	for (i = 0; i < cigar.n>>1; ++i) { // reverse to the input order
		uint32_t t = cigar.cigar[i];
		cigar.cigar[i] = cigar.cigar[cigar.n - i - 1];
		cigar.cigar[cigar.n - i - 1] = t;
	}
	*n_cigar = cigar.n;
	return cigar.cigar;
}

// pts and pqs MUST BE padded with wf_pad_str()
static void mwf_wfa_core(void *km, const mwf_opt_t *opt, int32_t tl, const char *pts, int32_t ql, const char *pqs, int32_t n_seg, wf_chkpt_t *seg, mwf_rst_t *r)
{
	int32_t max_pen, sid, is_tb = !!(opt->flag&MWF_F_CIGAR), last_state = 0, stopped = 0;
	wf_stripe_t *wf;
	wf_tb_t tb = {0,0,0};
	void *km_tb, *km_st;

	memset(r, 0, sizeof(*r));
	km_tb = is_tb && !(opt->flag&MWF_F_NO_KALLOC)? km_init2(km, 0) : 0;
	km_st = !(opt->flag&MWF_F_NO_KALLOC)? km_init2(km, 0) : 0;
	max_pen = opt->x;
	max_pen = max_pen > opt->o1 + opt->e1? max_pen : opt->o1 + opt->e1;
	max_pen = max_pen > opt->o2 + opt->e2? max_pen : opt->o2 + opt->e2;
	wf = wf_stripe_init(km_st, max_pen);
	assert(pts);

	sid = 0;
	while (1) {
		wf_slice_t *p = &wf->a[wf->top];
		int32_t d, lo, hi, *H = p->H;
		for (d = p->lo; d <= p->hi; ++d) {
			int32_t k = 0;
			if (H[d] < -1 || d + H[d] < -1 || H[d] >= tl || d + H[d] >= ql) continue;
			k = wf_extend1_padded(pts, pqs, H[d], d);
			//fprintf(stderr, "[s=%d] [%d,%d]:%d %d->%d,%d,%d,%d,%d\n", wf->s, p->lo, p->hi, d, H[d], k, wf->a[wf->top].E1[d], wf->a[wf->top].F1[d], wf->a[wf->top].E2[d], wf->a[wf->top].F2[d]);
			if (k == tl - 1 && d + k == ql - 1) {
				if (k == H[d] && is_tb)
					last_state = tb.a[tb.n-1].x[d - tb.a[tb.n-1].lo] & 7;
				break;
			}
			H[d] = k;
		}
		if (d <= p->hi) break;
		if (is_tb && seg && sid < n_seg && seg[sid].s == wf->s) {
			assert(seg[sid].d >= wf->lo && seg[sid].d <= wf->hi);
			wf->lo = wf->hi = seg[sid++].d;
		}
		lo = wf->lo > -tl? wf->lo - 1 : -tl;
		hi = wf->hi <  ql? wf->hi + 1 :  ql;
		wf_next_basic(km_st, km_tb, opt, wf, is_tb? &tb : 0, lo, hi);
		if ((wf->s&0xff) == 0) wf_stripe_shrink(wf, tl, ql);
		r->n_iter += hi - lo + 1;
		if ((opt->max_iter > 0 && r->n_iter > opt->max_iter) || (opt->max_s > 0 && wf->s > opt->max_s)) {
			stopped = 1;
			break;
		}
	}
	r->s = stopped? -1 : wf->s;
	if (is_tb && !stopped)
		r->cigar = wf_traceback(km, opt, &tb, tl-1, pts, ql-1, pqs, last_state, &r->n_cigar);
	if (km_st == 0) wf_stripe_destroy(km_st, wf);
	else km_destroy(km_st);
	km_destroy(km_tb);
	if (is_tb && !stopped)
		r->cigar = (uint32_t*)krelocate(km, r->cigar, r->n_cigar * sizeof(*r->cigar));
}

/*
 * Low-memory mode
 */
typedef struct {
	int32_t n, n_intv, max_s;
	int32_t *x;
	uint64_t *intv;
} wf_ss_t; // snapshot

typedef struct {
	int32_t n, m;
	wf_ss_t *a;
} wf_sss_t;

static void wf_snapshot1(void *km, wf_stripe_t *sf, wf_ss_t *ss)
{
	int32_t j, k, t;
	ss->n = 0, ss->max_s = sf->s;
	for (j = 0; j < sf->n; ++j)
		ss->n += 5 * (sf->a[j].hi - sf->a[j].lo + 1);
	ss->x = Kmalloc(km, int32_t, ss->n);
	ss->n_intv = sf->n;
	ss->intv = Kmalloc(km, uint64_t, ss->n_intv);
	for (j = 0, t = 0; j < sf->n; ++j) {
		wf_slice_t *p;
		k = (sf->top + 1 + j) % sf->n;
		p = &sf->a[k];
		ss->intv[j] = (uint64_t)p->lo << 32 | (p->hi - p->lo + 1) * 5;
		for (k = p->lo; k <= p->hi; ++k) {
			ss->x[t] = p->H[k],  p->H[k]  = t++;
			ss->x[t] = p->E1[k], p->E1[k] = t++;
			ss->x[t] = p->F1[k], p->F1[k] = t++;
			ss->x[t] = p->E2[k], p->E2[k] = t++;
			ss->x[t] = p->F2[k], p->F2[k] = t++;
		}
	}
	assert(t == ss->n);
}

static void wf_snapshot(void *km, wf_sss_t *sss, wf_stripe_t *sf)
{
	if (sss->n == sss->m) {
		sss->m += (sss->m>>1) + 8;
		sss->a = Krealloc(km, wf_ss_t, sss->a, sss->m);
	}
	wf_snapshot1(km, sf, &sss->a[sss->n++]);
}

static void wf_snapshot_free(void *km, wf_sss_t *sss)
{
	int32_t j;
	for (j = 0; j < sss->n; ++j) {
		kfree(km, sss->a[j].x);
		kfree(km, sss->a[j].intv);
	}
	kfree(km, sss->a);
}

static void wf_next_seg(void *km, const mwf_opt_t *opt, uint8_t *xbuf, wf_stripe_t *wf, wf_stripe_t *sf, int32_t lo, int32_t hi)
{
	int32_t d, *H, *E1, *E2, *F1, *F2;
	const int32_t *pHx, *pHo1, *pHo2, *pE1, *pE2, *pF1, *pF2;
	uint8_t *ax = xbuf - lo;

	wf_next_prep(km, opt, wf, lo, hi, &H, &E1, &F1, &E2, &F2, &pHx, &pHo1, &pHo2, &pE1, &pF1, &pE2, &pF2);
	wf_next_tb(lo, hi, H, E1, F1, E2, F2, ax, pHx, pHo1, pHo2, pE1, pF1, pE2, pF2);
	wf_next_prep(km, opt, sf, lo, hi, &H, &E1, &F1, &E2, &F2, &pHx, &pHo1, &pHo2, &pE1, &pF1, &pE2, &pF2);
	PRAGMA_LOOP_VECTORIZE
	for (d = lo; d <= hi; ++d) { // FIXME: merge this loop into the loop in wf_next_tb(). I tried but couldn't make clang vectorize.
		uint8_t x = ax[d];
		int32_t a, b, e1, f1, e2, f2, h;
		a = pHo1[d-1], b = pE1[d-1];
		e1 = E1[d] = (x&0x08) == 0? a : b;
		a = pHo1[d+1], b = pF1[d+1];
		f1 = F1[d] = (x&0x10) == 0? a : b;
		a = pHo2[d-1], b = pE2[d-1];
		e2 = E2[d] = (x&0x20) == 0? a : b;
		a = pHo2[d+1], b = pF2[d+1];
		f2 = F2[d] = (x&0x40) == 0? a : b;
		x &= 7;
		h = pHx[d];
		h = x == 1? e1 : h;
		h = x == 2? f1 : h;
		h = x == 3? e2 : h;
		h = x == 4? f2 : h;
		H[d] = h;
	}
	if (H[lo] >= -1 || E1[lo] >= -1 || F1[lo] >= -1 || E2[lo] >= -1 || F2[lo] >= -1) wf->lo = lo;
	if (H[hi] >= -1 || E1[hi] >= -1 || F1[hi] >= -1 || E2[hi] >= -1 || F2[hi] >= -1) wf->hi = hi;
}

static wf_chkpt_t *wf_traceback_seg(void *km, wf_sss_t *sss, int32_t last, int32_t *n_seg)
{
	int32_t j;
	wf_chkpt_t *seg;
	*n_seg = sss->n;
	seg = Kmalloc(km, wf_chkpt_t, sss->n);
	for (j = sss->n - 1; j >= 0; --j) {
		int32_t k, m;
		wf_ss_t *p = &sss->a[j];
		for (k = 0, m = 0; k < p->n_intv; ++k) {
			if (last >= m && last < m + (int32_t)p->intv[k])
				break;
			m += (int32_t)p->intv[k];
		}
		assert(k < p->n_intv);
		seg[j].s = p->max_s - (p->n_intv - k - 1);
		seg[j].d = (int32_t)(p->intv[k]>>32) + (last - m) / 5;
		last = p->x[last];
	}
	assert(last == -1);
	return seg;
}

wf_chkpt_t *mwf_wfa_seg(void *km, const mwf_opt_t *opt, int32_t tl, const char *pts, int32_t ql, const char *pqs, int32_t *n_seg_)
{
	int32_t max_pen, last, n_seg;
	wf_stripe_t *wf, *sf;
	wf_sss_t sss = {0,0,0};
	uint8_t *xbuf;
	wf_chkpt_t *seg;
	void *km_st;

	km_st = !(opt->flag&MWF_F_NO_KALLOC)? km_init2(km, 0) : 0;
	max_pen = opt->x;
	max_pen = max_pen > opt->o1 + opt->e1? max_pen : opt->o1 + opt->e1;
	max_pen = max_pen > opt->o2 + opt->e2? max_pen : opt->o2 + opt->e2;
	xbuf = Kcalloc(km_st, uint8_t, tl + ql + 1);
	wf = wf_stripe_init(km_st, max_pen);
	sf = wf_stripe_init(km_st, max_pen);
	assert(pts);

	while (1) {
		wf_slice_t *p = &wf->a[wf->top];
		int32_t d, lo, hi, *H = p->H;
		for (d = p->lo; d <= p->hi; ++d) {
			int32_t k;
			if (H[d] < -1 || d + H[d] < -1 || H[d] >= tl || d + H[d] >= ql) continue;
			k = wf_extend1_padded(pts, pqs, H[d], d);
			if (k == tl - 1 && d + k == ql - 1) {
				last = sf->a[sf->top].H[d];
				break;
			}
			H[d] = k;
		}
		if (d <= p->hi) break;
		lo = wf->lo > -tl? wf->lo - 1 : -tl;
		hi = wf->hi <  ql? wf->hi + 1 :  ql;
		if ((wf->s + 1) % opt->step == 0)
			wf_snapshot(km_st, &sss, sf);
		wf_next_seg(km_st, opt, xbuf, wf, sf, lo, hi);
		if ((wf->s&0xff) == 0) wf_stripe_shrink(wf, tl, ql);
	}
	seg = wf_traceback_seg(km, &sss, last, &n_seg);
	if (km_st == 0) {
		wf_snapshot_free(km_st, &sss);
		wf_stripe_destroy(km_st, wf);
		wf_stripe_destroy(km_st, sf);
		kfree(km_st, xbuf);
	} else km_destroy(km_st);

	seg = (wf_chkpt_t*)krelocate(km, seg, n_seg * sizeof(*seg));
	*n_seg_ = n_seg;
	return seg;
}

void mwf_wfa_exact(void *km, const mwf_opt_t *opt, int32_t tl, const char *ts, int32_t ql, const char *qs, mwf_rst_t *r)
{
	int32_t n_seg = 0;
	wf_chkpt_t *seg = 0;
	char *pts, *pqs;

	wf_pad_str(km, tl, ts, ql, qs, &pts, &pqs);
	if (opt->step > 0)
		seg = mwf_wfa_seg(km, opt, tl, pts, ql, pqs, &n_seg);
	mwf_wfa_core(km, opt, tl, pts, ql, pqs, n_seg, seg, r);
	kfree(km, seg);
	kfree(km, pts);
}

/*
 * Heuristics
 */
static int32_t mg_lis_64(void *km, int32_t n, const uint64_t *a, int32_t *b)
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

extern void radix_sort_gfa64(uint64_t*, uint64_t*);
extern unsigned char seq_nt4_table[256];

static int32_t mg_fc_kmer(int32_t len, const char *seq, int32_t rid, int32_t k, uint64_t *a)
{
	int32_t i, l, n;
	uint64_t x, mask = (1ULL<<k*2) - 1;
	for (i = l = 0, x = 0, n = 0; i < len; ++i) {
		int32_t c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) {
			x = (x << 2 | c) & mask;
			if (++l >= k) a[n++] = (x<<1|rid) << 32 | i;
		} else l = 0, x = 0;
	}
	return n;
}

static uint64_t *mg_chain(void *km, int32_t l1, const char *s1, int32_t l2, const char *s2, int32_t k, int32_t max_occ, int32_t *n_lis_)
{
	int32_t i, n_a, n_b, m_b, i0, n_lis, *lis;
	uint64_t *a, *b;

	*n_lis_ = 0;
	if (l1 < k || l2 < k) return 0;
	assert(k >= 2 && k <= 15);

	// collect k-mers
	KMALLOC(km, a, l1 + l2);
	n_a = mg_fc_kmer(l1, s1, 0, k, a);
	n_a += mg_fc_kmer(l2, s2, 1, k, &a[n_a]);
	radix_sort_gfa64(a, a + n_a);

	// collect k-mer matches
	n_b = m_b = 0, b = 0;
	for (i0 = 0, i = 1; i <= n_a; ++i) {
		if (i == n_a || a[i0]>>33 != a[i]>>33) {
			if (i - i0 >= 2) {
				int32_t j, s, t;
				for (j = i0; j < i && (a[j]>>32&1) == 0; ++j) {}
				if (j > i0 && j < i && j - i0 <= max_occ && i - j <= max_occ) {
					for (s = i0; s < j; ++s)
						for (t = j; t < i; ++t) {
							if (n_b == m_b) KEXPAND(km, b, m_b);
							b[n_b++] = a[s]<<32 | (uint32_t)a[t];
						}
				}
			}
			i0 = i;
		}
	}
	kfree(km, a);

	// find co-linear chain with LIS
	radix_sort_gfa64(b, b + n_b);
	for (i = 0; i < n_b; ++i)
		b[i] = b[i]>>32 | b[i]<<32;
	KMALLOC(km, lis, n_b);
	n_lis = mg_lis_64(km, n_b, b, lis);
	a = Kmalloc(km, uint64_t, n_lis);
	for (i = 0; i < n_lis; ++i) a[i] = b[lis[i]];
	kfree(km, lis);
	kfree(km, b);
	b = Kmalloc(km, uint64_t, n_lis);
	memcpy(b, a, sizeof(uint64_t) * n_lis);
	kfree(km, a);
	*n_lis_ = n_lis;
	for (i = 0; i < n_lis; ++i) // switch back, such that seq1 on the high bits
		b[i] = b[i]>>32 | b[i]<<32;
	return b;
}

static double mwf_ksim(void *km, int32_t l1, const char *s1, int32_t l2, const char *s2, int32_t k)
{
	int32_t i, i0, j, n_a, n1 = 0, n2 = 0, t1 = 0, t2 = 0;
	double p1, p2;
	uint64_t *a;
	if (l1 < k || l2 < k) return 0;
	assert(k >= 2 && k <= 15);
	KMALLOC(km, a, l1 + l2);
	n_a = mg_fc_kmer(l1, s1, 0, k, a);
	n_a += mg_fc_kmer(l2, s2, 1, k, &a[n_a]);
	radix_sort_gfa64(a, a + n_a);
	for (i0 = 0, i = 1; i <= n_a; ++i) {
		if (i == n_a || a[i0]>>33 != a[i]>>33) {
			int32_t m1, m2, min;
			for (j = i0; j < i && (a[j]>>32&1) == 0; ++j) {}
			m1 = j - i0, m2 = i - j;
			min = m1 < m2? m1 : m2;
			n1 += m1, n2 += m2;
			if (m1 > 0 && m2 > 0)
				t1 += min, t2 += min;
			i0 = i;
		}
	}
	kfree(km, a);
	p1 = (double)t1 / n1, p2 = (double)t2 / n2;
	return p1 > p2? p1 : p2;
}

#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

static void wf_cigar_push(void *km, wf_cigar_t *c, int32_t n_cigar, const uint32_t *cigar)
{
	if (n_cigar == 0) return;
	wf_cigar_push1(km, c, cigar[0]&0xf, cigar[0]>>4);
	if (c->n + n_cigar - 1 > c->m) {
		c->m = c->n + n_cigar - 1;
		kroundup32(c->m);
		c->cigar = Krealloc(km, uint32_t, c->cigar, c->m);
	}
	memcpy(&c->cigar[c->n], &cigar[1], sizeof(*cigar) * (n_cigar - 1));
	c->n += n_cigar - 1;
}

static int32_t wf_anchor_filter(int32_t n, uint64_t *a, int32_t tl, int32_t ql, int32_t k, int32_t min_l)
{
	int32_t i, st, x0, y0, x1, y1, j, l, m;
	for (i = 0, x0 = y0 = x1 = y1 = 0, st = -1, l = 0; i <= n; ++i) {
		int32_t x, y;
		if (i == n) x = tl, y = ql;
		else x = (int32_t)(a[i]>>32) + 1, y = (int32_t)a[i] + 1;
		if (x - x0 != y - y0) {
			//fprintf(stderr, "X\t%d\t(%d,%d) -> (%d,%d)\n", l, x0, y0, x, y);
			if (l < min_l)
				for (j = st > 0? st : 0; j < i; ++j)
					a[j] = 0;
			x0 = x, y0 = y, st = i, l = k;
		} else l += x - x1;
		x1 = x, y1 = y;
	}
	for (i = 0, m = 0; i < n; ++i)
		if (a[i] != 0) a[m++] = a[i];
	return m;
}

void mwf_wfa_chain(void *km, const mwf_opt_t *opt, int32_t tl, const char *ts, int32_t ql, const char *qs, mwf_rst_t *r)
{
	int32_t n_a, i, x0, y0;
	uint64_t *a;
	void *km_wfa;
	wf_cigar_t c = {0,0,0};

	km_wfa = !(opt->flag&MWF_F_NO_KALLOC)? km_init2(km, 0) : 0;
	a = mg_chain(km_wfa, tl, ts, ql, qs, opt->kmer, opt->max_occ, &n_a);
	n_a = wf_anchor_filter(n_a, a, tl, ql, opt->kmer, opt->min_len);
	r->s = 0;
	for (i = 0, x0 = y0 = 0; i <= n_a; ++i) {
		int32_t x1, y1;
		if (i == n_a) x1 = tl, y1 = ql;
		else x1 = (int32_t)(a[i]>>32) + 1, y1 = (int32_t)a[i] + 1;
		if (i < n_a && x1 - x0 == y1 - y0 && x1 - x0 <= opt->kmer) {
			if (opt->flag&MWF_F_CIGAR)
				wf_cigar_push1(km, &c, 7, x1 - x0);
		} else if (x0 < x1 && y0 < y1) {
			if (x1 - x0 >= 10000 && y1 - y0 >= 10000 && mwf_ksim(km, x1 - x0, &ts[x0], y1 - y0, &qs[y0], opt->kmer) < 0.02) {
				if (opt->flag&MWF_F_CIGAR) {
					wf_cigar_push1(km, &c, 2, x1 - x0);
					wf_cigar_push1(km, &c, 1, y1 - y0);
				}
				r->s += opt->o2 * 2 + opt->e2 * ((x1 - x0) + (y1 - y0));
			} else {
				mwf_rst_t q;
				mwf_wfa_exact(km_wfa, opt, x1 - x0, &ts[x0], y1 - y0, &qs[y0], &q);
				if (opt->flag&MWF_F_CIGAR)
					wf_cigar_push(km, &c, q.n_cigar, q.cigar);
				r->s += q.s;
				kfree(km_wfa, q.cigar);
			}
		} else if (x0 < x1) {
			wf_cigar_push1(km, &c, 2, x1 - x0);
			r->s += opt->o2 + (x1 - x0) * opt->e2 < opt->o1 + (x1 - x0) * opt->e1? opt->o2 + (x1 - x0) * opt->e2 : opt->o1 + (x1 - x0) * opt->e1;
		} else if (y0 < y1) {
			wf_cigar_push1(km, &c, 1, y1 - y0);
			r->s += opt->o2 + (y1 - y0) * opt->e2 < opt->o1 + (y1 - y0) * opt->e1? opt->o2 + (y1 - y0) * opt->e2 : opt->o1 + (y1 - y0) * opt->e1;
		}
		x0 = x1, y0 = y1;
	}
	if (km_wfa == 0) kfree(km_wfa, a);
	km_destroy(km_wfa);
	r->n_cigar = c.n, r->cigar = c.cigar;
	r->cigar = (uint32_t*)krelocate(km, r->cigar, r->n_cigar * sizeof(*r->cigar));
}

void mwf_wfa_auto(void *km, const mwf_opt_t *opt0, int32_t tl, const char *ts, int32_t ql, const char *qs, mwf_rst_t *r)
{
	mwf_opt_t opt = *opt0;
	opt.step = 0, opt.max_iter = 100000000;
	mwf_wfa_exact(km, &opt, tl, ts, ql, qs, r);
	if (r->s < 0) {
		if (opt.flag & MWF_F_CIGAR) opt.step = 5000;
		opt.max_iter = -1;
		mwf_wfa_chain(km, &opt, tl, ts, ql, qs, r);
	}
}
