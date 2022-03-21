#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "kalloc.h"

typedef struct {
	int32_t d;
	int32_t k:30, p:2;
} wf_diag_t;

typedef struct {
	int32_t n, d0;
	uint64_t *a;
} wf_tb1_t;

typedef struct {
	int32_t m, n;
	wf_tb1_t *a;
} wf_tb_t;

static void wf_tb_add(void *km, wf_tb_t *tb, int32_t m)
{
	wf_tb1_t *p;
	if (tb->n == tb->m) {
		tb->m = tb->m + (tb->m>>1) + 4;
		KREALLOC(km, tb->a, tb->m);
	}
	p = &tb->a[tb->n++];
	p->n = m;
	KCALLOC(km, p->a, (m + 31) >> 5);
}

static int32_t wf_step(void *km, wf_tb_t *tb, int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t n, wf_diag_t *a, int32_t *t_end, int32_t *q_end)
{
	int32_t j;
	wf_tb1_t *q;
	wf_diag_t *b = a + n + 2; // temporary array

	// wfa_extend
	*t_end = -1;
	for (j = 0; j < n; ++j) {
		wf_diag_t *p = &a[j];
		int32_t k = p->k, max_k;
		const char *ts_, *qs_;
		uint64_t cmp = 0;
		if (k >= tl || k + p->d >= ql) continue;
		max_k = (ql - p->d < tl? ql - p->d : tl) - 1;
		ts_ = ts + 1;
		qs_ = qs + p->d + 1;
		while (k + 7 < max_k) {
			uint64_t x = *(uint64_t*)(ts_ + k); // warning: unaligned memory access
			uint64_t y = *(uint64_t*)(qs_ + k);
			cmp = x ^ y;
			if (cmp == 0) k += 8;
			else break;
		}
		if (cmp)
			k += __builtin_ctzl(cmp) >> 3; // on x86, this is done via the BSR instruction: https://www.felixcloutier.com/x86/bsr
		else if (k + 7 >= max_k)
			while (k < max_k && *(ts_ + k) == *(qs_ + k)) // use this for generic CPUs. It is slightly faster than the unoptimized version
				++k;
		if (k + p->d == ql - 1 || k == tl - 1) {
			*t_end = k, *q_end = k + p->d;
			return -1;
		}
		p->k = k;
	}

	// wfa_next
	b[0].d = a[0].d - 1;
	b[0].p = 1;
	b[0].k = a[0].k + 1;
	b[1].d = a[0].d;
	b[1].p =  n == 1 || a[0].k > a[1].k? 0 : 1;
	b[1].k = (n == 1 || a[0].k > a[1].k? a[0].k : a[1].k) + 1;
	for (j = 1; j < n - 1; ++j) {
		int32_t k = a[j-1].k, p = -1;
		p = k > a[j+1].k + 1? p : 1;
		k = k > a[j+1].k + 1? k : a[j+1].k + 1;
		p = k > a[j].k + 1? p : 0;
		k = k > a[j].k + 1? k : a[j].k + 1;
		b[j+1].d = a[j].d, b[j+1].k = k, b[j+1].p = p;
	}
	if (n >= 2) {
		b[n].d = a[n-1].d;
		b[n].p = a[n-2].k > a[n-1].k + 1? -1 : 0;
		b[n].k = a[n-2].k > a[n-1].k + 1? a[n-2].k : a[n-1].k + 1;
	}
	b[n+1].d = a[n-1].d + 1;
	b[n+1].p = -1;
	b[n+1].k = a[n-1].k;

	// keep traceback information
	wf_tb_add(km, tb, n + 2);
	q = &tb->a[tb->n - 1];
	q->d0 = b[0].d;
	for (j = 0; j < n + 2; ++j)
		q->a[j>>5] |= (uint64_t)(b[j].p + 1) << (j&0x1f)*2;

	// drop out-of-bound cells
	memcpy(a, b, (n + 2) * sizeof(*a));
	return n + 2;
}

typedef struct {
	int32_t m, n;
	uint32_t *cigar;
} wf_cigar_t;

static void wf_cigar_push(void *km, wf_cigar_t *c, int32_t op, int32_t len)
{
	if (c->n && op == (c->cigar[c->n-1]&0xf)) {
		c->cigar[c->n-1] += len<<4;
	} else {
		if (c->n == c->m) {
			c->m = c->m + (c->m>>1) + 4;
			KREALLOC(km, c->cigar, c->m);
		}
		c->cigar[c->n++] = len<<4 | op;
	}
}

static uint32_t *wf_traceback(void *km, int32_t t_end, const char *ts, int32_t q_end, const char *qs, wf_tb_t *tb, int32_t *n_cigar)
{
	wf_cigar_t cigar = {0,0,0};
	int32_t i = q_end, k = t_end, s = tb->n - 1;
	for (;;) {
		int32_t k0 = k, j, pre;
		while (i >= 0 && k >= 0 && qs[i] == ts[k])
			--i, --k;
		if (k0 - k > 0)	
			wf_cigar_push(km, &cigar, 7, k0 - k);
		if (i < 0 || k < 0) break;
		if (s < 0) fprintf(stderr, "i=%d, k=%d, s=%d\n", i, k, tb->n);
		assert(s >= 0);
		j = i - k - tb->a[s].d0;
		assert(j < tb->a[s].n);
		pre = (int32_t)(tb->a[s].a[j>>5] >> (j&0x1f)*2 & 0x3) - 1;
		if (pre == 0) {
			if (ts[k] == qs[i]) wf_cigar_push(km, &cigar, 7, 1);
			else wf_cigar_push(km, &cigar, 8, 1);
			--i, --k;
		} else if (pre < 0) {
			wf_cigar_push(km, &cigar, 1, 1);
			--i;
		} else {
			wf_cigar_push(km, &cigar, 2, 1);
			--k;
		}
		--s;
	}
	if (i > 0) wf_cigar_push(km, &cigar, 1, i);
	else if (k > 0) wf_cigar_push(km, &cigar, 2, i);
	for (i = 0; i < cigar.n>>1; ++i) {
		uint32_t t = cigar.cigar[i];
		cigar.cigar[i] = cigar.cigar[cigar.n - i - 1];
		cigar.cigar[cigar.n - i - 1] = t;
	}
	*n_cigar = cigar.n;
	return cigar.cigar;
}

uint32_t *lv_ed_semi_cigar(void *km, int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t *score, int32_t *t_endl, int32_t *n_cigar)
{
	int32_t s = 0, n = 1, t_end = -1, q_end = -1, i;
	wf_diag_t *a;
	uint32_t *cigar;
	wf_tb_t tb = {0,0,0};
	KMALLOC(km, a, 2 * (tl + ql + 1));
	a[0].d = 0, a[0].k = -1;
	while (1) {
		n = wf_step(km, &tb, tl, ts, ql, qs, n, a, &t_end, &q_end);
		if (n < 0) break;
		++s;
	}
	kfree(km, a);
	cigar = wf_traceback(km, t_end, ts, q_end, qs, &tb, n_cigar);
	for (i = 0; i < tb.n; ++i) kfree(km, tb.a[i].a);
	kfree(km, tb.a);
	*score = s, *t_endl = t_end + 1;
	return cigar;
}
