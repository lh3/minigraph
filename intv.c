#include <assert.h>
#include "mgpriv.h"
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

