#include <assert.h>
#include "mgpriv.h"
#include "khashl.h"
#include "kthread.h"
#include "kvec-km.h"
#include "sys.h"

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASHL_MAP_INIT(KH_LOCAL, idxhash_t, mg_hidx, uint64_t, uint64_t, idx_hash, idx_eq)

typedef struct mg_idx_bucket_s {
	mg128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mg_idx_bucket_t;

mg_idx_t *mg_idx_init(int k, int w, int b)
{
	mg_idx_t *gi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	KCALLOC(0, gi, 1);
	gi->w = w, gi->k = k, gi->b = b;
	KCALLOC(0, gi->B, 1<<b);
	return gi;
}

void mg_idx_destroy(mg_idx_t *gi)
{
	uint32_t i;
	if (gi == 0) return;
	if (gi->B) {
		for (i = 0; i < 1U<<gi->b; ++i) {
			free(gi->B[i].p);
			free(gi->B[i].a.a);
			mg_hidx_destroy((idxhash_t*)gi->B[i].h);
		}
		free(gi->B);
	}
	gfa_edseq_destroy(gi->n_seg, gi->es);
	free(gi);
}

/****************
 * Index access *
 ****************/

const uint64_t *mg_idx_hget(const void *h_, const uint64_t *q, int suflen, uint64_t minier, int *n)
{
	khint_t k;
	const idxhash_t *h = (const idxhash_t*)h_;
	*n = 0;
	if (h == 0) return 0;
	k = mg_hidx_get(h, minier>>suflen<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &q[kh_val(h, k)>>32];
	}
}

const uint64_t *mg_idx_get(const mg_idx_t *gi, uint64_t minier, int *n)
{
	int mask = (1<<gi->b) - 1;
	mg_idx_bucket_t *b = &gi->B[minier&mask];
	return mg_idx_hget(b->h, b->p, gi->b, minier, n);
}

void mg_idx_cal_quantile(const mg_idx_t *gi, int32_t m, float f[], int32_t q[])
{
	int32_t i;
	uint64_t n = 0;
	khint_t *a, k;
	for (i = 0; i < 1<<gi->b; ++i)
		if (gi->B[i].h) n += kh_size((idxhash_t*)gi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = 0, n = 0; i < 1<<gi->b; ++i) {
		idxhash_t *h = (idxhash_t*)gi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	for (i = 0; i < m; ++i)
		q[i] = ks_ksmall_uint32_t(n, a, (size_t)((1.0 - (double)f[i]) * n));
	free(a);
}

/***************
 * Index build *
 ***************/

static void mg_idx_add(mg_idx_t *gi, int n, const mg128_t *a)
{
	int i, mask = (1<<gi->b) - 1;
	for (i = 0; i < n; ++i) {
		mg128_v *p = &gi->B[a[i].x>>8&mask].a;
		kv_push(mg128_t, 0, *p, a[i]);
	}
}

void mg_idx_hfree(void *h_)
{
	idxhash_t *h = (idxhash_t*)h_;
	if (h == 0) return;
	mg_hidx_destroy(h);
}

void *mg_idx_a2h(void *km, int32_t n_a, mg128_t *a, int suflen, uint64_t **q_, int32_t *n_)
{
	int32_t N, n, n_keys;
	int32_t j, start_a, start_q;
	idxhash_t *h;
	uint64_t *q;

	*q_ = 0, *n_ = 0;
	if (n_a == 0) return 0;

	// sort by minimizer
	radix_sort_128x(a, a + n_a);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, N = 0; j <= n_a; ++j) {
		if (j == n_a || a[j].x>>8 != a[j-1].x>>8) {
			++n_keys;
			if (n > 1) N += n;
			n = 1;
		} else ++n;
	}
	h = mg_hidx_init2(km);
	mg_hidx_resize(h, n_keys);
	KCALLOC(km, q, N);
	*q_ = q, *n_ = N;

	// create the hash table
	for (j = 1, n = 1, start_a = start_q = 0; j <= n_a; ++j) {
		if (j == n_a || a[j].x>>8 != a[j-1].x>>8) {
			khint_t itr;
			int absent;
			mg128_t *p = &a[j-1];
			itr = mg_hidx_put(h, p->x>>8>>suflen<<1, &absent);
			assert(absent && j == start_a + n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					q[start_q + k] = a[start_a + k].y;
				radix_sort_gfa64(&q[start_q], &q[start_q + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				kh_val(h, itr) = (uint64_t)start_q<<32 | n;
				start_q += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	assert(N == start_q);
	return h;
}

static void worker_post(void *g, long i, int tid)
{
	mg_idx_t *gi = (mg_idx_t*)g;
	mg_idx_bucket_t *b = &gi->B[i];
	if (b->a.n == 0) return;
	b->h = (idxhash_t*)mg_idx_a2h(0, b->a.n, b->a.a, gi->b, &b->p, &b->n);
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}

int mg_gfa_overlap(const gfa_t *g)
{
	int64_t i;
	for (i = 0; i < g->n_arc; ++i) // non-zero overlap
		if (g->arc[i].ov != 0 || g->arc[i].ow != 0)
			return 1;
	return 0;
}

mg_idx_t *mg_index_core(gfa_t *g, int k, int w, int b, int n_threads)
{
	mg_idx_t *gi;
	mg128_v a = {0,0,0};
	int i;

	if (mg_gfa_overlap(g)) {
		if (mg_verbose >= 1)
			fprintf(stderr, "[E::%s] minigraph doesn't work with graphs containing overlapping segments\n", __func__);
		return 0;
	}
	gi = mg_idx_init(k, w, b);
	gi->g = g;

	for (i = 0; i < g->n_seg; ++i) {
		gfa_seg_t *s = &g->seg[i];
		a.n = 0;
		mg_sketch(0, s->seq, s->len, w, k, i, &a); // TODO: this can be parallelized
		mg_idx_add(gi, a.n, a.a);
	}
	free(a.a);
	kt_for(n_threads, worker_post, gi, 1<<gi->b);
	return gi;
}

mg_idx_t *mg_index(gfa_t *g, const mg_idxopt_t *io, int n_threads, mg_mapopt_t *mo)
{
	int32_t i, j;
	mg_idx_t *gi;
	for (i = 0; i < g->n_seg; ++i) { // uppercase
		gfa_seg_t *s = &g->seg[i];
		for (j = 0; j < s->len; ++j)
			if (s->seq[j] >= 'a' && s->seq[j] <= 'z')
				s->seq[j] -= 32;
	}
	gi = mg_index_core(g, io->k, io->w, io->bucket_bits, n_threads);
	if (gi == 0) return 0;
	gi->es = gfa_edseq_init(gi->g);
	gi->n_seg = g->n_seg;
	if (mg_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] indexed the graph\n", __func__,
				realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));
	if (mo) mg_opt_update(gi, mo, 0);
	return gi;
}
