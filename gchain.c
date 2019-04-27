#include <string.h>
#include "mgpriv.h"
#include "kvec.h"

static char comp_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

typedef struct {
	int32_t n;
	uint64_t *q;
	void *h;
} sub_idx_t;

static void mg_sub_idx(void *km, const gfa_t *g, const gfa_sub_t *s, int k, int w, int is_hpc, int32_t max_dist_x, sub_idx_t *si)
{
	int32_t i, max_len = 0, sum_len = 0;
	mg128_v a = {0,0,0};
	char *t;
	for (i = 0; i < s->n_v; ++i) {
		int32_t len;
		gfa_subv_t *v = &s->v[i];
		if (v->d > max_dist_x) continue;
		len = g->seg[v->v>>1].len;
		len = v->d + len < max_dist_x? len : max_dist_x - v->d;
		max_len = max_len > len? max_len : len;
		sum_len += len;
	}

	t = KMALLOC(km, char, max_len);
	kv_resize(mg128_t, km, a, sum_len / w * 2);
	for (i = 0; i < s->n_v; ++i) {
		int32_t len0, len;
		gfa_subv_t *v = &s->v[i];
		if (v->d > max_dist_x) continue;
		len = len0 = g->seg[v->v>>1].len;
		len = v->d + len < max_dist_x? len : max_dist_x - v->d;
		if ((v->v&1) == 0) { // forward strand
			memcpy(t, g->seg[v->v>>1].seq, len);
		} else { // reverse strand
			int32_t j, k;
			char *s0 = g->seg[v->v>>1].seq;
			for (j = len0 - 1, k = 0; j >= len0 - len; --j)
				t[k++] = (uint8_t)s0[j] < 128? comp_tab[(uint8_t)s0[j]] : s0[j];
		}
		mg_sketch(km, t, len, w, k, i, is_hpc, &a);
	}
	kfree(km, t);
	si->h = mg_idx_a2h(km, a.n, a.a, 0, &si->q, &si->n);
	kfree(km, a.a);
}

typedef struct {
	int n;
	uint32_t qpos, q_span;
	const uint64_t *r;
} sub_seeds_t;

void mg_gchain(void *km, const gfa_t *g, gfa_sub_t *s, const char *qseq, int32_t qlen, int32_t max_dist_s, int32_t max_dist_g, int k, int w, int is_hpc)
{
	int32_t i, st;
	sub_idx_t si;
	sub_seeds_t *ss;
	mg128_v qm = {0,0,0};

	mg_sketch(km, qseq, qlen, w, k, 0, is_hpc, &qm);
	mg_sub_idx(km, g, s, k, w, is_hpc, max_dist_g, &si);
	ss = KCALLOC(km, sub_seeds_t, qm.n);

	for (i = 0; i < qm.n; ++i) {
		mg128_t *p = &qm.a[i];
		sub_seeds_t *q = &ss[i];
		q->qpos = (uint32_t)p->y, q->q_span = p->x & 0xff;
		q->r = mg_idx_hget(si.h, si.q, 0, p->x>>8, &q->n);
	}

	for (i = st = 0; i < qm.n; ++i) {
		sub_seeds_t *qi = &ss[i];
		int32_t j, ki, kj;
		while (st < i && qi->qpos > ss[st].qpos + max_dist_s) ++st;
		for (ki = 0; ki < qi->n; ++ki) {
			for (j = i - 1; j >= st; --j) {
				sub_seeds_t *qj = &ss[j];
				for (kj = 0; kj < qj->n; ++kj) {
					int32_t dg, ds = qj->qpos - qi->qpos;
				}
			}
		}
	}

	// free
	kfree(km, qm.a);
	kfree(km, si.q);
	// destroy si->h!
	kfree(km, ss);
}
