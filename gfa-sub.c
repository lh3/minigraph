#include <assert.h>
#include <stdio.h>
#include "gfa.h"
#include "kalloc.h"
#include "kavl.h"
#include "khash.h"

typedef struct tnode_s {
	uint64_t nd;
	uint32_t v, in_tree;
	KAVL_HEAD(struct tnode_s) head;
} tnode_t;

typedef tnode_t *tnode_p;

#define tn_lt(a, b) ((a)->nd < (b)->nd || ((a)->nd == (b)->nd && (a)->v < (b)->v))
#define tn_cmp(a, b) (tn_lt(b, a) - tn_lt(a, b))

KAVL_INIT(v, tnode_t, head, tn_cmp)
KHASH_MAP_INIT_INT(v, tnode_p)

static inline tnode_t *gen_tnode(void *km, const gfa_t *g, uint32_t v, int32_t d)
{
	tnode_t *p;
	p = KMALLOC(km, tnode_t, 1);
	p->v = v, p->in_tree = 1;
	p->nd = (uint64_t)gfa_arc_n(g, v^1) << 32 | d;
	return p;
}

/* Extract a subgraph extended from a vertex within a radius. If the subgraph
 * is DAG, vertices are in the topological sorting order. The algorithm is
 * modified from Kahn's algorithm.
 */
gfa_sub_t *gfa_sub_from(void *km0, const gfa_t *g, uint32_t v0, int32_t max_dist)
{
	void *km;
	tnode_t *p, *root = 0, **L = 0;
	khash_t(v) *h;
	khint_t k;
	int32_t j, n_L = 0, m_L = 0, n_arc = 0, off;
	int absent;
	gfa_sub_t *sub = 0;

	km = km_init2(km0, 0x10000);
	h = kh_init2(v, km);

	k = kh_put(v, h, v0, &absent);
	p = kh_val(h, k) = gen_tnode(km, g, v0, 0);
	kavl_insert(v, &root, p, 0);

	while (kavl_size(head, root) > 0) {
		tnode_t *q;
		int32_t i, nv, d;
		gfa_arc_t *av;

		q = kavl_erase_first(v, &root); // take out the "smallest" vertex
		q->in_tree = 0;
		if (n_L == m_L) KEXPAND(km, L, m_L);
		L[n_L++] = q;

		d = (uint32_t)q->nd;
		nv = gfa_arc_n(g, q->v);
		av = gfa_arc_a(g, q->v);
		for (i = 0; i < nv; ++i) {
			gfa_arc_t *avi = &av[i];
			int32_t dt = d + (uint32_t)avi->v_lv;
			if (dt > max_dist) continue;
			++n_arc;
			k = kh_put(v, h, avi->w, &absent);
			if (absent) { // a vertex that hasn't been visited before
				p = kh_val(h, k) = gen_tnode(km, g, avi->w, dt);
			} else { // visited before; then update the info
				p = kh_val(h, k);
				if (!p->in_tree) continue; // when there is a cycle, a vertex may be added to L[] already
				kavl_erase(v, &root, p, 0);
				if (dt < (uint32_t)p->nd)
					p->nd = p->nd>>32<<32 | dt;
			}
			assert(p->nd>>32 > 0);
			p->nd -= 1ULL<<32;
			kavl_insert(v, &root, p, 0); // insert/re-insert to the tree
		}
	}
	assert(kh_size(h) == n_L);

	sub = KCALLOC(km0, gfa_sub_t, 1);
	sub->km = km0;
	sub->n_v = n_L;
	sub->n_a = n_arc;
	sub->v = KCALLOC(sub->km, gfa_subv_t, n_L);
	sub->a = KCALLOC(sub->km, int32_t, n_arc);

	for (j = 0; j < n_L; ++j) L[j]->in_tree = j;
	for (j = 0, off = 0; j < sub->n_v; ++j) {
		int32_t i, nv, o0 = off;
		gfa_arc_t *av;
		nv = gfa_arc_n(g, L[j]->v);
		av = gfa_arc_a(g, L[j]->v);
		for (i = 0; i < nv; ++i) {
			gfa_arc_t *avi = &av[i];
			k = kh_get(v, h, avi->w);
			if (k == kh_end(h)) continue;
			sub->a[off++] = kh_val(h, k)->in_tree;
		}
		sub->v[j].v = L[j]->v;
		sub->v[j].off = o0;
		sub->v[j].n = off - o0;
	}
	assert(off == n_arc);

	km_destroy(km);
	return sub;
}

void gfa_sub_destroy(gfa_sub_t *sub)
{
	void *km;
	if (sub == 0) return;
	km = sub->km;
	kfree(km, sub->v); kfree(km, sub->a); kfree(km, sub);
}

void gfa_sub_print(FILE *fp, const gfa_t *g, const gfa_sub_t *sub)
{
	int32_t i, j;
	for (i = 0; i < sub->n_v; ++i) {
		gfa_subv_t *p = &sub->v[i];
		fprintf(fp, "[%d]\t%d\t%s\t%d", i, p->v, g->seg[p->v>>1].name, p->n);
		if (p->n > 0) {
			fputc('\t', fp);
			for (j = 0; j < p->n; ++j) {
				if (j) fputc(',', fp);
				fprintf(fp, "%d", sub->a[p->off + j]);
			}
		}
		fputc('\n', fp);
	}
}
