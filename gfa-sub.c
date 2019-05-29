#include <assert.h>
#include <stdio.h>
#include "gfa.h"
#include "kalloc.h"
#include "kavl.h"
#include "khash.h"
#include "ksort.h"

/*********************************************
 * Extract a subgraph starting from a vertex *
 *********************************************/

#define generic_key(x) (x)
KRADIX_SORT_INIT(32, int32_t, generic_key, 4)

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
	sub->is_dag = 1;

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
		sub->v[j].d = (uint32_t)L[j]->nd;
		sub->v[j].off = o0;
		sub->v[j].n = off - o0;
		radix_sort_32(&sub->a[o0], &sub->a[off]);
		if (sub->a[o0] <= j) sub->is_dag = 0;
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
		fprintf(fp, "[%d]\t%d\t%s\t%d\t%d", i, p->v, g->seg[p->v>>1].name, p->d, p->n);
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

/********************
 * k shortest paths *
 ********************/

typedef struct sp_node_s {
	uint64_t di; // dist<<32 | unique_id
	uint32_t v;
	int32_t pre;
	KAVL_HEAD(struct sp_node_s) head;
} sp_node_t, *sp_node_p;

#define sp_node_cmp(a, b) (((a)->di > (b)->di) - ((a)->di < (b)->di))
KAVL_INIT(sp, sp_node_t, head, sp_node_cmp)

#define sp_node_lt(a, b) ((a)->di < (b)->di)
KSORT_INIT(sp, sp_node_p, sp_node_lt)

typedef struct {
	int32_t k;
	sp_node_t *p[GFA_MAX_SHORT_K]; // this forms a max-heap
} sp_topk_t;

KHASH_MAP_INIT_INT(sp, sp_topk_t)

static inline sp_node_t *gen_sp_node(void *km, const gfa_t *g, uint32_t v, int32_t d, int32_t id)
{
	sp_node_t *p;
	p = KMALLOC(km, sp_node_t, 1);
	p->v = v, p->di = (uint64_t)d<<32 | id<<1, p->pre = -1;
	return p;
}

gfa_pathv_t *gfa_sub_shortest_k(void *km0, const gfa_t *g, uint32_t src, uint32_t dst, int32_t max_dist, int32_t max_k, int32_t target_dist, int32_t *n_pathv)
{
	sp_node_t *p, *root = 0, **out, *out_dst[GFA_MAX_SHORT_K];
	sp_topk_t *q;
	khash_t(sp) *h;
	void *km;
	khint_t k;
	int absent;
	uint32_t id, n_dst, n_out, m_out;
	gfa_pathv_t *ret = 0;

	*n_pathv = 0;
	if (max_k > GFA_MAX_SHORT_K) max_k = GFA_MAX_SHORT_K;
	km = km_init2(km0, 0x10000);
	h = kh_init2(sp, km);
	kh_resize(sp, h, 16);
	m_out = 16, n_out = 0;
	out = KMALLOC(km, sp_node_t*, m_out);

	id = 0;
	p = gen_sp_node(km, g, src, 0, id++);
	kavl_insert(sp, &root, p, 0);
	k = kh_put(sp, h, src, &absent);
	q = &kh_val(h, k);
	q->k = 1, q->p[0] = p;

	n_dst = 0;
	while (kavl_size(head, root) > 0) {
		int32_t i, nv;
		gfa_arc_t *av;
		sp_node_t *r;

		r = kavl_erase_first(sp, &root);
		//fprintf(stderr, "*** %u, %s\n", kavl_size(head, root), g->seg[r->v>>1].name);
		if (n_out == m_out) KEXPAND(km, out, m_out);
		r->di = r->di>>32<<32 | n_out; // lower 32 bits now for position in the out[] array
		out[n_out++] = r;

		if (r->v == dst) { // reached the dst vertex
			out_dst[n_dst++] = r;
			if (n_dst >= max_k) break;
			if (target_dist >= 0 && r->di>>32 >= target_dist)
				break;
		}

		nv = gfa_arc_n(g, r->v);
		av = gfa_arc_a(g, r->v);
		for (i = 0; i < nv; ++i) {
			gfa_arc_t *ai = &av[i];
			int32_t d = (r->di>>32) + (uint32_t)ai->v_lv;
			if (d > max_dist) continue;
			k = kh_put(sp, h, ai->w, &absent);
			q = &kh_val(h, k);
			if (absent) q->k = 0;
			if (q->k < max_k) { // enough room
				p = gen_sp_node(km, g, ai->w, d, id++);
				p->pre = n_out - 1;
				kavl_insert(sp, &root, p, 0);
				q->p[q->k++] = p;
				ks_heapup_sp(q->k, q->p);
			} else if (q->p[0]->di>>32 > d) { // update the longest
				p = kavl_erase(sp, &root, q->p[0], 0);
				assert(p);
				p->di = (uint64_t)d<<32 | id++;
				p->pre = n_out - 1;
				ks_heapdown_sp(0, q->k, q->p);
				kavl_insert(sp, &root, p, 0);
			}
		}
	}

	if (n_dst > 0) {
		if (target_dist < 0) {
			int32_t i;
			*n_pathv = n_out;
			ret = KMALLOC(km0, gfa_pathv_t, n_out);
			for (i = 0; i < n_out; ++i)
				ret[i].v = out[i]->v, ret[i].d = out[i]->di>>32, ret[i].pre = out[i]->pre;
		} else {
			int32_t i, min_i = -1, min_diff = 0x7fffffff, n;
			for (i = 0; i < n_dst; ++i) {
				int32_t d = out_dst[i]->di >> 32;
				int32_t diff = d > target_dist? d - target_dist : target_dist - d;
				if (diff < min_diff) min_diff = diff, min_i = i;
			}
			assert(min_i >= 0);
			for (i = (int32_t)out_dst[min_i]->di, n = 0; i >= 0; i = (int32_t)out[i]->pre)
				++n;
			*n_pathv = n;
			ret = KMALLOC(km0, gfa_pathv_t, n_out);
			for (i = (int32_t)out_dst[min_i]->di; i >= 0; i = (int32_t)out[i]->pre) {
				--n;
				ret[n].v = out[i]->v, ret[n].d = out[i]->di>>32, ret[n].pre = n - 1;
			}
		}
	}

	km_destroy(km);
	return ret;
}

void gfa_sub_print_path(FILE *fp, const gfa_t *g, int32_t n, gfa_pathv_t *path)
{
	int32_t i;
	for (i = 0; i < n; ++i) {
		gfa_pathv_t *p = &path[i];
		fprintf(fp, "[%d]\t%d\t%s\t%d\t%d\n", i, p->v, g->seg[p->v>>1].name, p->d, p->pre);
	}
}
