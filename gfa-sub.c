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
KRADIX_SORT_INIT(gfa32, int32_t, generic_key, 4)

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
	KMALLOC(km, p, 1);
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

	KCALLOC(km0, sub, 1);
	sub->km = km0;
	sub->n_v = n_L;
	sub->n_a = n_arc;
	KCALLOC(sub->km, sub->v, n_L);
	KCALLOC(sub->km, sub->a, n_arc);
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
		radix_sort_gfa32(&sub->a[o0], &sub->a[off]);
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

#define generic_key(x) (x)
KRADIX_SORT_INIT(gfa64, uint64_t, generic_key, 8)

typedef struct {
	int32_t k;
	sp_node_t *p[GFA_MAX_SHORT_K]; // this forms a max-heap
} sp_topk_t;

KHASH_MAP_INIT_INT(sp, sp_topk_t)
KHASH_MAP_INIT_INT(sp2, uint64_t)

static inline sp_node_t *gen_sp_node(void *km, const gfa_t *g, uint32_t v, int32_t d, int32_t id)
{
	sp_node_t *p;
	KMALLOC(km, p, 1);
	p->v = v, p->di = (uint64_t)d<<32 | id, p->pre = -1;
	return p;
}

gfa_pathv_t *gfa_shortest_k(void *km0, const gfa_t *g, uint32_t src, int32_t n_dst, gfa_path_dst_t *dst, int32_t max_dist, int32_t max_k, int32_t *n_pathv)
{
	sp_node_t *p, *root = 0, **out;
	sp_topk_t *q;
	khash_t(sp) *h;
	khash_t(sp2) *h2;
	void *km;
	khint_t k;
	int absent;
	int32_t i, j, n_finished, n_found;
	uint32_t id, n_out, m_out;
	int8_t *dst_finish;
	gfa_pathv_t *ret = 0;
	uint64_t *dst_group;

	if (n_pathv) *n_pathv = 0;
	if (n_dst <= 0) return 0;
	for (i = 0; i < n_dst; ++i)
		dst[i].dist = -1, dst[i].n_path = 0, dst[i].path_end = -1;
	if (max_k > GFA_MAX_SHORT_K) max_k = GFA_MAX_SHORT_K;
	km = km_init2(km0, 0x10000);

	KCALLOC(km, dst_finish, n_dst);
	KMALLOC(km, dst_group, n_dst);
	for (i = 0; i < n_dst; ++i) // multiple dst[] may have the same dst[].v. We need to group them first.
		dst_group[i] = (uint64_t)dst[i].v<<32 | i;
	radix_sort_gfa64(dst_group, dst_group + n_dst);

	h2 = kh_init2(sp2, km); // this hash table keeps all destinations
	kh_resize(sp2, h2, n_dst * 2);
	for (i = 1, j = 0; i <= n_dst; ++i) {
		if (i == n_dst || dst_group[i]>>32 != dst_group[j]>>32) {
			k = kh_put(sp2, h2, dst_group[j]>>32, &absent);
			kh_val(h2, k) = (uint64_t)j << 32 | (i - j);
			assert(absent);
			j = i;
		}
	}

	h = kh_init2(sp, km); // this hash table keeps visited vertices
	kh_resize(sp, h, 16);
	m_out = 16, n_out = 0;
	KMALLOC(km, out, m_out);

	id = 0;
	p = gen_sp_node(km, g, src, 0, id++);
	kavl_insert(sp, &root, p, 0);
	k = kh_put(sp, h, src, &absent);
	q = &kh_val(h, k);
	q->k = 1, q->p[0] = p;

	n_finished = 0;
	while (kavl_size(head, root) > 0) {
		int32_t i, nv;
		gfa_arc_t *av;
		sp_node_t *r;

		r = kavl_erase_first(sp, &root); // take out the closest vertex in the heap (as a binary tree)
		//fprintf(stderr, "XX\t%d\t%d\t%d\t%c%s[%d]\t%d\n", n_out, kavl_size(head, root), n_finished, "><"[(r->v&1)^1], g->seg[r->v>>1].name, r->v, (int32_t)(r->di>>32));
		if (n_out == m_out) KEXPAND(km, out, m_out);
		r->di = r->di>>32<<32 | n_out; // lower 32 bits now for position in the out[] array
		out[n_out++] = r;

		k = kh_get(sp2, h2, r->v);
		if (k != kh_end(h2)) { // we have reached one dst vertex
			int32_t finished = 0, j;
			int32_t off = kh_val(h2, k) >> 32, cnt = (int32_t)kh_val(h2, k);
			for (j = 0; j < cnt; ++j) {
				gfa_path_dst_t *t = &dst[(int32_t)dst_group[off + j]];
				if (t->n_path == 0) { // TODO: when there is only one path, but the distance is smaller than target_dist, the dst won't be finished
					t->path_end = n_out - 1;
				} else if (t->target_dist >= 0) { // we have a target distance; choose the closest
					int32_t d0 = out[t->path_end]->di >> 32, d1 = r->di >> 32;
					d0 = d0 > t->target_dist? d0 - t->target_dist : t->target_dist - d0;
					d1 = d1 > t->target_dist? d1 - t->target_dist : t->target_dist - d1;
					if (d1 < d0) t->path_end = n_out - 1;
				}
				if (t->target_dist >= 0 && r->di>>32 >= t->target_dist) finished = 1;
				++t->n_path;
				if (t->n_path >= max_k) finished = 1;
				if (dst_finish[j] == 0 && finished)
					dst_finish[j] = 1, ++n_finished;
				if (n_finished == n_dst) break;
			}
		}

		nv = gfa_arc_n(g, r->v);
		av = gfa_arc_a(g, r->v);
		for (i = 0; i < nv; ++i) { // visit all neighbors
			gfa_arc_t *ai = &av[i];
			int32_t d = (r->di>>32) + (uint32_t)ai->v_lv;
			if (d > max_dist) continue; // don't probe vertices too far away
			k = kh_put(sp, h, ai->w, &absent);
			q = &kh_val(h, k);
			if (absent) q->k = 0;
			if (q->k < max_k) { // enough room: add to the heap
				p = gen_sp_node(km, g, ai->w, d, id++);
				p->pre = n_out - 1;
				kavl_insert(sp, &root, p, 0);
				q->p[q->k++] = p;
				ks_heapup_sp(q->k, q->p);
			} else if (q->p[0]->di>>32 > d) { // shorter than the longest path so far: replace the longest (TODO: this block is not well tested)
				p = kavl_erase(sp, &root, q->p[0], 0);
				if (p) {
					p->di = (uint64_t)d<<32 | (id++);
					p->pre = n_out - 1;
					kavl_insert(sp, &root, p, 0);
					ks_heapdown_sp(0, q->k, q->p);
				} else {
					fprintf(stderr, "Warning: logical bug in gfa_shortest_k(): q->k=%d,q->p[0]->{d,i}={%d,%d},d=%d,src=%u,max_dist=%d,n_dst=%d\n", q->k, (int32_t)(q->p[0]->di>>32), (int32_t)q->p[0]->di, d, src, max_dist, n_dst);
					km_destroy(km);
					return 0;
				}
			} // else: the path is longer than all the existing paths ended at ai->w
		}
	}

	n_found = 0;
	for (i = 0; i < n_dst; ++i) {
		gfa_path_dst_t *t = &dst[i];
		t->dist = t->n_path > 0? out[t->path_end]->di>>32 : -1;
		if (t->n_path) ++n_found;
	}

	if (n_found > 0 && n_pathv) { // then generate the backtrack array
		int32_t n, *trans;

		kh_destroy(sp, h);
		kfree(km, dst_finish);

		KCALLOC(km, trans, n_out); // used to squeeze unused elements in out[]
		for (i = 0; i < n_dst; ++i) { // mark dst vertices with a target distance
			gfa_path_dst_t *t = &dst[i];
			if (t->n_path > 0 && t->target_dist >= 0)
				trans[(int32_t)out[t->path_end]->di] = 1;
		}
		for (i = 0; i < n_out; ++i) { // mark dst vertices without a target distance
			k = kh_get(sp2, h2, out[i]->v);
			if (k != kh_end(h2)) { // TODO: check if this is correct!
				int32_t off = kh_val(h2, k)>>32, cnt = (int32_t)kh_val(h2, k);
				for (j = off; j < off + cnt; ++j)
					if (dst[j].target_dist < 0)
						trans[i] = 1;
			}
		}
		for (i = n_out - 1; i >= 0; --i) // mark all predecessors
			if (trans[i] && out[i]->pre >= 0)
				trans[out[i]->pre] = 1;
		for (i = n = 0; i < n_out; ++i) // generate coordinate translations
			if (trans[i]) trans[i] = n++;
			else trans[i] = -1;

		*n_pathv = n;
		KMALLOC(km0, ret, n);
		for (i = 0; i < n_out; ++i) { // generate the backtrack array
			gfa_pathv_t *p;
			if (trans[i] < 0) continue;
			p = &ret[trans[i]];
			p->v = out[i]->v, p->d = out[i]->di >> 32;
			p->pre = out[i]->pre < 0? out[i]->pre : trans[out[i]->pre];
		}
		for (i = 0; i < n_dst; ++i) // translate "path_end"
			if (dst[i].path_end >= 0)
				dst[i].path_end = trans[dst[i].path_end];
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
