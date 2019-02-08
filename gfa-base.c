#include <stdlib.h>
#include <ctype.h>
#include "gfa.h"

#include "khash.h"
KHASH_MAP_INIT_STR(seg, uint32_t)
typedef khash_t(seg) seghash_t;

#include "ksort.h"
#define gfa_arc_key(a) ((a).v_lv)
KRADIX_SORT_INIT(arc, gfa_arc_t, gfa_arc_key, 8)

int gfa_verbose = 2;

gfa_t *gfa_init(void)
{
	gfa_t *g;
	g = (gfa_t*)calloc(1, sizeof(gfa_t));
	g->h_names = kh_init(seg);
	return g;
}

void gfa_destroy(gfa_t *g)
{
	uint32_t i, j;
	uint64_t k;
	if (g == 0) return;
	kh_destroy(seg, (seghash_t*)g->h_names);
	for (i = 0; i < g->n_seg; ++i) {
		gfa_seg_t *s = &g->seg[i];
		free(s->name);
		free(s->aux.aux);
		for (j = 0; j < s->utg.n; ++j)
			free(s->utg.name[j]);
		free(s->utg.name);
		free(s->utg.a);
	}
	for (k = 0; k < g->n_arc; ++k)
		free(g->arc_aux[k].aux);
	free(g->idx); free(g->seg); free(g->arc); free(g->arc_aux);
	free(g);
}

int32_t gfa_add_seg(gfa_t *g, const char *name)
{
	khint_t k;
	int absent;
	seghash_t *h = (seghash_t*)g->h_names;
	k = kh_put(seg, h, name, &absent);
	if (absent) {
		gfa_seg_t *s;
		if (g->n_seg == g->m_seg) {
			uint32_t old_m = g->m_seg;
			g->m_seg = g->m_seg? g->m_seg<<1 : 16;
			g->seg = (gfa_seg_t*)realloc(g->seg, g->m_seg * sizeof(gfa_seg_t));
			memset(&g->seg[old_m], 0, (g->m_seg - old_m) * sizeof(gfa_seg_t));
		}
		s = &g->seg[g->n_seg++];
		kh_key(h, k) = s->name = strdup(name);
		s->del = s->len = 0;
		kh_val(h, k) = g->n_seg - 1;
	}
	return kh_val(h, k);
}

int32_t gfa_name2id(const gfa_t *g, const char *name)
{
	seghash_t *h = (seghash_t*)g->h_names;
	khint_t k;
	k = kh_get(seg, h, name);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

uint64_t gfa_add_arc1(gfa_t *g, uint32_t v, uint32_t w, int32_t ov, int32_t ow, int64_t link_id, int comp)
{
	gfa_arc_t *a;
	if (g->m_arc == g->n_arc) {
		uint64_t old_m = g->m_arc;
		g->m_arc = g->m_arc? g->m_arc<<1 : 16;
		g->arc = (gfa_arc_t*)realloc(g->arc, g->m_arc * sizeof(gfa_arc_t));
		memset(&g->arc[old_m], 0, (g->m_arc - old_m) * sizeof(gfa_arc_t));
		g->arc_aux = (gfa_aux_t*)realloc(g->arc_aux, g->m_arc * sizeof(gfa_aux_t));
		memset(&g->arc_aux[old_m], 0, (g->m_arc - old_m) * sizeof(gfa_aux_t));
	}
	a = &g->arc[g->n_arc++];
	a->v_lv = (uint64_t)v << 32;
	a->w = w, a->ov = ov, a->ow = ow, a->lw = 0;
	a->link_id = link_id >= 0? link_id : g->n_arc - 1;
	a->del = 0;
	a->comp = comp;
	return a->link_id;
}

void gfa_arc_sort(gfa_t *g)
{
	radix_sort_arc(g->arc, g->arc + g->n_arc);
}

uint64_t *gfa_arc_index_core(size_t max_seq, size_t n, const gfa_arc_t *a)
{
	size_t i, last;
	uint64_t *idx;
	idx = (uint64_t*)calloc(max_seq * 2, 8);
	for (i = 1, last = 0; i <= n; ++i)
		if (i == n || gfa_arc_head(a[i-1]) != gfa_arc_head(a[i]))
			idx[gfa_arc_head(a[i-1])] = (uint64_t)last<<32 | (i - last), last = i;
	return idx;
}

void gfa_arc_index(gfa_t *g)
{
	if (g->idx) free(g->idx);
	g->idx = gfa_arc_index_core(g->n_seg, g->n_arc, g->arc);
}

/********************
 * Fix graph issues *
 ********************/

uint32_t gfa_fix_no_seg(gfa_t *g)
{
	uint32_t i, n_err = 0;
	for (i = 0; i < g->n_seg; ++i) {
		gfa_seg_t *s = &g->seg[i];
		if (s->len == 0) {
			++n_err, s->del = 1;
			if (gfa_verbose >= 2)
				fprintf(stderr, "[W] segment '%s' is used on an L-line but not defined on an S-line\n", s->name);
		}
	}
	return n_err;
}

void gfa_fix_arc_len(gfa_t *g)
{
	uint64_t k;
	for (k = 0; k < g->n_arc; ++k) {
		gfa_arc_t *a = &g->arc[k];
		uint32_t v = gfa_arc_head(*a), w = gfa_arc_tail(*a);
		if (g->seg[v>>1].del || g->seg[w>>1].del) {
			a->del = 1;
		} else {
			a->v_lv |= g->seg[v>>1].len - a->ov;
			if (a->ow != INT32_MAX)
				a->lw = g->seg[w>>1].len - a->ow;
		}
	}
}

uint32_t gfa_fix_semi_arc(gfa_t *g)
{
	uint32_t n_err = 0, v, n_vtx = gfa_n_vtx(g);
	int i, j;
	for (v = 0; v < n_vtx; ++v) {
		int nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			if (!av[i].del && (av[i].ow == INT32_MAX || av[i].ov == INT32_MAX)) { // overlap length is missing
				uint32_t w = av[i].w^1;
				int is_multi = 0, c, jv = -1, nw = gfa_arc_n(g, w);
				gfa_arc_t *aw = gfa_arc_a(g, w);
				for (j = 0, c = 0; j < nw; ++j)
					if (!aw[j].del && aw[j].w == (v^1)) ++c, jv = j;
				if (c == 1) {
					if (av[i].ov != INT32_MAX && aw[jv].ow != INT32_MAX && av[i].ov != aw[jv].ow) is_multi = 1;
					if (av[i].ow != INT32_MAX && aw[jv].ov != INT32_MAX && av[i].ow != aw[jv].ov) is_multi = 1;
				}
				if (c == 1 && !is_multi) {
					if (aw[jv].ov != INT32_MAX) av[i].ow = aw[jv].ov;
					if (aw[jv].ow != INT32_MAX) av[i].ov = aw[jv].ow;
				} else {
					if (gfa_verbose >= 2)
						fprintf(stderr, "[W] can't infer overlap length for %s%c -> %s%c\n",
								g->seg[v>>1].name, "+-"[v&1], g->seg[w>>1].name, "+-"[(w^1)&1]);
					++n_err;
					av[i].del = 1;
				}
			}
		}
	}
	return n_err;
}

uint32_t gfa_fix_symm(gfa_t *g)
{
	uint32_t n_err = 0, v, n_vtx = gfa_n_vtx(g);
	int i;
	for (v = 0; v < n_vtx; ++v) {
		int nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			int j, nw;
			gfa_arc_t *aw, *avi = &av[i];
			if (avi->del || avi->comp) continue;
			nw = gfa_arc_n(g, avi->w^1);
			aw = gfa_arc_a(g, avi->w^1);
			for (j = 0; j < nw; ++j) {
				gfa_arc_t *awj = &aw[j];
				if (awj->del || awj->comp) continue;
				if (awj->w == (v^1) && awj->ov == avi->ow && awj->ow == avi->ov) { // complement found
					awj->comp = 1;
					awj->link_id = avi->link_id;
					break;
				}
			}
			if (j == nw) {
				gfa_arc_t *arc_old = g->arc;
				gfa_add_arc1(g, avi->w^1, v^1, avi->ow, avi->ov, avi->link_id, 1);
				if (arc_old != g->arc) av = gfa_arc_a(g, v); // g->arc may be reallocated
			}
		}
	}
	if (n_vtx < gfa_n_vtx(g)) {
		gfa_arc_sort(g);
		gfa_arc_index(g);
	}
	return n_err;
}

void gfa_arc_rm(gfa_t *g)
{
	uint32_t e, n;
	for (e = n = 0; e < g->n_arc; ++e) {
		uint32_t u = g->arc[e].v_lv>>32, v = g->arc[e].w;
		if (!g->arc[e].del && !g->seg[u>>1].del && !g->seg[v>>1].del)
			g->arc[n++] = g->arc[e];
	}
	if (n < g->n_arc) { // arc index is out of sync
		if (g->idx) free(g->idx);
		g->idx = 0;
	}
	g->n_arc = n;
}

void gfa_cleanup(gfa_t *g)
{
	gfa_arc_rm(g);
	if (!g->is_srt) {
		gfa_arc_sort(g);
		g->is_srt = 1;
		if (g->idx) free(g->idx);
		g->idx = 0;
	}
	if (g->idx == 0) gfa_arc_index(g);
}

void gfa_finalize(gfa_t *g)
{
	gfa_fix_no_seg(g);
	gfa_arc_sort(g);
	gfa_arc_index(g);
	gfa_fix_semi_arc(g);
	gfa_fix_symm(g);
	gfa_fix_arc_len(g);
	gfa_cleanup(g);
}

/********************
 * Tag manipulation *
 ********************/

static inline int gfa_aux_type2size(int x)
{
	if (x == 'C' || x == 'c' || x == 'A') return 1;
	else if (x == 'S' || x == 's') return 2;
	else if (x == 'I' || x == 'i' || x == 'f') return 4;
	else return 0;
}

#define __skip_tag(s) do { \
		int type = toupper(*(s)); \
		++(s); \
		if (type == 'Z') { while (*(s)) ++(s); ++(s); } \
		else if (type == 'B') (s) += 5 + gfa_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
		else (s) += gfa_aux_type2size(type); \
	} while(0)

uint8_t *gfa_aux_get(int l_data, const uint8_t *data, const char tag[2])
{
	const uint8_t *s = data;
	int y = tag[0]<<8 | tag[1];
	while (s < data + l_data) {
		int x = (int)s[0]<<8 | s[1];
		s += 2;
		if (x == y) return (uint8_t*)s;
		__skip_tag(s);
	}
	return 0;
}

// s MUST BE returned by gfa_aux_get()
int gfa_aux_del(int l_data, uint8_t *data, uint8_t *s)
{
	uint8_t *p;
	p = s - 2;
	__skip_tag(s);
	memmove(p, s, l_data - (s - data));
	return l_data - (s - p);
}
