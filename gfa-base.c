#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "gfa-priv.h"
#include "kstring.h"

#include "khash.h"
KHASH_MAP_INIT_STR(seg, uint32_t)
typedef khash_t(seg) seghash_t;

#include "ksort.h"
#define gfa_arc_key(a) ((a).v_lv)
KRADIX_SORT_INIT(arc, gfa_arc_t, gfa_arc_key, 8)

#define generic_key(x) (x)
KRADIX_SORT_INIT(gfa64, uint64_t, generic_key, 8)

int gfa_verbose = 2;

gfa_t *gfa_init(void)
{
	gfa_t *g;
	g = (gfa_t*)calloc(1, sizeof(gfa_t));
	g->h_names = kh_init(seg);
	g->h_snames = kh_init(seg);
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
		free(s->seq);
		free(s->aux.aux);
		if (s->utg) {
			for (j = 0; j < s->utg->n; ++j)
				free(s->utg->name[j]);
			free(s->utg->name);
			free(s->utg->a);
			free(s->utg);
		}
	}
	for (i = 0; i < g->n_sseq; ++i) free(g->sseq[i].name);
	kh_destroy(seg, (seghash_t*)g->h_snames);
	if (g->link_aux)
		for (k = 0; k < g->n_arc; ++k)
			free(g->link_aux[k].aux);
	free(g->idx); free(g->seg); free(g->arc); free(g->link_aux); free(g->sseq);
	free(g);
}

char *gfa_strdup(const char *src)
{
	int32_t len;
	char *dst;
	len = strlen(src);
	GFA_MALLOC(dst, len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}

char *gfa_strndup(const char *src, size_t n)
{
	char *dst;
	GFA_MALLOC(dst, n + 1);
	strncpy(dst, src, n);
	dst[n] = 0;
	return dst;
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
		kh_key(h, k) = s->name = gfa_strdup(name);
		s->del = s->len = 0;
		s->snid = s->soff = s->rank = -1;
		kh_val(h, k) = g->n_seg - 1;
	}
	return kh_val(h, k);
}

int32_t gfa_sseq_add(gfa_t *g, const char *sname)
{
	khash_t(seg) *h = (khash_t(seg)*)g->h_snames;
	khint_t k;
	int absent;
	k = kh_put(seg, h, sname, &absent);
	if (absent) {
		gfa_sseq_t *ss;
		if (g->n_sseq == g->m_sseq) GFA_EXPAND(g->sseq, g->m_sseq);
		ss = &g->sseq[g->n_sseq++];
		kh_val(h, k) = g->n_sseq - 1;
		kh_key(h, k) = ss->name = gfa_strdup(sname);
		ss->min = -1, ss->max = -1, ss->rank = -1;
	}
	return kh_val(h, k);
}

int32_t gfa_sseq_get(const gfa_t *g, const char *sname)
{
	khash_t(seg) *h = (khash_t(seg)*)g->h_snames;
	khint_t k;
	k = kh_get(seg, h, sname);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

void gfa_sseq_update(gfa_t *g, const gfa_seg_t *s)
{
	gfa_sseq_t *ps;
	if (s->snid < 0 || s->snid >= g->n_sseq) return;
	ps = &g->sseq[s->snid];
	if (ps->min < 0 || s->soff < ps->min) ps->min = s->soff;
	if (ps->max < 0 || s->soff + s->len > ps->max) ps->max = s->soff + s->len;
	if (ps->rank < 0) ps->rank = s->rank;
	else if (ps->rank != s->rank) {
		if (gfa_verbose >= 2)
			fprintf(stderr, "[W] stable sequence '%s' associated with different ranks on segment '%s': %d != %d\n", ps->name, s->name, ps->rank, s->rank);
	}
}

int32_t gfa_name2id(const gfa_t *g, const char *name)
{
	seghash_t *h = (seghash_t*)g->h_names;
	khint_t k;
	k = kh_get(seg, h, name);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

gfa_arc_t *gfa_add_arc1(gfa_t *g, uint32_t v, uint32_t w, int32_t ov, int32_t ow, int64_t link_id, int comp)
{
	gfa_arc_t *a;
	if (g->m_arc == g->n_arc) {
		uint64_t old_m = g->m_arc;
		g->m_arc = g->m_arc? g->m_arc<<1 : 16;
		g->arc = (gfa_arc_t*)realloc(g->arc, g->m_arc * sizeof(gfa_arc_t));
		memset(&g->arc[old_m], 0, (g->m_arc - old_m) * sizeof(gfa_arc_t));
		g->link_aux = (gfa_aux_t*)realloc(g->link_aux, g->m_arc * sizeof(gfa_aux_t));
		memset(&g->link_aux[old_m], 0, (g->m_arc - old_m) * sizeof(gfa_aux_t));
	}
	a = &g->arc[g->n_arc++];
	a->v_lv = (uint64_t)v << 32;
	a->w = w, a->ov = ov, a->ow = ow, a->rank = -1;
	a->link_id = link_id >= 0? link_id : g->n_arc - 1;
	if (link_id >= 0) a->rank = g->arc[link_id].rank; // TODO: this is not always correct!
	a->del = a->strong = 0;
	a->comp = comp;
	return a;
}

int gfa_arc_is_sorted(const gfa_t *g)
{
	uint64_t e;
	for (e = 1; e < g->n_arc; ++e)
		if (g->arc[e-1].v_lv > g->arc[e].v_lv)
			break;
	return (e == g->n_arc);
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
		const gfa_seg_t *sv = &g->seg[v>>1];
		if (!sv->del && sv->len < a->ov) {
			if (gfa_verbose >= 2)
				fprintf(stderr, "[W] overlap length longer than segment length for '%s': %d > %d\n", sv->name, a->ov, sv->len);
			a->ov = sv->len;
		}
		if (sv->del || g->seg[w>>1].del) {
			a->del = 1;
		} else {
			a->v_lv |= sv->len - a->ov;
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

uint32_t gfa_fix_symm_add(gfa_t *g)
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
				gfa_arc_t *arc_old = g->arc, *arc_new;
				arc_new = gfa_add_arc1(g, avi->w^1, v^1, avi->ow, avi->ov, avi->link_id, 1);
				if (arc_old != g->arc) av = gfa_arc_a(g, v); // g->arc may be reallocated
				arc_new->rank = av[i].rank;
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
		else {
			gfa_aux_t *aux = g->arc[e].link_id < g->n_arc? &g->link_aux[g->arc[e].link_id] : 0;
			if (aux) {
				free(aux->aux);
				aux->aux = 0, aux->l_aux = aux->m_aux = 0;
			}
		}
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
	if (!gfa_arc_is_sorted(g)) {
		gfa_arc_sort(g);
		if (g->idx) free(g->idx);
		g->idx = 0;
	}
	if (g->idx == 0) gfa_arc_index(g);
}

int32_t gfa_check_multi(const gfa_t *g)
{
	uint32_t v, n_vtx = gfa_n_vtx(g);
	int32_t max_nv = -1, n_multi = 0;
	uint64_t *buf; // actually, uint32_t is enough
	for (v = 0; v < n_vtx; ++v) {
		int32_t nv = gfa_arc_n(g, v);
		max_nv = max_nv > nv? max_nv : nv;
	}
	if (max_nv == 1 || max_nv < 0) return 0;
	GFA_MALLOC(buf, max_nv);
	for (v = 0; v < n_vtx; ++v) {
		int32_t i, s, nv = gfa_arc_n(g, v);
		const gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i) buf[i] = av[i].w;
		radix_sort_gfa64(buf, buf + nv);
		for (s = 0, i = 1; i <= nv; ++i)
			if (i == nv || buf[i] != buf[s])
				n_multi += i - s - 1, s = i;
	}
	free(buf);
	return n_multi;
}

uint32_t gfa_fix_multi(gfa_t *g)
{
	uint32_t v, n_vtx = gfa_n_vtx(g), n_rm = 0;
	int32_t max_nv = -1;
	uint64_t *buf; // actually, uint32_t is enough
	for (v = 0; v < n_vtx; ++v) {
		int32_t nv = gfa_arc_n(g, v);
		max_nv = max_nv > nv? max_nv : nv;
	}
	if (max_nv == 1) return 0;
	GFA_MALLOC(buf, max_nv);
	for (v = 0; v < n_vtx; ++v) {
		int32_t i, j, s, nv = gfa_arc_n(g, v), nb;
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = j = 0; i < nv; ++i)
			if (!av[i].del) buf[j++] = (uint64_t)av[i].w<<32 | i;
		nb = j;
		if (nb < 1) continue;
		radix_sort_gfa64(buf, buf + nb);
		for (s = 0, i = 1; i <= nb; ++i) {
			if (i == nv || buf[i]>>32 != buf[s]>>32) {
				if (i - s > 1) {
					int32_t k = (int32_t)buf[s], min_rank = av[k].rank; // prefer longest overlap
					for (j = s + 1; j < i; ++j) { // rank has higher priority
						int32_t t = (int32_t)buf[j];
						if (av[t].rank >= 0 && av[t].rank < min_rank)
							min_rank = av[t].rank, k = t;
					}
					if (av[k].w == (v^1)) { // a weird loop
						if (gfa_verbose >= 2)
							fprintf(stderr, "[W::%s] can't fix multiple edges due to '>v -- <v' involving segment %s\n", __func__, g->seg[v>>1].name);
					} else {
						int32_t nw = gfa_arc_n(g, av[k].w^1), n_wdel;
						gfa_arc_t *aw = gfa_arc_a(g, av[k].w^1);
						uint64_t link_id = av[k].link_id;
						n_rm += i - s - 1;
						for (j = s + 1; j < i; ++j)
							av[(int32_t)buf[j]].del = 1;
						for (j = 0, n_wdel = 0; j < nw; ++j)
							if (aw[j].w == (v^1) && aw[j].link_id != link_id)
								aw[j].del = 1, ++n_wdel;
						assert(n_wdel == i - s - 1);
					}
				}
				s = i;
			}
		}
	}
	free(buf);
	if (n_rm > 0) {
		if (gfa_verbose >= 2)
			fprintf(stderr, "[W::%s] removed %d multiple link(s)\n", __func__, n_rm);
		gfa_arc_rm(g);
		gfa_arc_index(g);
	}
	return n_rm;
}

void gfa_finalize(gfa_t *g)
{
	gfa_fix_no_seg(g);
	gfa_arc_sort(g);
	gfa_arc_index(g);
	gfa_fix_semi_arc(g);
	gfa_fix_symm_add(g);
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
		int type = *(s); \
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

void gfa_aux_update_f(gfa_aux_t *a, const char tag[2], float x)
{
	uint8_t *p = 0;
	if (a->l_aux > 0)
		p = gfa_aux_get(a->l_aux, a->aux, "cv");
	if (p) {
		memcpy(p + 1, &x, 4);
	} else {
		kstring_t str;
		str.l = a->l_aux, str.m = a->m_aux, str.s = (char*)a->aux;
		ks_resize(&str, str.l + 7);
		kputsn_(tag, 2, &str);
		kputc_('f', &str);
		kputsn_(&x, 4, &str);
		a->l_aux = str.l, a->m_aux = str.m, a->aux = (uint8_t*)str.s;
	}
}

void gfa_aux_update_cv(gfa_t *g, const char *tag, const double *cov_seg, const double *cov_link)
{
	int64_t i;
	if (cov_seg)
		for (i = 0; i < g->n_seg; ++i)
			gfa_aux_update_f(&g->seg[i].aux, tag, cov_seg[i]);
	if (cov_link)
		for (i = 0; i < g->n_arc; ++i)
			if (g->arc[i].comp == 0)
				gfa_aux_update_f(&g->link_aux[g->arc[i].link_id], tag, cov_link[i]);
}

/*********************
 * Translation table *
 *********************/

unsigned char gfa_comp_table[256] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 96, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};
