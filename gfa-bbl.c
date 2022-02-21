#include <assert.h>
#include <stdio.h>
#include "gfa-priv.h"
#include "kalloc.h"
#include "ksort.h"
#include "kvec.h"

#define generic_key(x) (x)
KRADIX_SORT_INIT(gfa32, uint32_t, generic_key, 4)

void gfa_sort_ref_arc(gfa_t *g)
{
	uint32_t v, n_vtx = gfa_n_vtx(g);
	for (v = 0; v < n_vtx; ++v) {
		gfa_seg_t *s = &g->seg[v>>1];
		int32_t i, nv;
		gfa_arc_t *av, b;
		if (s->rank != 0) continue;
		nv = gfa_arc_n(g, v);
		av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			uint32_t w = av[i].w;
			gfa_seg_t *t = &g->seg[w>>1];
			if (t->rank == 0 && t->snid == s->snid && (v&1) == (w&1)) {
				if (((v&1) == 0 && s->soff + s->len == t->soff) || ((v&1) == 1 && t->soff + t->len == s->soff))
					break;
			}
		}
		if (nv > 0 && i == nv) fprintf(stderr, "X\t%c%s\t%d\t%s\t%d\n", "><"[v&1], s->name, i, g->sseq[s->snid].name, s->soff);
		assert(nv == 0 || i < nv);
		if (i > 0 && i < nv) b = av[i], av[i] = av[0], av[0] = b;
	}
}

void gfa_sub_print(FILE *fp, const gfa_t *g, const gfa_sub_t *sub)
{
	int32_t i, j;
	for (i = 0; i < sub->n_v; ++i) {
		gfa_subv_t *p = &sub->v[i];
		fprintf(fp, "[%d]\t%d\t%c%s\t%d\t%d", i, p->v, "><"[p->v&1], g->seg[p->v>>1].name, p->d, p->n);
		if (p->n > 0) {
			fputc('\t', fp);
			for (j = 0; j < p->n; ++j) {
				if (j) fputc(',', fp);
				fprintf(fp, "%d", (uint32_t)(sub->a[p->off + j]>>32));
			}
		}
		fputc('\n', fp);
	}
}

/****************
 * Tarjan's SCC *
 ****************/

typedef struct {
	uint32_t index, low:31, stack:1;
	uint32_t i;     // index in gfa_sub_t::v[]; a temporary field
	uint32_t start; // starting vertex
} gfa_scinfo_t;

struct gfa_scbuf_s {
	uint32_t index;
	gfa_scinfo_t *a;     // node information
	kvec_t(uint32_t) ts; // Tarjan's stack
	kvec_t(uint64_t) ds; // DFS stack
};

gfa_scbuf_t *gfa_scbuf_init(const gfa_t *g)
{
	uint32_t v, n_vtx = gfa_n_vtx(g);
	gfa_scbuf_t *b;
	GFA_CALLOC(b, 1);
	GFA_CALLOC(b->a, n_vtx);
	for (v = 0; v < n_vtx; ++v)
		b->a[v].index = b->a[v].start = (uint32_t)-1;
	return b;
}

void gfa_scbuf_destroy(gfa_scbuf_t *b)
{
	free(b->a); free(b->ts.a); free(b->ds.a); free(b);
}

gfa_sub_t *gfa_scc1(void *km0, const gfa_t *g, gfa_scbuf_t *b, uint32_t v0)
{
	gfa_sub_t *sub;
	uint32_t k, off, m_v = 0;

	KCALLOC(km0, sub, 1);
	sub->km = km0;

	kv_push(uint64_t, b->ds, (uint64_t)v0<<32);
	while (b->ds.n > 0) {
		uint64_t x = kv_pop(b->ds);
		uint32_t i = (uint32_t)x, v = x>>32, nv;
		if (i == 0) { // i is the number of outgoing edges already visited
			b->a[v].low = b->a[v].index = b->index++;
			b->a[v].stack = 1;
			kv_push(uint32_t, b->ts, v);
		}
		nv = gfa_arc_n(g, v);
		if (i == nv) { // done with v
			if (b->a[v].low == b->a[v].index) {
				int32_t i, j = b->ts.n - 1;
				while (b->ts.a[j] != v) --j;
				for (i = b->ts.n - 1; i >= j; --i) {
					uint32_t w = b->ts.a[i];
					gfa_subv_t *p;
					//fprintf(stderr, "V\t%c%s\t%d\t%c%s\t%d\t%d\n", "><"[v&1], g->seg[v>>1].name, i, "><"[w&1], g->seg[w>>1].name, b->a[w^1].stack, b->a[w].index);
					if (sub->n_v == m_v) KEXPAND(sub->km, sub->v, m_v);
					p = &sub->v[sub->n_v++];
					p->v = w;
					b->a[w].stack = 0;
				}
				b->ts.n = j;
			}
			if (b->ds.n > 0) { // if the DFS stack is not empty, update the top element
				uint32_t w = v;
				v = b->ds.a[b->ds.n - 1] >> 32;
				b->a[v].low = b->a[v].low < b->a[w].low? b->a[v].low : b->a[w].low;
			}
		} else { // process v's neighbor av[i].w
			gfa_arc_t *av = gfa_arc_a(g, v);
			uint32_t w = av[i].w;
			kv_push(uint64_t, b->ds, (uint64_t)v<<32 | (i+1)); // update the old top of the stack
			if (b->a[w].index == (uint32_t)-1 && b->a[w^1].stack == 0)
				kv_push(uint64_t, b->ds, (uint64_t)w<<32);
			else if (b->a[w].stack)
				b->a[v].low = b->a[v].low < b->a[w].index? b->a[v].low : b->a[w].index;
		}
	}

	// reverse the vertex array
	for (k = 0; k < sub->n_v>>1; ++k) {
		gfa_subv_t x;
		x = sub->v[k], sub->v[k] = sub->v[sub->n_v - k - 1], sub->v[sub->n_v - k - 1] = x;
	}

	// fill other fields in sub
	for (k = 0; k < sub->n_v; ++k)
		b->a[sub->v[k].v].start = v0, b->a[sub->v[k].v].i = k;
	for (k = 0, off = 0; k < sub->n_v; ++k) { // precompute the length of gfa_sub_t::a[]
		uint32_t v = sub->v[k].v;
		int32_t i, nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i)
			if (b->a[av[i].w].start == v0)
				++off;
	}
	sub->n_a = off;
	KCALLOC(sub->km, sub->a, sub->n_a);
	for (k = 0, off = 0; k < sub->n_v; ++k) {
		uint32_t o0, v = sub->v[k].v;
		int32_t i, nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0, o0 = off; i < nv; ++i)
			if (b->a[av[i].w].start == v0)
				sub->a[off++] = (uint64_t)b->a[av[i].w].i << 32 | (&av[i] - g->arc);
		sub->v[k].d = 0;
		sub->v[k].off = o0;
		sub->v[k].n = off - o0;
		if (o0 < off) {
			radix_sort_gfa64(&sub->a[o0], &sub->a[off]);
			if (sub->a[o0]>>32 <= k) sub->is_dag = 0;
		}
	}
	return sub;
}

void gfa_scc_all(const gfa_t *g)
{
	uint32_t v, n_vtx = gfa_n_vtx(g);
	gfa_scbuf_t *b;
	b = gfa_scbuf_init(g);
	for (v = 0; v < n_vtx; ++v)
		if (b->a[v].index == (uint32_t)-1 && b->a[v^1].index == (uint32_t)-1) {
			gfa_sub_t *sub;
			sub = gfa_scc1(0, g, b, v);
			gfa_sub_print(stderr, g, sub);
			gfa_sub_destroy(sub);
		}
	gfa_scbuf_destroy(b);
}

void gfa_sub_destroy(gfa_sub_t *sub)
{
	void *km;
	if (sub == 0) return;
	km = sub->km;
	kfree(km, sub->v); kfree(km, sub->a); kfree(km, sub);
}

/******************
 * Bubble calling *
 ******************/

typedef struct {
	int32_t ld, sd, rd;
	int32_t lp, sp;
	float lf, sf, rf;
} bb_aux_t;

static void bb_write_seq(const gfa_t *g, int32_t n, const uint32_t *v, int32_t l_seq, char *seq)
{
	int32_t k, l;
	for (k = n - 1, l = 0; k >= 0; --k) {
		const gfa_seg_t *s = &g->seg[v[k]>>1];
		if (v[k]&1) {
			int32_t p;
			for (p = s->len - 1; p >= 0; --p)
				seq[l++] = gfa_comp_table[(uint8_t)s->seq[p]];
		} else {
			memcpy(&seq[l], s->seq, s->len);
			l += s->len;
		}
	}
	assert(l == l_seq);
	seq[l] = 0;
}

static int32_t bb_n_paths(const gfa_t *g, const gfa_sub_t *sub, int32_t js, int32_t je)
{
	int32_t j, k;
	int64_t *cnt, c;
	GFA_CALLOC(cnt, je - js + 1);
	cnt[0] = 1;
	for (j = js; j < je; ++j) {
		const gfa_subv_t *t = &sub->v[j];
		for (k = 0; k < t->n; ++k) {
			uint64_t a = sub->a[t->off + k];
			int32_t jv = (int32_t)(a>>32);
			if (jv <= j || jv > je) continue;
			if (cnt[jv - js] + cnt[j - js] > INT32_MAX)
				cnt[jv - js] = INT32_MAX;
			else cnt[jv - js] += cnt[j - js];
		}
	}
	c = cnt[je - js];
	free(cnt);
	return c < INT32_MAX? c : INT32_MAX;
}

gfa_bubble_t *gfa_bubble(const gfa_t *g, int32_t *n_bb_)
{
	uint32_t i, *vs, *vmin, *vtmp = 0;
	int32_t n_bb = 0, m_bb = 0, m_vtmp = 0;
	gfa_bubble_t *bb = 0;
	gfa_scbuf_t *scbuf;

	GFA_MALLOC(vs, g->n_sseq);
	GFA_MALLOC(vmin, g->n_sseq);
	for (i = 0; i < g->n_sseq; ++i)
		vs[i] = (uint32_t)-1, vmin[i] = UINT32_MAX;
	for (i = 0; i < g->n_seg; ++i) {
		const gfa_seg_t *s = &g->seg[i];
		if (s->rank != 0 || s->snid < 0) continue;
		if ((uint32_t)s->soff < vmin[s->snid])
			vmin[s->snid] = s->soff, vs[s->snid] = i<<1;
	}
	free(vmin);

	scbuf = gfa_scbuf_init(g);
	for (i = 0; i < g->n_sseq; ++i) {
		gfa_sub_t *sub;
		int32_t j, jst, max_a, max_soff;
		bb_aux_t *ba;

		if (vs[i] == (uint32_t)-1) continue;
		#if 0
		sub = gfa_sub_from(0, g, vs[i], 0);
		#else
		sub = gfa_scc1(0, g, scbuf, vs[i]);
		#endif
		//gfa_sub_print(stderr, g, sub);
		GFA_CALLOC(ba, sub->n_v);
		for (j = 0; j < sub->n_v; ++j)
			ba[j].sd = INT32_MAX, ba[j].lp = ba[j].sp = -1;
		ba[0].sd = 0;
		for (j = 0; j < sub->n_v; ++j) {
			gfa_subv_t *t = &sub->v[j];
			int32_t k;
			for (k = 0; k < t->n; ++k) {
				uint64_t a = sub->a[t->off + k];
				int32_t jv = (int32_t)(a>>32);
				int32_t l = (int32_t)g->arc[(uint32_t)a].v_lv;
				if (jv <= j) continue; // skip loop or cycle
				if (ba[jv].sd >= ba[j].sd + l)
					ba[jv].sd = ba[j].sd + l, ba[jv].sp = j;
				if (ba[jv].ld < ba[j].ld + l)
					ba[jv].ld = ba[j].ld + l, ba[jv].lp = j;
			}
		}
		for (j = 0, jst = 0, max_a = max_soff = -1; j < sub->n_v; ++j) {
			gfa_subv_t *t = &sub->v[j];
			int32_t k;
			if (j == max_a && g->seg[t->v>>1].soff > max_soff) {
				const gfa_seg_t *sst = &g->seg[sub->v[jst].v>>1];
				const gfa_seg_t *sen = &g->seg[t->v>>1];
				if (sst->snid == i && sen->snid == i) {
					int32_t n, l;
					uint32_t *v;
					gfa_bubble_t *b;

					// basic information
					if (n_bb == m_bb) GFA_EXPAND(bb, m_bb);
					b = &bb[n_bb++];
					b->snid = i;
					b->vs = sub->v[jst].v;
					b->ve = t->v;
					b->ss = sst->soff + sst->len;
					b->se = sen->soff;
					b->len_min = ba[j].sd - ba[jst].sd - sst->len;
					b->len_max = ba[j].ld - ba[jst].ld - sst->len;
					b->n_paths = bb_n_paths(g, sub, jst, j);
					//fprintf(stderr, "X\t%s[%d]\tvs=%c%s\tve=%c%s\tlen_min=%d\n", g->sseq[i].name, i, "><"[b->vs&1], g->seg[b->vs>>1].name, "><"[b->ve&1], g->seg[b->ve>>1].name, b->len_min);
					assert(b->len_min >= 0);
					assert(b->len_max >= 0 && b->len_max >= b->len_min);
					b->n_seg = j - jst + 1;
					l = (b->len_min + 1) + (b->len_max + 1);
					l = (l + 3) / 4 + b->n_seg;
					GFA_CALLOC(b->v, l);
					b->seq_min = (char*)(b->v + b->n_seg);
					b->seq_max = b->seq_min + b->len_min + 1;
					for (k = jst; k <= j; ++k)
						b->v[k - jst] = sub->v[k].v;

					// test bubble involving both strands (mostly inversions)
					if (b->n_seg > m_vtmp) {
						m_vtmp = b->n_seg;
						kroundup32(m_vtmp);
						GFA_REALLOC(vtmp, m_vtmp);
					}
					for (k = 0; k < b->n_seg; ++k) vtmp[k] = b->v[k]>>1;
					radix_sort_gfa32(vtmp, vtmp + b->n_seg);
					for (k = 1; k < b->n_seg; ++k)
						if (vtmp[k] == vtmp[k-1]) break;
					b->is_bidir = (k < b->n_seg);

					// generate sequences and cf_min/cf_max
					GFA_MALLOC(v, j - jst);
					k = j, n = 0;
					while (k > jst) {
						if (k < j) v[n++] = sub->v[k].v;
						k = ba[k].sp;
					}
					bb_write_seq(g, n, v, b->len_min, b->seq_min);
					k = j, n = 0;
					while (k > jst) {
						if (k < j) v[n++] = sub->v[k].v;
						k = ba[k].lp;
					}
					bb_write_seq(g, n, v, b->len_max, b->seq_max);
					free(v);
				} // ~if(sst->snid==i&&sen->snid==i)
				max_a = max_soff = -1, jst = j;
			} // ~if(j==max_a)
			for (k = 0; k < t->n; ++k)
				if ((int32_t)(sub->a[t->off + k]>>32) > max_a)
					max_a = sub->a[t->off + k]>>32;
			if (g->seg[t->v>>1].snid == i && g->seg[t->v>>1].soff > max_soff)
				max_soff = g->seg[t->v>>1].soff;
		}
		free(ba);
		gfa_sub_destroy(sub);
	}
	free(vtmp);
	gfa_scbuf_destroy(scbuf);
	free(vs);
	*n_bb_ = n_bb;
	return bb;
}
