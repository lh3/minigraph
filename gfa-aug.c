#include <assert.h>
#include "gfa-priv.h"
#include "ksort.h"

typedef struct {
	uint32_t side;
	uint32_t ins:31, end:1;
} gfa_split_t;

#define split_key(p) ((p).side)
KRADIX_SORT_INIT(split, gfa_split_t, split_key, 4)

static inline void create_first_arc(gfa_t *g, const gfa_seg_t *seg, uint32_t v, uint32_t w)
{
	gfa_arc_t *a;
	if (g->n_arc == g->m_arc) GFA_EXPAND(g->arc, g->m_arc);
	a = &g->arc[g->n_arc++];
	a->v_lv = (uint64_t)v<<32 | seg[v>>1].len;
	a->w = w;
	a->lw = seg[w>>1].len;
	a->ov = a->ow = 0;
	a->link_id = g->n_arc - 1;
	a->del = a->comp = 0;
}

void gfa_augment(gfa_t *g, int32_t n_ins, const gfa_ins_t *ins, int32_t n_ctg, const char *const* name, const char *const* seq)
{
	int32_t i, j, k, *scnt, *soff, n_ctg_seg, n_old_seg, n_seg;
	gfa_split_t *sp;
	gfa_seg_t *seg;
	char buf[16];
	uint32_t *old2new;
	uint64_t t, n_old_arc = g->n_arc, *ins_side;

	if (n_ins <= 0 || n_ctg <= 0) return;

	// set soff[]
	GFA_CALLOC(scnt, g->n_seg);
	for (i = 0; i < n_ins; ++i)
		++scnt[ins[i].v[0]>>1], ++scnt[ins[i].v[1]>>1];
	GFA_MALLOC(soff, g->n_seg + 1);
	for (j = 1, soff[0] = 0; j <= g->n_seg; ++j)
		soff[j] = soff[j-1] + scnt[j-1];

	// populate sp[]
	GFA_MALLOC(sp, soff[g->n_seg]);
	GFA_BZERO(scnt, g->n_seg);
	for (i = 0, n_ctg_seg = 0; i < n_ins; ++i) {
		const gfa_ins_t *p = &ins[i];
		for (k = 0; k < 2; ++k) {
			uint32_t vlen = g->seg[p->v[k]>>1].len;
			gfa_split_t *q = &sp[soff[p->v[k]>>1] + scnt[p->v[k]>>1]];
			q->ins = i, q->end = k;
			q->side = (p->v[k]&1? vlen - p->voff[k] : p->voff[k]) << 1 | ((p->v[k]&1) ^ k);
			assert(q->side != (0<<1|0) && q->side != (vlen<<1|1)); // not possible to link such sides
			++scnt[p->v[k]>>1];
		}
		if (p->coff[1] > p->coff[0])
			++n_ctg_seg;
	}
	free(scnt);

	// sort sp[]
	for (j = 0, n_old_seg = 0; j < g->n_seg; ++j)
		if (soff[j+1] - soff[j] > 1)
			radix_sort_split(&sp[soff[j]], &sp[soff[j+1]]);

	// precompute the number of segments after split
	for (j = 0, n_old_seg = 0; j < g->n_seg; ++j) {
		int32_t i0;
		for (i0 = soff[j], i = i0 + 1, k = 0; i <= soff[j+1]; ++i)
			if (i == soff[j+1] || sp[i0].side>>1 != sp[i].side>>1) {
				if (sp[i0].side>>1 != 0 && sp[i0].side>>1 != g->seg[j].len) // otherwise no new segment will be created
					++k;
				i0 = i;
			}
		n_old_seg += k + 1;
	}

	// compute ins_side[] and split old segments
	n_seg = n_old_seg + n_ctg_seg;
	GFA_CALLOC(seg, n_seg);
	g->is_srt = g->is_symm = 0;
	GFA_CALLOC(ins_side, n_ins);
	GFA_MALLOC(old2new, g->n_seg * 2);
	for (j = 0, k = 0; j < g->n_seg; ++j) {
		int32_t i0, l, off = 0, k0 = k;
		gfa_seg_t *s = &g->seg[j];
		gfa_seg_t *t = &seg[k]; // this is so far a placeholder
		// create the first half of a new segment
		snprintf(buf, 15, "v%d", k + 1);
		t->name = strdup(buf);
		t->pnid = s->pnid, t->ppos = s->ppos, t->rank = s->rank;
		// iterate over splits
		for (i0 = soff[j], i = i0 + 1; i <= soff[j+1]; ++i) {
			if (i == soff[j+1] || sp[i].side>>1 != sp[i0].side>>1) {
				gfa_split_t *q0 = &sp[i0];
				for (l = i0; l < i; ++l) {
					gfa_split_t *q = &sp[l];
					int32_t shift = q->end == 0? 32 : 0; // first end on the higher 32 bits
					int32_t side = q->side & 1;
					ins_side[q->ins] |= (uint64_t)((uint32_t)(k + side) << 1 | (side^q->end)) << shift;
				}
				if (q0->side>>1 != 0 && q0->side>>1 != g->seg[j].len) { // create a new segment
					t->len = (q0->side>>1) - off;
					GFA_MALLOC(t->seq, t->len + 1);
					memcpy(t->seq, &s->seq[off], t->len);
					t->seq[t->len] = 0;
					off += t->len;
					t = &seg[++k]; // create a new segment
					snprintf(buf, 15, "v%d", k + 1);
					t->name = strdup(buf);
					t->pnid = s->pnid, t->ppos = s->ppos + off, t->rank = s->rank;
				}
				i0 = i;
			}
		}
		// finish the last segment
		t->len = s->len - off;
		GFA_MALLOC(t->seq, t->len + 1);
		memcpy(t->seq, &s->seq[off], t->len);
		t->seq[t->len] = 0;
		++k;
		// translate table for old and new vertices
		old2new[j<<1|0] = (uint32_t)(k - 1) << 1 | 0;
		old2new[j<<1|1] = (uint32_t)k0 << 1 | 1;
		// add new arcs between newly created segments
		for (i = 0; i < k - k0 - 1; ++i)
			create_first_arc(g, seg, (uint32_t)(k0+i)<<1, (uint32_t)(k0+i+1)<<1);
	}
	assert(k == n_old_seg);
	free(soff);
	free(sp);

	// update existing g->arc[]
	for (t = 0; t < n_old_arc; ++t) {
		gfa_arc_t *a = &g->arc[t];
		uint32_t v = old2new[a->v_lv>>32];
		a->v_lv = (uint64_t)v << 32 | seg[v>>1].len;
		a->w = old2new[a->w];
		a->lw = seg[a->w>>1].len;
	}
	free(old2new);

	// create newly inserted segments
	for (i = 0, k = n_old_seg; i < n_ins; ++i) {
		const gfa_ins_t *p = &ins[i];
		if (p->coff[0] < p->coff[1]) {
			gfa_seg_t *t = &seg[k];
			snprintf(buf, 15, "v%d", k + 1);
			t->name = strdup(buf);
			GFA_MALLOC(t->seq, p->coff[1] - p->coff[0] + 1);
			for (j = 0; j < p->coff[1] - p->coff[0]; ++j)
				t->seq[j] = seq[p->ctg][p->coff[0] + j];
			t->seq[j] = 0;
			t->len = j;
			t->pnid = gfa_add_pname(g, name[p->ctg]);
			t->ppos = p->coff[0];
			t->rank = g->max_rank + 1; // TODO: to deal with SN/SS/SR tags somewhere
			create_first_arc(g, seg, ins_side[i]>>32, (uint32_t)k<<1);
			create_first_arc(g, seg, (uint32_t)k<<1, (uint32_t)ins_side[i]);
			++k;
		} else {
			create_first_arc(g, seg, ins_side[i]>>32, (uint32_t)ins_side[i]);
		}
	}
	free(ins_side);

	// update *g
	for (j = 0; j < g->n_seg; ++j) {
		free(g->seg[j].name);
		free(g->seg[j].seq);
		free(g->seg[j].aux.aux);
	}
	free(g->seg);
	g->seg = seg, g->n_seg = g->m_seg = n_seg;
	GFA_REALLOC(g->arc_aux, g->m_arc);
	GFA_BZERO(&g->arc_aux[n_old_arc], g->m_arc - n_old_arc);
	gfa_arc_sort(g);
	gfa_arc_index(g);
	gfa_fix_symm(g);
}
