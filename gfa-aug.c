#include <assert.h>
#include "gfa-priv.h"

typedef struct {
	uint32_t side;
	uint32_t ins:31, end:1;
} gfa_split_t;

#define split_key(p) ((p).side)
KRADIX_SORT_INIT(split, gfa_split_t, split_key, 4)

void gfa_augment(gfa_t *g, int32_t n_ins, const gfa_ins_t *ins, int32_t n_ctg, const char **name, const char **seq)
{
	int32_t i, j, k, *scnt, *soff, *pos, n_ctg_seg, n_old_seg, n_seg, *ins2seg;
	uint32_t *o2n;
	gfa_split_t *sp;
	gfa_seg_t *seg;

	if (n_ins <= 0 || n_ctg <= 0) return 0;

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
			q->side = (p->v[k]&1? vlen - q->voff[k] : q->voff[k]) << 1 | ((p->v[k]&1) ^ k);
			assert((q->side != 0<<1|0) && (q->side != vlen<<1|1)); // not possible to link such sides
			++scnt[q->v[k]>>1];
		}
		if (p->coff[1] > p->coff[0])
			++n_ctg_seg;
	}

	// sort
	for (j = 0, n_old_seg = 0; j < g->n_seg; ++j)
		if (soff[j+1] - soff[j] > 1)
			radix_sort(split, &sp[soff[j]], &sp[soff[j+1]]);

	// precompute the number of segments after split
	GFA_BZERO(scnt, g->n_seg);
	for (j = 0, n_old_seg = 0; j < g->n_seg; ++j) {
		int32_t i0;
		for (i0 = soff[j], i = i0 + 1, k = 0; i <= soff[j+1]; ++i)
			if (i == soff[j+1] || sp[i0].side>>1 != sp[i].side>>1) {
				if (sp[i0].side>>1 != 0 && sp[i0].side>>1 != g->seg[j].len) // otherwise no new segment will be created
					++k;
				i0 = i;
			}
		scnt[j] = k;
		n_old_seg += k + 1;
	}

	// create newly inserted segments
	n_seg = n_old_seg + n_ctg_seg;
	GFA_CALLOC(seg, n_seg);
	GFA_MALLOC(ins2seg, n_ins);
	for (i = 0, k = n_old_seg; i < n_ins; ++i) {
		const gfa_ins_t *p = &ins[i];
		char buf[16];
		gfa_seg_t *s;
		ins2seg[i] = -1;
		if (p->coff[1] <= p->coff[0]) continue; // no new segment created
		ins2seq[i] = k;
		s = &seg[k++];
		snprintf(buf, 15, "v%d", k);
		s->name = strdup(buf);
		GFA_MALLOC(s->seq, p->coff[1] - p->coff[0] + 1);
		for (j = 0; j < p->coff[1] - p->coff[0]; ++j)
			s->seq[j] = seq[i][p->coff[0] + j];
		s->seq[j] = 0;
		s->len = j;
		s->pnid = gfa_add_pname(g, name[i]);
		s->ppos = p->coff[0];
		s->rank = g->max_rank + 1; // TODO: to deal with SN/SS/SR tags somewhere
	}

	// populate seg[]
	GFA_MALLOC(o2n, g->n_seg * 2);
	for (j = 0, k = 0; j < g->n_seg; ++j) {
		int32_t i0, l, off = 0, k0 = k;
		for (i0 = soff[j], i = i0 + 1; i <= soff[j+1]; ++i) {
			if (i == soff[j+1] || sp[i].side>>1 != sp[i0].side>>1) {
				gfa_seg_t *s = &seg[k];
				for (l = i0; l < i; ++l) {
				}
				i0 = i;
			}
		}
		o2n[(uint32_t)j<<1] = (uint32_t)n_old_seg<<1;
		o2n[(uint32_t)j<<1|1] = (uint32_t)(n_old_seg-1)<<1|1;
	}
	assert(k == n_old_seg);

	free(sp); free(soff); free(scnt);
}
