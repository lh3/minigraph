#include <assert.h>
#include "gfa.h"

typedef struct {
	uint32_t pos;
	uint32_t ins:31, end:1;
} gfa_split_t;

#define split_key(p) ((p).pos)
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
			gfa_split_t *q = &sp[soff[p->v[k]>>1] + scnt[p->v[k]>>1]];
			q->ins = i, q->end = k;
			q->pos = q->vpos[k]<<1 | ((p->v[k]&1) ^ k);
			++scnt[q->v[k]>>1];
		}
		if (p->cpos[1] > p->cpos[0])
			++n_ctg_seg;
	}

	// count the number of distinct split and save to scnt[]
	for (j = 0, n_old_seg = 0; j < g->n_seg; ++j) {
		if (soff[j+1] - soff[j] > 1) {
			radix_sort(split, &sp[soff[j]], &sp[soff[j+1]]);
			for (i = soff[j], k = 0; i < soff[j+1]; ++i)
				if (i + 1 == soff[j+1] || sp[i].pos>>1 != sp[i+1].pos>>1)
					++k;
			scnt[j] = k + 1;
		} else scnt[j] = soff[j+1] - soff[j];
		n_old_seg += scnt[j] + 1;
	}

	// populate seg[]
	n_seg = n_old_seg + n_ctg_seg;
	GFA_MALLOC(o2n, g->n_seg * 2);
	GFA_CALLOC(seg, n_seg);
	GFA_MALLOC(ins2seg, n_ins);
	for (i = 0; i < n_ins; ++i) ins2seg[i] = -1;
	for (j = 0, n_old_seg = n_ctg_seg = 0; j < g->n_seg; ++j) {
		int32_t i0;
		o2n[(uint32_t)j<<1] = (uint32_t)n_old_seg<<1;
		for (i0 = soff[j], i = i0 + 1; i <= soff[j+1]; ++i) {
			if (i == soff[j+1] || sp[i].pos>>1 != sp[i0].pos>>1) {
				gfa_seg_t *s = &seg[n_old_seg];
				for (k = i0; k < i; ++k) {
				}
				++n_old_seg;
				i0 = i;
			}
		}
		o2n[(uint32_t)j<<1|1] = (uint32_t)(n_old_seg-1)<<1|1;
	}

	free(sp); free(soff); free(scnt);
}
