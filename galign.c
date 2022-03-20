#include <assert.h>
#include <string.h>
#include "mgpriv.h"
#include "kalloc.h"

static void append_cigar1(void *km, mg32_v *c, int32_t op, int32_t len)
{
	if (c->n > 0 && (c->a[c->n - 1]&0xf) == op) {
		c->a[c->n - 1] += len<<4;
	} else {
		if (c->n == c->m) {
			c->m += (c->m>>1) + 16;
			KREALLOC(km, c->a, c->m);
		}
		c->a[c->n++] = len<<4 | op;
	}
}

void mg_gchain_cigar(void *km, const gfa_t *g, const gfa_edseq_t *es, const char *qseq, mg_gchains_t *gt)
{
	int32_t i, l_seq = 0, m_seq = 0;
	char *seq = 0;
	mg32_v cigar = {0,0,0};
	for (i = 0; i < gt->n_gc; ++i) {
		mg_gchain_t *gc = &gt->gc[i];
		int32_t l0 = gc->off;
		int32_t off_a0 = gt->lc[l0].off;
		int32_t j, j0 = 0, l;
		cigar.n = 0;
		//append_cigar1(km, &cigar, 7, gt->a[off_a0].y>>32&0xff);
		for (j = 1; j < gc->n_anchor; ++j) {
			const mg128_t *q, *p = &gt->a[off_a0 + j];
			if ((p->y & MG_SEED_IGNORE) && j != gc->n_anchor - 1) continue;
			q = &gt->a[off_a0 + j0];
			// find the lchain that contains the anchor
			for (l = l0; l < gc->off + gc->cnt; ++l) {
				mg_llchain_t *r = &gt->lc[l];
				if (off_a0 + j >= r->off && off_a0 + j < r->off + r->cnt)
					break;
			}
			assert(l < gc->off + gc->cnt);
			// do the alignment
			if (l == l0) { // on the same vertex
				l_seq = (int32_t)p->x - (int32_t)q->x;
				if (l_seq + 1 > m_seq) {
					m_seq = l_seq + 1;
					kroundup32(m_seq);
					KREALLOC(km, seq, m_seq);
				}
				memcpy(seq, &es[gt->lc[l].v].seq[(int32_t)q->x + 1], l_seq);
			} else {
			}
			j0 = j, l0 = l;
		}
	}
	kfree(km, seq);
}
