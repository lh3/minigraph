#include <assert.h>
#include "mgpriv.h"
#include "gfa-priv.h"
#include "algo.h"

int32_t mg_gc_index(void *km, int min_mapq, int min_map_len, int min_depth_len, const gfa_t *g, int32_t n_seq, mg_gchains_t *const* gcs,
					double *a_dens, int32_t **soff_, int32_t **qoff_, mg_intv_t **sintv_, mg_intv_t **qintv_);

typedef struct {
	int32_t bid;
	uint8_t is_stem;
} callaux_t;

void mg_call_asm(const gfa_t *g, int32_t n_seq, mg_gchains_t *const *gcs, int32_t min_mapq, int32_t min_blen)
{
	int32_t i, j, max_acnt, *soff, *qoff, n_bb;
	mg_intv_t *sintv, *qintv;
	double a_dens;
	gfa_bubble_t *bb;
	callaux_t *ca;

	max_acnt = mg_gc_index(0, min_mapq, min_blen>>1, min_blen, g, n_seq, gcs, &a_dens, &soff, &qoff, &sintv, &qintv);
	if (max_acnt == 0) return;

	bb = gfa_bubble(g, &n_bb);
	GFA_CALLOC(ca, gfa_n_vtx(g));
	for (i = 0; i < n_bb; ++i) {
		gfa_bubble_t *b = &bb[i];
		assert(b->n_seg >= 2);
		for (j = 0; j < b->n_seg; ++j)
			ca[b->v[j]].bid = ca[b->v[j]^1].bid = i;
		ca[b->v[0]].is_stem = ca[b->v[b->n_seg-1]].is_stem = 1;
	}

	for (i = 0; i < n_bb; ++i) {
		gfa_bubble_t *b = &bb[i];
		printf("%s\t%d\t%d\n", g->sseq[b->snid].name, b->ss, b->se);
	}

	free(ca);
	free(soff); free(qoff); free(sintv); free(qintv);
	for (i = 0; i < n_bb; ++i) free(bb[i].v);
	free(bb);
}
