#include "mgpriv.h"
#include "bseq.h"

gfa_t *mg_ggsimple(void *km, const mg_ggopt_t *opt, const gfa_t *g, int32_t n_seq, const mg_bseq1_t *seq, const mg_gchain_t **gcs)
{
	gfa_t *h = 0;
	int32_t *vcnt;
	vcnt = KCALLOC(km, int32_t, g->n_seg);
	kfree(km, vcnt);
	return h;
}
