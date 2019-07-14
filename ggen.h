#ifndef MG_GGEN_H
#define MG_GGEN_H

#include "minigraph.h"
#include "bseq.h"

#ifdef __cplusplus
extern "C" {
#endif

int32_t mg_fastcmp(void *km, int32_t l1, const char *s1, int32_t l2, const char *s2, int32_t k, int32_t max_occ);
int32_t mg_path2seq(void *km, const gfa_t *g, const mg_gchains_t *gcs, int32_t ls, int32_t le, int32_t voff[2], char **seq_, int32_t *cap_);
void mg_ggsimple(void *km, const mg_ggopt_t *opt, gfa_t *g, int32_t n_seq, const mg_bseq1_t *seq, mg_gchains_t *const* gcs);

#ifdef __cplusplus
}
#endif

#endif
