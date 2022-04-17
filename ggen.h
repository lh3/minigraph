#ifndef MG_GGEN_H
#define MG_GGEN_H

#include "minigraph.h"
#include "bseq.h"

#ifdef __cplusplus
extern "C" {
#endif

int32_t mg_path2seq(void *km, const gfa_t *g, const mg_gchains_t *gcs, int32_t ls, int32_t le, int32_t voff[2], char **seq_, int32_t *cap_);
void mg_ggsimple(void *km, const mg_ggopt_t *opt, gfa_t *g, int32_t n_seq, const mg_bseq1_t *seq, mg_gchains_t *const* gcs);
void mg_ggsimple_cigar(void *km, const mg_ggopt_t *opt, gfa_t *g, int32_t n_seq, const mg_bseq1_t *seq, mg_gchains_t *const* gcs);

void mg_call_asm(const gfa_t *g, int32_t n_seq, const mg_bseq1_t *seq, mg_gchains_t *const *gcs, int32_t min_mapq, int32_t min_blen);

#ifdef __cplusplus
}
#endif

#endif
