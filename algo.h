#ifndef MG_ALGO_H
#define MG_ALGO_H

#include <stdint.h>

#define MG_MSS_TYPE int32_t
#define MG_LIS_TYPE uint64_t

typedef struct {
    int32_t st, en;
    MG_MSS_TYPE sc;
} mg_msseg_t;

typedef struct {
	uint32_t st, en:31, rev:1;
	int32_t far, i;
} mg_intv_t;

#ifdef __cplusplus
extern "C" {
#endif

mg_msseg_t *mg_mss_all(void *km, int32_t n, const MG_MSS_TYPE *S, MG_MSS_TYPE min_sc, MG_MSS_TYPE xdrop, int32_t *n_seg);
int32_t mg_intv_index(int32_t n, mg_intv_t *a);
int32_t mg_intv_overlap(void *km, int32_t n_a, const mg_intv_t *a, int32_t st, int32_t en, int32_t **b_, int32_t *m_b_);
void radix_sort_mg_intv(mg_intv_t *st, mg_intv_t *en);
int32_t mg_wfa_cmp(void *km, int32_t l1, const char *s1, int32_t l2, const char *s2, int32_t max_pen, int32_t *mlen, int32_t *blen);

#ifdef __cplusplus
}
#endif

#endif
