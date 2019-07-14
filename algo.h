#ifndef MG_ALGO_H
#define MG_ALGO_H

#include <stdint.h>

#define MSS_TYPE int32_t
#define LIS_TYPE uint64_t

typedef struct {
    int32_t st, en;
    MSS_TYPE sc;
} mg_msseg_t;

typedef struct {
	uint32_t st, en:31, rev:1;
	int32_t far, i;
} mg_intv_t;

#ifdef __cplusplus
extern "C" {
#endif

int32_t mg_lis_64(void *km, int32_t n, const LIS_TYPE *a, int32_t *b);
mg_msseg_t *mg_mss_all(void *km, int32_t n, const MSS_TYPE *S, MSS_TYPE min_sc, MSS_TYPE xdrop, int32_t *n_seg);
int32_t mg_intv_index(int32_t n, mg_intv_t *a);
int32_t mg_intv_overlap(void *km, int32_t n_a, const mg_intv_t *a, int32_t st, int32_t en, int32_t **b_, int32_t *m_b_);

#ifdef __cplusplus
}
#endif

#endif
