#ifndef MG_ALGO_H
#define MG_ALGO_H

#include <stdint.h>

#define MSS_TYPE int32_t
#define LIS_TYPE uint64_t

typedef struct {
    int32_t st, en;
    MSS_TYPE sc;
} msseg_t;

#ifdef __cplusplus
extern "C" {
#endif

int32_t ks_lis_64(void *km, int32_t n, const LIS_TYPE *a, int32_t *b);
msseg_t *mss_find_all(void *km, int32_t n, const MSS_TYPE *S, MSS_TYPE min_sc, MSS_TYPE xdrop, int32_t *n_seg);

#ifdef __cplusplus
}
#endif

#endif
