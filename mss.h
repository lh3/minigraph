#ifndef MSS_H
#define MSS_H

#include <stdint.h>
#define __STDC_LIMIT_MACROS
#include "kalloc.h"

#define MSS_TYPE    int32_t
#define MSS_NEG_INF INT32_MIN

typedef struct {
    int32_t st, en;
    MSS_TYPE sc;
} msseg_t;

msseg_t *mss_find_all(void *km, int32_t n, const MSS_TYPE *S, MSS_TYPE min_sc, MSS_TYPE xdrop, int32_t *n_seg);

#endif
