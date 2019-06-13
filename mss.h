#ifndef MSS_H
#define MSS_H

#include <stdint.h>
#define __STDC_LIMIT_MACROS
#include "kalloc.h"

#define MSS_TYPE    int64_t
#define MSS_NEG_INF INT64_MIN

typedef struct {
    int st, en;
    MSS_TYPE sc;
} msseg_t;

msseg_t *mss_find_all(void *km, int n, const MSS_TYPE *S, MSS_TYPE min_sc, MSS_TYPE xdrop, int *n_seg);

#endif
