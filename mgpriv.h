#ifndef MGPRIV_H
#define MGPRIV_H

#include <stdlib.h>
#include "kalloc.h"
#include "minigraph.h"

#define MG_SEED_SEG_SHIFT  48
#define MG_SEED_SEG_MASK   (0xffULL<<(MG_SEED_SEG_SHIFT))

#ifdef __cplusplus
extern "C" {
#endif

void mg_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mg128_v *p);

mg128_t *mg_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int min_cnt, int min_sc, int is_cdna, int n_segs, int64_t n, mg128_t *a, int *n_u_, uint64_t **_u, void *km);

void radix_sort_128x(mg128_t *beg, mg128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);

double realtime(void);
double cputime(void);
long peakrss(void);

#ifdef __cplusplus
}
#endif

#endif
