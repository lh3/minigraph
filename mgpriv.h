#ifndef MGPRIV_H
#define MGPRIV_H

#include <stdlib.h>
#include "kalloc.h"
#include "minigraph.h"

#ifdef __cplusplus
extern "C" {
#endif

void mg_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mg128_v *p);

void radix_sort_128x(mg128_t *beg, mg128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);

double realtime(void);
double cputime(void);
long peakrss(void);

#ifdef __cplusplus
}
#endif

#endif
