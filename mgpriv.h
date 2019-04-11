#ifndef MGPRIV_H
#define MGPRIV_H

#include <stdlib.h>
#include "kalloc.h"
#include "minigraph.h"

#define MG_DBG_NO_KALLOC   0x1
#define MG_DBG_PRINT_QNAME 0x2
#define MG_DBG_PRINT_SEED  0x4

#define MG_SEED_IGNORE     (1ULL<<41)
#define MG_SEED_TANDEM     (1ULL<<42)

#define MG_SEED_SEG_SHIFT  48
#define MG_SEED_SEG_MASK   (0xffULL<<(MG_SEED_SEG_SHIFT))

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	unsigned l, m;
	char *s;
} kstring_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned char seq_nt4_table[256];

void mg_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mg128_v *p);

const uint64_t *mg_idx_get(const mg_idx_t *gi, uint64_t minier, int *n);

mg128_t *mg_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int min_cnt, int min_sc, int is_cdna, int n_segs, int64_t n, mg128_t *a, int *n_u_, uint64_t **_u, void *km);
mg_lchain1_t *mg_lchain_gen(void *km, uint32_t hash, int qlen, int n_u, uint64_t *u, mg128_t *a);

void radix_sort_128x(mg128_t *beg, mg128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);

double realtime(void);
double cputime(void);
long peakrss(void);

#ifdef __cplusplus
}
#endif

#endif
