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

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int mg_ilog2_32(uint32_t v)
{
	uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

extern unsigned char seq_nt4_table[256];

void mg_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, mg128_v *p);

void *mg_idx_a2h(void *km, int32_t n_a, mg128_t *a, int suflen, uint64_t **q_, int32_t *n_);
const uint64_t *mg_idx_hget(const void *h_, const uint64_t *q, int suflen, uint64_t minier, int *n);

const uint64_t *mg_idx_get(const mg_idx_t *gi, uint64_t minier, int *n);

uint64_t *mg_chain_backtrack(void *km, int64_t n, const int32_t *f, const int32_t *p, int32_t *v, int32_t *t, int32_t min_cnt, int32_t min_sc, int32_t extra_u, int32_t *n_u_, int32_t *n_v_);
mg128_t *mg_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int min_cnt, int min_sc, int is_cdna, int n_segs, int64_t n, mg128_t *a, int *n_u_, uint64_t **_u, void *km);
mg_lchain_t *mg_lchain_gen(void *km, uint32_t hash, int qlen, int n_u, uint64_t *u, mg128_t *a);
void mg_update_anchors(int32_t n_a, mg128_t *a, int32_t n, const int32_t *mini_pos);

int32_t mg_gchain1_dp(void *km, const gfa_t *g, int32_t n_frag, mg_lchain_t *frag, int32_t qlen, int32_t max_dist_g, int32_t max_dist_q, int32_t bw, uint64_t **u_);
mg_gchains_t *mg_gchain_gen(void *km_dst, void *km, const gfa_t *g, int32_t n_u, const uint64_t *u, const mg_lchain_t *lc, const mg128_t *a, uint32_t hash, int32_t min_gc_cnt, int32_t min_gc_score);
void mg_gchain_free(mg_gchains_t *gs);

void mg_print_lchain(FILE *fp, const mg_idx_t *gi, int n_lc0, const mg_lchain_t *lc, const mg128_t *a, const char *qname);
void mg_write_paf(kstring_t *s, const gfa_t *g, const mg_gchains_t *gs, int32_t qlen, const char *qname, void *km);
void mg_print_paf(FILE *fp, const gfa_t *g, const mg_gchains_t *gs, int32_t qlen, const char *qname, void *km);

void radix_sort_128x(mg128_t *beg, mg128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);

void mg_err_fputs(const char *str, FILE *fp);

double realtime(void);
double cputime(void);
long peakrss(void);

#ifdef __cplusplus
}
#endif

#endif
