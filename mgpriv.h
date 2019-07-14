#ifndef MGPRIV_H
#define MGPRIV_H

#include <stdlib.h>
#include "kalloc.h"
#include "bseq.h"
#include "minigraph.h"

#define MG_DBG_NO_KALLOC   0x1
#define MG_DBG_QNAME       0x2
#define MG_DBG_SEED        0x4
#define MG_DBG_LCHAIN      0x8
#define MG_DBG_INSERT      0x10

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

typedef struct {
	uint32_t st, en:31, rev:1;
	int32_t far, i;
} mg_intv_t;

// shortest path
typedef struct {
	uint32_t v;
	int32_t target_dist;
	int32_t dist, n_path, path_end;
	int32_t meta;
	uint32_t target_hash, hash;
} mg_path_dst_t;

typedef struct {
	uint32_t v, d;
	int32_t pre;
} mg_pathv_t;

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
int32_t mg_idx_cal_max_occ(const mg_idx_t *gi, float f);

uint64_t *mg_chain_backtrack(void *km, int64_t n, const int32_t *f, const int32_t *p, int32_t *v, int32_t *t, int32_t min_cnt, int32_t min_sc, int32_t extra_u, int32_t *n_u_, int32_t *n_v_);
mg128_t *mg_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, int is_cdna, int n_segs, int64_t n, mg128_t *a, int *n_u_, uint64_t **_u, void *km);
mg_lchain_t *mg_lchain_gen(void *km, uint32_t hash, int qlen, int n_u, uint64_t *u, mg128_t *a);
void mg_update_anchors(int32_t n_a, mg128_t *a, int32_t n, const int32_t *mini_pos);

mg_pathv_t *mg_shortest_k(void *km0, const gfa_t *g, uint32_t src, int32_t n_dst, mg_path_dst_t *dst, int32_t max_dist, int32_t max_k, int32_t ql, const char *qs, int32_t *n_pathv);
int32_t mg_gchain1_dp(void *km, const gfa_t *g, int32_t *n_lc_, mg_lchain_t *lc, int32_t qlen, int32_t max_dist_g, int32_t max_dist_q, int32_t bw,
					  const char *qseq, uint64_t **u_);
mg_gchains_t *mg_gchain_gen(void *km_dst, void *km, const gfa_t *g, int32_t n_u, const uint64_t *u, const mg_lchain_t *lc, const mg128_t *a,
							uint32_t hash, int32_t min_gc_cnt, int32_t min_gc_score);
void mg_gchain_free(mg_gchains_t *gs);

void mg_gchain_restore_order(void *km, mg_gchains_t *gcs);
void mg_gchain_restore_offset(mg_gchains_t *gcs);
void mg_gchain_sort_by_score(void *km, mg_gchains_t *gcs);
void mg_gchain_set_parent(void *km, float mask_level, int n, mg_gchain_t *r, int sub_diff, int hard_mask_level);
int mg_gchain_flt_sub(float pri_ratio, int min_diff, int best_n, int n, mg_gchain_t *r);
void mg_gchain_drop_flt(void *km, mg_gchains_t *gcs);
void mg_gchain_set_mapq(void *km, mg_gchains_t *gcs, int qlen, int max_mini, int min_gc_score);

void mg_print_lchain(FILE *fp, const mg_idx_t *gi, int n_lc0, const mg_lchain_t *lc, const mg128_t *a, const char *qname);
void mg_write_paf(kstring_t *s, const gfa_t *g, const mg_gchains_t *gs, int32_t qlen, const char *qname, uint64_t flag, void *km);

int32_t mg_intv_index(int32_t n, mg_intv_t *a);
int32_t mg_intv_overlap(void *km, int32_t n_a, const mg_intv_t *a, int32_t st, int32_t en, int32_t **b_, int32_t *m_b_);

int32_t mg_fastcmp(void *km, int32_t l1, const char *s1, int32_t l2, const char *s2, int32_t k, int32_t max_occ);
int32_t mg_path2seq(void *km, const gfa_t *g, const mg_gchains_t *gcs, int32_t ls, int32_t le, int32_t voff[2], char **seq_, int32_t *cap_);
void mg_ggsimple(void *km, const mg_ggopt_t *opt, gfa_t *g, int32_t n_seq, const mg_bseq1_t *seq, mg_gchains_t *const* gcs);

void radix_sort_128x(mg128_t *beg, mg128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

#ifdef __cplusplus
}
#endif

#endif
