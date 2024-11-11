#ifndef MGPRIV_H
#define MGPRIV_H

#include <stdlib.h>
#include "minigraph.h"

#define MG_DBG_NO_KALLOC   0x1
#define MG_DBG_QNAME       0x2
#define MG_DBG_SEED        0x4
#define MG_DBG_LCHAIN      0x8
#define MG_DBG_INSERT      0x10
#define MG_DBG_SHORTK      0x20
#define MG_DBG_GC1         0x40
#define MG_DBG_LC_PROF     0x80
#define MG_DBG_MINIWFA     0x100
#define MG_DBG_MWF_SEQ     0x200

#define MG_SEED_IGNORE     (1ULL<<41)
#define MG_SEED_TANDEM     (1ULL<<42)
#define MG_SEED_FIXED      (1ULL<<43)

#define MG_MAX_SEG        255
#define MG_SEED_SEG_SHIFT  48
#define MG_SEED_SEG_MASK   (0xffULL<<(MG_SEED_SEG_SHIFT))
#define mg_seg_id(a) ((int32_t)(((a).y&MG_SEED_SEG_MASK) >> MG_SEED_SEG_SHIFT))

#define MG_SEED_OCC_SHIFT  56

#define MG_MAX_SHORT_K  15

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	unsigned l, m;
	char *s;
} kstring_t;
#endif

// shortest path
typedef struct {
	// input
	uint32_t v;
	int32_t target_dist;
	uint32_t target_hash;
	uint32_t meta:30, check_hash:1, inner:1;
	int32_t qlen;
	// output
	uint32_t n_path:31, is_0:1;
	int32_t path_end;
	int32_t dist;
	uint32_t hash;
} mg_path_dst_t;

typedef struct {
	uint32_t v, d;
	int32_t pre;
} mg_pathv_t;

#ifdef __cplusplus
extern "C" {
#endif

static inline float mg_log2(float x) // NB: this doesn't work when x<2
{
	union { float f; uint32_t i; } z = { x };
	float log_2 = ((z.i >> 23) & 255) - 128;
	z.i &= ~(255 << 23);
	z.i += 127 << 23;
	log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
	return log_2;
}

extern unsigned char seq_nt4_table[256];

void mg_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, mg128_v *p);

void *mg_idx_a2h(void *km, int32_t n_a, mg128_t *a, int suflen, uint64_t **q_, int32_t *n_);
const uint64_t *mg_idx_hget(const void *h_, const uint64_t *q, int suflen, uint64_t minier, int *n);
void mg_idx_hfree(void *h_);

const uint64_t *mg_idx_get(const mg_idx_t *gi, uint64_t minier, int *n);
void mg_idx_cal_quantile(const mg_idx_t *gi, int32_t m, float f[], int32_t q[]);

uint64_t *mg_chain_backtrack(void *km, int64_t n, const int32_t *f, const int64_t *p, int32_t *v, int32_t *t, int32_t min_cnt, int32_t min_sc, int32_t max_drop,
							 int32_t extra_u, int32_t *n_u_, int32_t *n_v_);
mg128_t *mg_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					  int is_cdna, int n_seg, int64_t n, mg128_t *a, int *n_u_, uint64_t **_u, void *km);
mg128_t *mg_lchain_rmq(int max_dist, int max_dist_inner, int bw, int max_chn_skip, int cap_rmq_size, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					   int64_t n, mg128_t *a, int *n_u_, uint64_t **_u, void *km);
mg_lchain_t *mg_lchain_gen(void *km, uint32_t hash, int qlen, int n_u, uint64_t *u, mg128_t *a);
void mg_update_anchors(int32_t n_a, mg128_t *a, int32_t n, const int32_t *mini_pos);

mg_pathv_t *mg_shortest_k(void *km0, const gfa_t *g, uint32_t src, int32_t n_dst, mg_path_dst_t *dst, int32_t max_dist, int32_t max_k, int32_t *n_pathv);
int32_t mg_gchain1_dp(void *km, const gfa_t *g, int32_t *n_lc_, mg_lchain_t *lc, int32_t qlen, int32_t max_dist_g, int32_t max_dist_q, int32_t bw, int32_t max_skip,
					  int32_t ref_bonus, float chn_pen_gap, float chn_pen_skip, float mask_level, const mg128_t *an, uint64_t **u_);
mg_gchains_t *mg_gchain_gen(void *km_dst, void *km, const gfa_t *g, const gfa_edseq_t *es, int32_t n_u, const uint64_t *u,
							mg_lchain_t *lc, const mg128_t *a, uint32_t hash, int32_t min_gc_cnt, int32_t min_gc_score,
							int32_t gdp_max_ed, int32_t n_seg, const char *qseq);
void mg_gchain_cigar(void *km, const gfa_t *g, const gfa_edseq_t *es, const char *qseq, mg_gchains_t *gt, const char *qname);
void mg_gchain_gen_ds(void *km, const gfa_t *g, const gfa_edseq_t *es, const char *qseq, mg_gchains_t *gt);
void mg_gchain_free(mg_gchains_t *gs);

uint32_t *lv_ed_unified(void *km, int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t is_ext, int32_t *score, int32_t *t_endl, int32_t *q_endl, int32_t *n_cigar);

void mg_gchain_restore_order(void *km, mg_gchains_t *gcs);
void mg_gchain_restore_offset(mg_gchains_t *gcs);
void mg_gchain_sort_by_score(void *km, mg_gchains_t *gcs);
void mg_gchain_set_parent(void *km, float mask_level, int n, mg_gchain_t *r, int sub_diff, int hard_mask_level);
int mg_gchain_flt_sub(float pri_ratio, int min_diff, int best_n, int n, mg_gchain_t *r);
void mg_gchain_drop_flt(void *km, mg_gchains_t *gcs);
void mg_gchain_set_mapq(void *km, mg_gchains_t *gcs, int qlen, int max_mini, int min_gc_score);

void mg_cov_map(const gfa_t *g, const mg_gchains_t *gt, int32_t min_mapq, int32_t min_blen, double *c_seg, double *c_link, const char *qname);
void mg_cov_asm(const gfa_t *g, int32_t n_seq, mg_gchains_t *const *gcs, int32_t min_mapq, int32_t min_blen, double *cov_seg, double *cov_link);

void mg_print_lchain(FILE *fp, const mg_idx_t *gi, int n_lc0, const mg_lchain_t *lc, const mg128_t *a, const char *qname);
void mg_write_gaf(kstring_t *s, const gfa_t *g, const mg_gchains_t *gs, int32_t n_seg, const int32_t *qlens, const char *qname, uint64_t flag, void *km);

void mg_sprintf_lite(kstring_t *s, const char *fmt, ...);
void mg_sprintf_km(void *km, kstring_t *s, const char *fmt, ...);
void mg_str_write(void *km, kstring_t *s, int32_t len, char *str);
void mg_str_reserve(void *km, kstring_t *s, int32_t len);

void radix_sort_128x(mg128_t *beg, mg128_t *end);
void radix_sort_gfa64(uint64_t *beg, uint64_t *end);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

#ifdef __cplusplus
}
#endif

#endif
