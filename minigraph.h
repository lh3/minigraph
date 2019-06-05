#ifndef MINIGRAPH_H
#define MINIGRAPH_H

#include <stdint.h>
#include "gfa.h"

#define MG_VERSION "r1"

#define MG_M_SPLICE       0x10
#define MG_M_SR           0x20
#define MG_M_FRAG_MODE    0x40
#define MG_M_FOR_ONLY     0x100
#define MG_M_REV_ONLY     0x200
#define MG_M_HEAP_SORT    0x400
#define MG_M_COPY_COMMENT 0x10000
#define MG_M_INDEPEND_SEG 0x20000
#define MG_M_NO_QUAL      0x40000
#define MG_M_2_IO_THREADS 0x80000

#define MG_MAX_SEG        255

typedef struct { uint64_t x, y; } mg128_t;
typedef struct { size_t n, m; mg128_t *a; } mg128_v;

typedef struct {
	uint32_t flag;
	int w, k;
	int bucket_bits;
} mg_idxopt_t;

typedef struct {
	uint64_t flag;
	int seed;
	int mini_batch_size;
	int max_qlen;
	int pe_ori;
	int mid_occ, max_occ;
	int bw, max_gap, max_gap_ref, max_frag_len;
	int max_chain_skip;
	int min_lc_cnt, min_lc_score;
	int min_gc_cnt, min_gc_score;
} mg_mapopt_t;

typedef struct {
	int32_t b, w, k, flag;
	gfa_t *g;
	struct mg_idx_bucket_s *B; // index (hidden)
} mg_idx_t;

typedef struct {
	int32_t as, cnt;
	uint32_t v;
	int32_t rs, re, qs, qe;
	int32_t score;
} mg_lchain_t;

typedef struct {
	int32_t as, cnt;
	uint32_t v;
	int32_t score;
} mg_llchain_t;

typedef struct {
	int32_t ls, cnt;
	int32_t n_anchor, score;
} mg_gchain_t;

typedef struct {
	void *km;
	int32_t n_g, n_l, n_a;
	mg_gchain_t *g;
	mg_llchain_t *l;
	mg128_t *a;
} mg_gchains_t;

typedef struct mg_tbuf_s mg_tbuf_t;

extern int mg_verbose, mg_dbg_flag;
extern double mg_realtime0;

#ifdef __cplusplus
extern "C" {
#endif

int mg_opt_set(const char *preset, mg_idxopt_t *io, mg_mapopt_t *mo);
int mg_opt_check(const mg_idxopt_t *io, const mg_mapopt_t *mo);

mg_idx_t *mg_index_gfa(gfa_t *g, int k, int w, int b, int flag, int n_threads);
mg_idx_t *mg_index_file(const char *fn, int k, int w, int b, int flag, int n_threads);
void mg_idx_destroy(mg_idx_t *gi);

int mg_map_file(const mg_idx_t *idx, const char *fn, const mg_mapopt_t *opt, int n_threads);

#ifdef __cplusplus
}
#endif

#endif
