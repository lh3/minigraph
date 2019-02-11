#ifndef MINIGRAPH_H
#define MINIGRAPH_H

#include <stdint.h>
#include "gfa.h"

#define MG_VERSION "r1"

#define MG_I_HPC 0x1

typedef struct { uint64_t x, y; } mg128_t;
typedef struct { size_t n, m; mg128_t *a; } mg128_v;

typedef struct {
	uint32_t flag;
	int w, k;
	int bucket_bits;
} mg_idxopt_t;

typedef struct {
	uint64_t flag;
	int bw;
} mg_mapopt_t;

typedef struct {
	int32_t b, w, k, flag;
	gfa_t *g;
	struct mg_idx_bucket_s *B; // index (hidden)
} mg_idx_t;

typedef struct {
	int32_t id;
	int32_t as, cnt;
	int32_t rid;
	int32_t rs, re, qs, qe;
	int32_t sc_chain;
	int32_t mlen, blen;
	uint32_t hash;
	uint32_t rev:1, dummy:31;
} mg_lchain1_t;

extern int mg_verbose;
extern double mg_realtime0;

#ifdef __cplusplus
extern "C" {
#endif

int mg_opt_set(const char *preset, mg_idxopt_t *io, mg_mapopt_t *mo);
int mg_opt_check(const mg_idxopt_t *io, const mg_mapopt_t *mo);

mg_idx_t *mg_index_gfa(gfa_t *g, int k, int w, int b, int flag, int n_threads);
mg_idx_t *mg_index_file(const char *fn, int k, int w, int b, int flag, int n_threads);
void mg_idx_destroy(mg_idx_t *gi);

#ifdef __cplusplus
}
#endif

#endif
