#ifndef MINIGRAPH_H
#define MINIGRAPH_H

#include "gfa.h"

#define MG_VERSION "r1"

#define MG_I_HPC 0x1

typedef struct { uint64_t x, y; } mg128_t;
typedef struct { size_t n, m; mg128_t *a; } mg128_v;

typedef struct {
	uint64_t flag;
	int w, k;
	int bucket_bits;
} mg_idxopt_t;

typedef struct {
} mg_mapopt_t;

typedef struct {
	int32_t b, w, k, flag;
	gfa_t *g;
	struct mg_idx_bucket_s *B; // index (hidden)
} mg_idx_t;

#endif
