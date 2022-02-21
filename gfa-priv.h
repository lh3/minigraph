#ifndef __GFA_PRIV_H__
#define __GFA_PRIV_H__

#include <stdlib.h>
#include "gfa.h"

#define GFA_MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define GFA_CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define GFA_REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
#define GFA_BZERO(ptr, len) memset((ptr), 0, (len) * sizeof(*(ptr)))
#define GFA_EXPAND(a, m) do { \
		(m) = (m)? (m) + ((m)>>1) : 16; \
		GFA_REALLOC((a), (m)); \
	} while (0)

typedef struct { uint64_t x, y; } gfa128_t;

// linearized subgraphs

typedef struct {
	uint32_t v, d;
	int32_t off, n;
} gfa_subv_t;

typedef struct {
	int32_t n_v, n_a, is_dag;
	gfa_subv_t *v;
	uint64_t *a; // high 32 bits: point to the neighbor; low 32 bit: arc index in the graph
	void *km;
} gfa_sub_t;

typedef struct {
	int32_t snid, ss, se;
	uint32_t vs, ve;
	int32_t is_bidir, n_seg, len_max, len_min;
	uint32_t *v, n_paths;
	char *seq_max, *seq_min; // seq_max and seq_min point to v[]
} gfa_bubble_t;

struct gfa_scbuf_s;
typedef struct gfa_scbuf_s gfa_scbuf_t;

#ifdef __cplusplus
extern "C" {
#endif

char *gfa_strdup(const char *src);
char *gfa_strndup(const char *src, size_t n);
void radix_sort_gfa64(uint64_t *st, uint64_t *en);

// add/delete one segment/arc/stable sequence
int32_t gfa_add_seg(gfa_t *g, const char *name);
gfa_arc_t *gfa_add_arc1(gfa_t *g, uint32_t v, uint32_t w, int32_t ov, int32_t ow, int64_t link_id, int comp);
int32_t gfa_sseq_get(const gfa_t *g, const char *sname);
int32_t gfa_sseq_add(gfa_t *g, const char *sname);
void gfa_sseq_update(gfa_t *g, const gfa_seg_t *s);

// whole graph operations
void gfa_arc_sort(gfa_t *g);
void gfa_arc_index(gfa_t *g);
uint32_t gfa_fix_symm_add(gfa_t *g);
void gfa_fix_symm_del(gfa_t *g); // delete multiple edges and restore skew-symmetry
void gfa_arc_rm(gfa_t *g);
void gfa_cleanup(gfa_t *g); // permanently delete arcs marked as deleted, sort and then index
void gfa_finalize(gfa_t *g);
int32_t gfa_check_multi(const gfa_t *g);
uint32_t gfa_fix_multi(gfa_t *g);

int gfa_arc_del_multi_risky(gfa_t *g);
int gfa_arc_del_asymm_risky(gfa_t *g);

// assembly related routines
int gfa_arc_del_trans(gfa_t *g, int fuzz); // transitive reduction
int gfa_arc_del_weak(gfa_t *g);
int gfa_arc_pair_strong(gfa_t *g);
int gfa_arc_del_short(gfa_t *g, int min_ovlp_len, float drop_ratio); // delete short arcs
int gfa_drop_tip(gfa_t *g, int tip_cnt, int tip_len); // cut tips
int gfa_drop_internal(gfa_t *g, int max_ext);
int gfa_cut_z(gfa_t *g, int32_t min_dist, int32_t max_dist);
int gfa_topocut(gfa_t *g, float drop_ratio, int32_t tip_cnt, int32_t tip_len);
int gfa_bub_simple(gfa_t *g, int min_side, int max_side);
int gfa_pop_bubble(gfa_t *g, int radius, int max_del, int protect_tip); // bubble popping
gfa_t *gfa_ug_gen(const gfa_t *g);
void gfa_scc_all(const gfa_t *g);

// subset, modifying the graph
void gfa_sub(gfa_t *g, int n, char *const* seg, int step);
char **gfa_query_by_reg(const gfa_t *g, int32_t n_bb, const gfa_bubble_t *bb, const char *reg, int *n_seg);

// subset, without modifying the graph
gfa_sub_t *gfa_sub_from(void *km0, const gfa_t *g, uint32_t v0, int32_t max_dist);
void gfa_sub_destroy(gfa_sub_t *sub);
void gfa_sub_print(FILE *fp, const gfa_t *g, const gfa_sub_t *sub);

gfa_scbuf_t *gfa_scbuf_init(const gfa_t *g);
gfa_sub_t *gfa_scc1(void *km0, const gfa_t *g, gfa_scbuf_t *b, uint32_t v0);
void gfa_scbuf_destroy(gfa_scbuf_t *b);

// graph augmentation
int gfa_ins_adj(const gfa_t *g, int min_len, gfa_ins_t *ins, const char *seq);
int32_t gfa_ins_filter(const gfa_t *g, int32_t n_ins, gfa_ins_t *ins);
void gfa_augment(gfa_t *g, int32_t n_ins, const gfa_ins_t *ins, int32_t n_ctg, const char *const* name, const char *const* seq);

gfa_sfa_t *gfa_gfa2sfa(const gfa_t *g, int32_t *n_sfa_, int32_t write_seq);

void gfa_sort_ref_arc(gfa_t *g);
gfa_bubble_t *gfa_bubble(const gfa_t *g, int32_t *n_); // FIXME: doesn't work with translocation

void gfa_gt_simple_print(const gfa_t *g, float min_dc, int32_t is_path); // FIXME: doesn't work with translocations

void gfa_aux_update_cv(gfa_t *g, const char *tag, const double *cov_seg, const double *cov_link);

void gfa_sql_write(FILE *fp, const gfa_t *g, int write_seq);

static inline int64_t gfa_find_arc(const gfa_t *g, uint32_t v, uint32_t w)
{
	uint32_t i, nv = gfa_arc_n(g, v), nw = 0, k = (uint32_t)-1;
	gfa_arc_t *av = gfa_arc_a(g, v);
	for (i = 0; i < nv; ++i)
		if (av[i].w == w) ++nw, k = i;
	return nw == 1? (int64_t)(&av[k] - g->arc) : nw == 0? -1 : -2;
}

#ifdef __cplusplus
}
#endif

#endif // ~__GFA_PRIV_H__
