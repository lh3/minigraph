#ifndef __GFA_PRIV_H__
#define __GFA_PRIV_H__

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
	int32_t *a;
	void *km;
} gfa_sub_t;

#ifdef __cplusplus
extern "C" {
#endif

char *gfa_strdup(const char *src);
char *gfa_strndup(const char *src, size_t n);

// add/delete one segment/arc/stable sequence
int32_t gfa_add_seg(gfa_t *g, const char *name);
gfa_arc_t *gfa_add_arc1(gfa_t *g, uint32_t v, uint32_t w, int32_t ov, int32_t ow, int64_t link_id, int comp);
int32_t gfa_sseq_add(gfa_t *g, const char *sname);
void gfa_sseq_update(gfa_t *g, const gfa_seg_t *s);

// whole graph operations
void gfa_arc_sort(gfa_t *g);
void gfa_arc_index(gfa_t *g);
uint32_t gfa_fix_symm(gfa_t *g);
void gfa_symm(gfa_t *g); // delete multiple edges and restore skew-symmetry
void gfa_arc_rm(gfa_t *g);
void gfa_cleanup(gfa_t *g); // permanently delete arcs marked as deleted, sort and then index
void gfa_finalize(gfa_t *g);

// assembly related routines
int gfa_arc_del_trans(gfa_t *g, int fuzz); // transitive reduction
int gfa_arc_del_short(gfa_t *g, float drop_ratio); // delete short arcs
int gfa_cut_tip(gfa_t *g, int max_ext); // cut tips
int gfa_cut_internal(gfa_t *g, int max_ext); // drop internal segments
int gfa_cut_biloop(gfa_t *g, int max_ext); // Hmm... I forgot... Some type of weird local topology
int gfa_pop_bubble(gfa_t *g, int max_dist); // bubble popping
gfa_t *gfa_ug_gen(const gfa_t *g);

// subset, modifying the graph
void gfa_sub(gfa_t *g, int n, char *const* seg, int step);

// subset, without modifying the graph
gfa_sub_t *gfa_sub_from(void *km0, const gfa_t *g, uint32_t v0, int32_t max_dist);
void gfa_sub_destroy(gfa_sub_t *sub);
void gfa_sub_print(FILE *fp, const gfa_t *g, const gfa_sub_t *sub);

// graph augmentation
int gfa_ins_adj(const gfa_t *g, int min_len, gfa_ins_t *ins, const char *seq);
int32_t gfa_ins_filter(const gfa_t *g, int32_t n_ins, gfa_ins_t *ins);
void gfa_augment(gfa_t *g, int32_t n_ins, const gfa_ins_t *ins, int32_t n_ctg, const char *const* name, const char *const* seq);

gfa_sfa_t *gfa_gfa2sfa(const gfa_t *g, int32_t *n_sfa_, int32_t write_seq);

void gfa_blacklist_print(const gfa_t *g, FILE *fp, int32_t min_len); // FIXME: doesn't work with translocations

#ifdef __cplusplus
}
#endif

#endif // ~__GFA_PRIV_H__
