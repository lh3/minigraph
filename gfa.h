#ifndef __GFA_H__
#define __GFA_H__

#include <stdio.h>
#include <stdint.h>

#define GFA_VERSION "r2"

#define GFA_MAX_SHORT_K  15

/*
  A segment is a sequence. A vertex is one side of a segment. In the code,
  segment_id is an integer, and vertex_id=segment_id<<1|orientation. The
  convention is to use variable u, v or w for a vertex, not for a segment. An
  arc is a directed edge between two vertices in the graph. Each arc has a
  complement arc. A link represents an arc and its complement. The following
  diagram shows an arc v->w, and the lengths used in the gfa_arc_t struct:

       |<--- lv --->|<-- ov -->|
    v: ------------------------>
                    ||overlap|||
                 w: -------------------------->
                    |<-- ow -->|<---- lw ---->|

  The graph topology is solely represented by an array of gfa_arc_t objects
  (see gfa_t::arc[]), where both an arc and its complement are present. The
  array is sorted by gfa_arc_t::v_lv and indexed by gfa_t::idx[] most of time.
  gfa_arc_a(g, v), of size gfa_arc_n(g, v), gives the array of arcs that leaves
  a vertex v in the graph g.
*/

typedef struct {
	uint64_t v_lv; // higher 32 bits: vertex_id; lower 32 bits: lv; packed together for sorting
	uint32_t w, lw;
	int32_t ov, ow;
	uint64_t link_id:62, del:1, comp:1;
} gfa_arc_t;

#define gfa_arc_head(a) ((uint32_t)((a).v_lv>>32))
#define gfa_arc_tail(a) ((a).w)
#define gfa_arc_len(a) ((uint32_t)(a).v_lv) // different from the original string graph

#define gfa_arc_n(g, v) ((uint32_t)(g)->idx[(v)])
#define gfa_arc_a(g, v) (&(g)->arc[(g)->idx[(v)]>>32])

typedef struct {
	uint32_t m_aux, l_aux;
	uint8_t *aux;
} gfa_aux_t;

typedef struct {
	uint32_t start, end; // start: starting vertex in the string graph; end: ending vertex
	uint32_t len2, dummy; // len_r: the other length of the unitig
	uint32_t m, n; // number of reads
	uint64_t *a; // list of reads
	char **name;
} gfa_utg_t;

typedef struct {
	int32_t len;
	uint32_t del:16, circ:16;
	int32_t pnid; // persistent name ID
	int32_t ppos; // persistent start position
	int32_t rank; // persistent rank
	char *name, *seq;
	gfa_utg_t *utg;
	gfa_aux_t aux;
} gfa_seg_t;

#define gfa_n_vtx(g) ((g)->n_seg << 1)

typedef struct {
	// segments
	uint32_t m_seg, n_seg, max_rank;
	gfa_seg_t *seg;
	void *h_names;
	// persistent names
	uint32_t m_pname, n_pname;
	char **pname;
	void *h_pnames;
	// links
	uint64_t m_arc, n_arc:62, is_srt:1, is_symm:1;
	gfa_arc_t *arc;
	gfa_aux_t *arc_aux;
	uint64_t *idx;
} gfa_t;

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

// shortest path

typedef struct {
	uint32_t v;
	int32_t target_dist;
	int32_t dist, n_path, path_end;
	int32_t meta;
} gfa_path_dst_t;

typedef struct {
	uint32_t v, d;
	int32_t pre;
} gfa_pathv_t;

// graph augmentation

typedef struct {
	uint32_t v[2];
	int32_t voff[2];
	int32_t coff[2], ctg;
} gfa_ins_t;

extern int gfa_verbose;
unsigned char gfa_comp_table[256];

#ifdef __cplusplus
extern "C" {
#endif

gfa_t *gfa_init(void);
int32_t gfa_add_seg(gfa_t *g, const char *name);
int32_t gfa_add_pname(gfa_t *g, const char *pname);
int32_t gfa_name2id(const gfa_t *g, const char *name);
uint64_t gfa_add_arc1(gfa_t *g, uint32_t v, uint32_t w, int32_t ov, int32_t ow, int64_t link_id, int comp);
void gfa_cleanup(gfa_t *g); // permanently delete arcs marked as deleted, sort and then index
void gfa_finalize(gfa_t *g);
void gfa_destroy(gfa_t *g);

gfa_t *gfa_read(const char *fn);
void gfa_print(const gfa_t *g, FILE *fp, int M_only);

void gfa_symm(gfa_t *g); // delete multiple edges and restore skew-symmetry
int gfa_arc_del_trans(gfa_t *g, int fuzz); // transitive reduction
int gfa_arc_del_short(gfa_t *g, float drop_ratio); // delete short arcs
int gfa_cut_tip(gfa_t *g, int max_ext); // cut tips
int gfa_cut_internal(gfa_t *g, int max_ext); // drop internal segments
int gfa_cut_biloop(gfa_t *g, int max_ext); // Hmm... I forgot... Some type of weird local topology
int gfa_pop_bubble(gfa_t *g, int max_dist); // bubble popping
gfa_t *gfa_ug_gen(const gfa_t *g);

uint8_t *gfa_aux_get(int l_data, const uint8_t *data, const char tag[2]);
int gfa_aux_del(int l_data, uint8_t *data, uint8_t *s);

void gfa_sub(gfa_t *g, int n, char *const* seg, int step);

gfa_sub_t *gfa_sub_from(void *km0, const gfa_t *g, uint32_t v0, int32_t max_dist);
void gfa_sub_destroy(gfa_sub_t *sub);
void gfa_sub_print(FILE *fp, const gfa_t *g, const gfa_sub_t *sub);
gfa_pathv_t *gfa_shortest_k(void *km0, const gfa_t *g, uint32_t src, int32_t n_dst, gfa_path_dst_t *dst, int32_t max_dist, int32_t max_k, int32_t *n_pathv);
void gfa_sub_print_path(FILE *fp, const gfa_t *g, int32_t n, gfa_pathv_t *path);

int gfa_ins_adj(const gfa_t *g, int min_len, gfa_ins_t *ins, const char *seq);
int32_t gfa_ins_filter(const gfa_t *g, int32_t n_ins, gfa_ins_t *ins);
void gfa_augment(gfa_t *g, int32_t n_ins, const gfa_ins_t *ins, int32_t n_ctg, const char *const* name, const char *const* seq);

#ifdef __cplusplus
}
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline void gfa_arc_del(gfa_t *g, uint32_t v, uint32_t w, int del)
{
	uint32_t i, nv = gfa_arc_n(g, v);
	gfa_arc_t *av = gfa_arc_a(g, v);
	for (i = 0; i < nv; ++i)
		if (av[i].w == w) av[i].del = !!del;
}

static inline void gfa_seg_del(gfa_t *g, uint32_t s)
{
	uint32_t k;
	g->seg[s].del = 1;
	for (k = 0; k < 2; ++k) {
		uint32_t i, v = s<<1 | k;
		uint32_t nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			av[i].del = 1;
			gfa_arc_del(g, av[i].w^1, v^1, 1);
		}
	}
}

#endif
