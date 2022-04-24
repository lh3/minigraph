#ifndef __GFA_H__
#define __GFA_H__

#include <stdio.h>
#include <stdint.h>

#define GFA_VERSION "0.5-r247-dirty"

#define GFA_O_OV_EXT   0x1
#define GFA_O_NO_SEQ   0x2

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
	uint32_t w;
	int32_t rank;
	int32_t ov, ow;
	uint64_t link_id:61, strong:1, del:1, comp:1; // link_id: a pair of dual arcs are supposed to have the same link_id
} gfa_arc_t;

#define gfa_arc_head(a) ((uint32_t)((a).v_lv>>32))
#define gfa_arc_tail(a) ((a).w)
#define gfa_arc_len(a) ((uint32_t)(a).v_lv) // different from the original string graph
#define gfa_arc_lw(g, a) ((g)->seg[(a).w>>1].len - (a).ow)

#define gfa_arc_n(g, v) ((uint32_t)(g)->idx[(v)])
#define gfa_arc_a(g, v) (&(g)->arc[(g)->idx[(v)]>>32])

typedef struct {
	uint32_t m_aux, l_aux;
	uint8_t *aux;
} gfa_aux_t;

typedef struct {
	uint32_t start, end; // start: starting vertex in the string graph; end: ending vertex
	uint32_t len_comp, dummy; // len_comp: the length of the complement unitig
	uint32_t m, n; // number of reads
	uint64_t *a; // list of reads
	uint64_t *r; // start and end on each read
	char **name;
} gfa_utg_t;

typedef struct {
	int32_t len;
	uint32_t del:16, circ:16;
	int32_t snid; // stable name ID
	int32_t soff; // stable start position
	int32_t rank; // stable rank
	char *name, *seq;
	gfa_utg_t *utg;
	gfa_aux_t aux;
} gfa_seg_t;

typedef struct {
	int32_t len, snid, soff, rank;
	uint64_t end[2];
	char *seq;
} gfa_sfa_t;

typedef struct {
	char *name;
	int32_t min, max, rank;
} gfa_sseq_t;

#define gfa_n_vtx(g) ((g)->n_seg << 1)

typedef struct {
	// segments
	uint32_t m_seg, n_seg, max_rank;
	gfa_seg_t *seg;
	void *h_names;
	// persistent names
	uint32_t m_sseq, n_sseq;
	gfa_sseq_t *sseq;
	void *h_snames;
	// links
	uint64_t m_arc, n_arc;
	gfa_arc_t *arc;
	gfa_aux_t *link_aux;
	uint64_t *idx;
} gfa_t;

typedef struct {
	const char *seq;
	int32_t len;
} gfa_edseq_t;

// graph augmentation

typedef struct {
	uint32_t v[2];
	int32_t voff[2];
	int32_t coff[2], ctg;
} gfa_ins_t;

extern int gfa_verbose;
extern unsigned char gfa_comp_table[256];

#ifdef __cplusplus
extern "C" {
#endif

gfa_t *gfa_init(void);
void gfa_destroy(gfa_t *g);
gfa_t *gfa_read(const char *fn);
void gfa_print(const gfa_t *g, FILE *fp, int M_only);

gfa_edseq_t *gfa_edseq_init(const gfa_t *g);
void gfa_edseq_destroy(int32_t n_seg, gfa_edseq_t *es);

int32_t gfa_name2id(const gfa_t *g, const char *name);
uint8_t *gfa_aux_get(int l_data, const uint8_t *data, const char tag[2]);
int gfa_aux_del(int l_data, uint8_t *data, uint8_t *s);

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
