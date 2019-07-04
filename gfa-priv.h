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

#ifdef __cplusplus
extern "C" {
#endif

char *gfa_strdup(const char *src);
void gfa_arc_sort(gfa_t *g);
void gfa_arc_index(gfa_t *g);
uint32_t gfa_fix_symm(gfa_t *g);
int32_t gfa_pseq_add(gfa_t *g, const char *pname);
void gfa_pseq_update(gfa_t *g, const gfa_seg_t *s);

#ifdef __cplusplus
}
#endif

#endif // ~__GFA_PRIV_H__
