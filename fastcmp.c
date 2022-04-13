#include <assert.h>
#include <string.h>
#include "mgpriv.h"
#include "algo.h"
#include "kalloc.h"
#define HAVE_KALLOC
#include "ksw2.h"

int32_t mg_nwcmp(void *km, int32_t l1, const char *s1, int32_t l2, const char *s2, const int8_t mat[25], int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2,
				 int32_t bw, int32_t *mlen_, int32_t *blen_)
{
	extern unsigned char seq_nt4_table[256]; // in sketch.c
	int32_t i, mlen, blen, x1, x2;
	ksw_extz_t ez;
	uint8_t *t1, *t2;
	*mlen_ = 0, *blen_ = l1 > l2? l1 : l2;
	if (l1 == 0 || l2 == 0) return 0;
	memset(&ez, 0, sizeof(ez));
	KMALLOC(km, t1, l1);
	KMALLOC(km, t2, l2);
	for (i = 0; i < l1; ++i)
		t1[i] = seq_nt4_table[(uint8_t)s1[i]];
	for (i = 0; i < l2; ++i)
		t2[i] = seq_nt4_table[(uint8_t)s2[i]];
	ksw_extd2_sse(km, l1, t1, l2, t2, 5, mat, gapo, gape, gapo2, gape2, bw, 1000, 0, 0, &ez); // 1000 for zdrop
	if (!ez.zdropped) {
		for (i = x1 = x2 = mlen = blen = 0; i < ez.n_cigar; ++i) {
			int32_t j, op = ez.cigar[i]&0xf, len = ez.cigar[i]>>4;
			if (op == 0) {
				for (j = 0; j < len; ++j)
					if (t1[x1+j] == t2[x2+j] && t1[x1+j] < 4)
						++mlen;
				x1 += len, x2 += len, blen += len;
			} else if (op == 1) blen += len, x1 += len;
			else if (op == 2) blen += len, x2 += len;
		}
		assert(x1 == l1 && x2 == l2);
	} else mlen = 0, blen = l1 > l2? l1 : l2;
	kfree(km, t1);
	kfree(km, t2);
	kfree(km, ez.cigar);
	*mlen_ = mlen, *blen_ = blen;
	return ez.score == KSW_NEG_INF? INT32_MIN : ez.score;
}
