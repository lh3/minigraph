#include <assert.h>
#include <string.h>
#include "mgpriv.h"
#include "algo.h"
#include "kalloc.h"
#define HAVE_KALLOC
#include "ksw2.h"

static int32_t mg_fc_kmer(int32_t len, const char *seq, int32_t rid, int32_t k, mg128_t *a)
{
	int32_t i, l, n;
	uint64_t x, mask = (1ULL<<k*2) - 1;
	for (i = l = 0, x = 0, n = 0; i < len; ++i) {
		int32_t c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) {
			x = (x << 2 | c) & mask;
			if (++l >= k) a[n].x = x<<1 | rid, a[n].y = i, ++n;
		} else l = 0, x = 0;
	}
	return n;
}

int32_t mg_anchor2mlen(void *km, int32_t k, int32_t n_b, uint64_t *b, int32_t *st_high, int32_t *en_high)
{
	int32_t i, mlen, n_lis, *lis;
	if (st_high) *st_high = -1;
	if (en_high) *en_high = -1;
	if (n_b == 0) return 0;
	radix_sort_gfa64(b, b + n_b);
	for (i = 0; i < n_b; ++i)
		b[i] = b[i]>>32 | b[i]<<32;
	KMALLOC(km, lis, n_b);
	n_lis = mg_lis_64(km, n_b, b, lis);
	for (i = 1, mlen = k; i < n_lis; ++i) {
		int32_t ll2 = (int32_t)(b[lis[i]]>>32) - (int32_t)(b[lis[i-1]]>>32);
		int32_t ll1 = (int32_t)b[lis[i]] - (int32_t)b[lis[i-1]];
		mlen += ll1 > k && ll2 > k? k : ll1 < ll2? ll1 : ll2;
	}
	if (st_high) *st_high = (uint32_t)b[lis[0]] + 1 - k;
	if (en_high) *en_high = (uint32_t)b[lis[n_lis-1]] + 1;
	kfree(km, lis);
	return mlen;
}

int32_t mg_fastcmp(void *km, int32_t l1, const char *s1, int32_t l2, const char *s2, int32_t k, int32_t max_occ)
{
	int32_t i, n_a, n_b, m_b, i0, mlen;
	mg128_t *a;
	uint64_t *b;

	if (l1 < k || l2 < k) return 0;

	KMALLOC(km, a, l1 + l2);
	n_a = mg_fc_kmer(l1, s1, 0, k, a);
	n_a += mg_fc_kmer(l2, s2, 1, k, &a[n_a]);
	radix_sort_128x(a, a + n_a);

	n_b = m_b = 0, b = 0;
	for (i0 = 0, i = 1; i <= n_a; ++i) {
		if (i == n_a || a[i0].x>>1 != a[i].x>>1) {
			if (i - i0 >= 2) {
				int32_t j, s, t;
				for (j = i0; j < i && (a[j].x&1) == 0; ++j) {}
				if (j > i0 && j < i && j - i0 <= max_occ && i - j <= max_occ) {
					for (s = i0; s < j; ++s)
						for (t = j; t < i; ++t) {
							if (n_b == m_b) KEXPAND(km, b, m_b);
							b[n_b++] = (uint64_t)a[s].y<<32 | a[t].y;
						}
				}
			}
			i0 = i;
		}
	}
	kfree(km, a);
	mlen = mg_anchor2mlen(km, k, n_b, b, 0, 0);
	kfree(km, b);
	return mlen;
}

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
