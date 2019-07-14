#include "mgpriv.h"
#include "algo.h"

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

int32_t mg_fastcmp(void *km, int32_t l1, const char *s1, int32_t l2, const char *s2, int32_t k, int32_t max_occ)
{
	int32_t i, n_a, n_b, m_b, i0;
	mg128_t *a;
	uint64_t *b;
	int32_t n_lis, *lis, mlen;

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

	radix_sort_64(b, b + n_b);
	for (i = 0; i < n_b; ++i)
		b[i] = b[i]>>32 | b[i]<<32;
	KMALLOC(km, lis, n_b);
	n_lis = ks_lis_64(km, n_b, b, lis);

	mlen = k;
	for (i = 1; i < n_lis; ++i) {
		int32_t ll2 = (b[lis[i]]>>32) - (b[lis[i-1]]>>32);
		int32_t ll1 = (int32_t)b[lis[i]] - (int32_t)b[lis[i-1]];
		mlen += ll1 > k && ll2 > k? k : ll1 < ll2? ll1 : ll2;
	}

	kfree(km, b);
	kfree(km, lis);
	return mlen;
}
