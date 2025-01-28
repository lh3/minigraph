#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "kalloc.h"
#include "mgpriv.h"
#include "khashl.h"
#include "sys.h"

struct mg_tbuf_s {
	void *km;
	int frag_gap;
};

mg_tbuf_t *mg_tbuf_init(void)
{
	mg_tbuf_t *b;
	b = (mg_tbuf_t*)calloc(1, sizeof(mg_tbuf_t));
	if (!(mg_dbg_flag & MG_DBG_NO_KALLOC)) b->km = km_init();
	return b;
}

void mg_tbuf_destroy(mg_tbuf_t *b)
{
	if (b == 0) return;
	if (b->km) km_destroy(b->km);
	free(b);
}

void *mg_tbuf_get_km(mg_tbuf_t *b)
{
	return b->km;
}

static void collect_minimizers(void *km, const mg_mapopt_t *opt, const mg_idx_t *gi, int n_segs, const int *qlens, const char **seqs, mg128_v *mv)
{
	int i, n, sum = 0;
	mv->n = 0;
	for (i = n = 0; i < n_segs; ++i) {
		size_t j;
		mg_sketch(km, seqs[i], qlens[i], gi->w, gi->k, i, mv);
		for (j = n; j < mv->n; ++j)
			mv->a[j].y += sum << 1;
		sum += qlens[i], n = mv->n;
	}
}

#include "ksort.h"
#define heap_lt(a, b) ((a).x > (b).x)
KSORT_INIT(heap, mg128_t, heap_lt)

typedef struct {
	uint32_t n;
	uint32_t q_pos, q_span;
	uint32_t seg_id:31, is_tandem:1;
	const uint64_t *cr;
} mg_match_t;

static mg_match_t *collect_matches(void *km, int *_n_m, int max_occ, const mg_idx_t *gi, const mg128_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos, int32_t **mini_pos)
{
	int rep_st = 0, rep_en = 0, n_m;
	size_t i;
	mg_match_t *m;
	*n_mini_pos = 0;
	KMALLOC(km, *mini_pos, mv->n);
	m = (mg_match_t*)kmalloc(km, mv->n * sizeof(mg_match_t));
	for (i = 0, n_m = 0, *rep_len = 0, *n_a = 0; i < mv->n; ++i) {
		const uint64_t *cr;
		mg128_t *p = &mv->a[i];
		uint32_t q_pos = (uint32_t)p->y, q_span = p->x & 0xff;
		int t;
		cr = mg_idx_get(gi, p->x>>8, &t);
		if (t >= max_occ) {
			int en = (q_pos >> 1) + 1, st = en - q_span;
			if (st > rep_en) {
				*rep_len += rep_en - rep_st;
				rep_st = st, rep_en = en;
			} else rep_en = en;
		} else {
			mg_match_t *q = &m[n_m++];
			q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> 32;
			q->is_tandem = 0;
			if (i > 0 && p->x>>8 == mv->a[i - 1].x>>8) q->is_tandem = 1;
			if (i < mv->n - 1 && p->x>>8 == mv->a[i + 1].x>>8) q->is_tandem = 1;
			*n_a += q->n;
			(*mini_pos)[(*n_mini_pos)++] = q_pos>>1;
		}
	}
	*rep_len += rep_en - rep_st;
	*_n_m = n_m;
	return m;
}

static mg128_t *collect_seed_hits_heap(void *km, const mg_mapopt_t *opt, int max_occ, const mg_idx_t *gi, const char *qname, const mg128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, int32_t **mini_pos)
{
	int i, n_m, heap_size = 0;
	int64_t n_for = 0, n_rev = 0;
	mg_match_t *m;
	mg128_t *a, *heap;

	m = collect_matches(km, &n_m, max_occ, gi, mv, n_a, rep_len, n_mini_pos, mini_pos);

	heap = (mg128_t*)kmalloc(km, n_m * sizeof(mg128_t));
	a = (mg128_t*)kmalloc(km, *n_a * sizeof(mg128_t));

	for (i = 0, heap_size = 0; i < n_m; ++i) {
		if (m[i].n > 0) {
			heap[heap_size].x = m[i].cr[0];
			heap[heap_size].y = (uint64_t)i<<32;
			++heap_size;
		}
	}
	ks_heapmake_heap(heap_size, heap);
	while (heap_size > 0) {
		mg_match_t *q = &m[heap->y>>32];
		mg128_t *p;
		uint64_t r = heap->x;
		int32_t rpos = (uint32_t)r >> 1;
		// TODO: skip anchor if MG_F_NO_DIAL
		if ((r&1) == (q->q_pos&1)) { // forward strand
			p = &a[n_for++];
			p->x = r>>32<<33 | rpos;
		} else { // reverse strand; TODO: more testing needed for this block
			p = &a[(*n_a) - (++n_rev)];
			p->x = r>>32<<33 | 1ULL<<32 | (gi->g->seg[r>>32].len - (rpos + 1 - q->q_span) - 1);
		}
		p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
		p->y |= (uint64_t)q->seg_id << MG_SEED_SEG_SHIFT;
		if (q->is_tandem) p->y |= MG_SEED_TANDEM;
		p->y |= (uint64_t)(q->n < 255? q->n : 255) << MG_SEED_OCC_SHIFT;
		// update the heap
		if ((uint32_t)heap->y < q->n - 1) {
			++heap[0].y;
			heap[0].x = m[heap[0].y>>32].cr[(uint32_t)heap[0].y];
		} else {
			heap[0] = heap[heap_size - 1];
			--heap_size;
		}
		ks_heapdown_heap(0, heap_size, heap);
	}
	kfree(km, m);
	kfree(km, heap);

	// reverse anchors on the reverse strand, as they are in the descending order
	if (*n_a > n_for + n_rev) {
		memmove(a + n_for, a + (*n_a) - n_rev, n_rev * sizeof(mg128_t));
		*n_a = n_for + n_rev;
	}
	return a;
}

static mg128_t *collect_seed_hits(void *km, const mg_mapopt_t *opt, int max_occ, const mg_idx_t *gi, const char *qname, const mg128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, int32_t **mini_pos)
{
	int i, n_m;
	mg_match_t *m;
	mg128_t *a;
	m = collect_matches(km, &n_m, max_occ, gi, mv, n_a, rep_len, n_mini_pos, mini_pos);
	a = (mg128_t*)kmalloc(km, *n_a * sizeof(mg128_t));
	for (i = 0, *n_a = 0; i < n_m; ++i) {
		mg_match_t *q = &m[i];
		const uint64_t *r = q->cr;
		uint32_t k;
		for (k = 0; k < q->n; ++k) {
			int32_t rpos = (uint32_t)r[k] >> 1;
			mg128_t *p;
			if (qname && (opt->flag & MG_M_NO_DIAG)) {
				const gfa_seg_t *s = &gi->g->seg[r[k]>>32];
				const char *gname = s->snid >= 0 && gi->g->sseq? gi->g->sseq[s->snid].name : s->name;
				int32_t g_pos;
				if (s->snid >= 0 && gi->g->sseq)
					gname = gi->g->sseq[s->snid].name, g_pos = s->soff + (uint32_t)r[k];
				else
					gname = s->name, g_pos = (uint32_t)r[k];
				if (g_pos == q->q_pos && strcmp(qname, gname) == 0)
					continue;
			}
			p = &a[(*n_a)++];
			if ((r[k]&1) == (q->q_pos&1)) // forward strand
				p->x = r[k]>>32<<33 | rpos;
			else // reverse strand
				p->x = r[k]>>32<<33 | 1ULL<<32 | (gi->g->seg[r[k]>>32].len - (rpos + 1 - q->q_span) - 1);
			p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			p->y |= (uint64_t)q->seg_id << MG_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MG_SEED_TANDEM;
			p->y |= (uint64_t)(q->n < 255? q->n : 255) << MG_SEED_OCC_SHIFT;
		}
	}
	kfree(km, m);
	radix_sort_128x(a, a + (*n_a));
	return a;
}

static void mm_fix_bad_ends(const mg128_t *a, int32_t lc_max_occ, int32_t lc_max_trim, int32_t *as, int32_t *cnt)
{
	int32_t i, k, as0 = *as, cnt0 = *cnt;
	for (i = as0 + cnt0 - 1, k = 0; k < lc_max_trim && k < cnt0; ++k, --i)
		if (a[i].y>>MG_SEED_OCC_SHIFT <= lc_max_occ)
			break;
	*cnt -= k;
	for (i = as0, k = 0; k < *cnt && k < lc_max_trim; ++i, ++k)
		if (a[i].y>>MG_SEED_OCC_SHIFT <= lc_max_occ)
			break;
	*as += k, *cnt -= k;
}

static void mm_fix_bad_ends_alt(const mg128_t *a, int32_t score, int bw, int min_match, int32_t *as, int32_t *cnt)
{
	int32_t i, l, m, as0 = *as, cnt0 = *cnt;
	if (cnt0 < 3) return;
	m = l = a[as0].y >> 32 & 0xff;
	for (i = as0 + 1; i < as0 + cnt0 - 1; ++i) {
		int32_t lq, lr, min, max;
		int32_t q_span = a[i].y >> 32 & 0xff;
		lr = (int32_t)a[i].x - (int32_t)a[i-1].x;
		lq = (int32_t)a[i].y - (int32_t)a[i-1].y;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) *as = i;
		l += min;
		m += min < q_span? min : q_span;
		if (l >= bw << 1 || (m >= min_match && m >= bw) || m >= score>>1) break;
	}
	*cnt = as0 + cnt0 - *as;
	m = l = a[as0 + cnt0 - 1].y >> 32 & 0xff;
	for (i = as0 + cnt0 - 2; i > *as; --i) {
		int32_t lq, lr, min, max;
		int32_t q_span = a[i+1].y >> 32 & 0xff;
		lr = (int32_t)a[i+1].x - (int32_t)a[i].x;
		lq = (int32_t)a[i+1].y - (int32_t)a[i].y;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) *cnt = i + 1 - *as;
		l += min;
		m += min < q_span? min : q_span;
		if (l >= bw << 1 || (m >= min_match && m >= bw) || m >= score>>1) break;
	}
}

static int *collect_long_gaps(void *km, int as1, int cnt1, mg128_t *a, int min_gap, int *n_)
{
	int i, n, *K;
	*n_ = 0;
	for (i = 1, n = 0; i < cnt1; ++i) { // count the number of gaps longer than min_gap
		int gap = ((int32_t)a[as1 + i].y - a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - a[as1 + i - 1].x);
		if (gap < -min_gap || gap > min_gap) ++n;
	}
	if (n <= 1) return 0;
	K = (int*)kmalloc(km, n * sizeof(int));
	for (i = 1, n = 0; i < cnt1; ++i) { // store the positions of long gaps
		int gap = ((int32_t)a[as1 + i].y - a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - a[as1 + i - 1].x);
		if (gap < -min_gap || gap > min_gap)
			K[n++] = i;
	}
	*n_ = n;
	return K;
}

static void mm_filter_bad_seeds(void *km, int as1, int cnt1, mg128_t *a, int min_gap, int diff_thres, int max_ext_len, int max_ext_cnt)
{
	int max_st, max_en, n, i, k, max, *K;
	K = collect_long_gaps(km, as1, cnt1, a, min_gap, &n);
	if (K == 0) return;
	max = 0, max_st = max_en = -1;
	for (k = 0;; ++k) { // traverse long gaps
		int gap, l, n_ins = 0, n_del = 0, qs, rs, max_diff = 0, max_diff_l = -1;
		if (k == n || k >= max_en) {
			if (max_en > 0)
				for (i = K[max_st]; i < K[max_en]; ++i)
					a[as1 + i].y |= MG_SEED_IGNORE;
			max = 0, max_st = max_en = -1;
			if (k == n) break;
		}
		i = K[k];
		gap = ((int32_t)a[as1 + i].y - (int32_t)a[as1 + i - 1].y) - (int32_t)(a[as1 + i].x - a[as1 + i - 1].x);
		if (gap > 0) n_ins += gap;
		else n_del += -gap;
		qs = (int32_t)a[as1 + i - 1].y;
		rs = (int32_t)a[as1 + i - 1].x;
		for (l = k + 1; l < n && l <= k + max_ext_cnt; ++l) {
			int j = K[l], diff;
			if ((int32_t)a[as1 + j].y - qs > max_ext_len || (int32_t)a[as1 + j].x - rs > max_ext_len) break;
			gap = ((int32_t)a[as1 + j].y - (int32_t)a[as1 + j - 1].y) - (int32_t)(a[as1 + j].x - a[as1 + j - 1].x);
			if (gap > 0) n_ins += gap;
			else n_del += -gap;
			diff = n_ins + n_del - abs(n_ins - n_del);
			if (max_diff < diff)
				max_diff = diff, max_diff_l = l;
		}
		if (max_diff > diff_thres && max_diff > max)
			max = max_diff, max_st = k, max_en = max_diff_l;
	}
	kfree(km, K);
}

static void mm_filter_bad_seeds_alt(void *km, int as1, int cnt1, mg128_t *a, int min_gap, int max_ext)
{
	int n, k, *K;
	K = collect_long_gaps(km, as1, cnt1, a, min_gap, &n);
	if (K == 0) return;
	for (k = 0; k < n;) {
		int i = K[k], l;
		int gap1 = ((int32_t)a[as1 + i].y - (int32_t)a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - (int32_t)a[as1 + i - 1].x);
		int re1 = (int32_t)a[as1 + i].x;
		int qe1 = (int32_t)a[as1 + i].y;
		gap1 = gap1 > 0? gap1 : -gap1;
		for (l = k + 1; l < n; ++l) {
			int j = K[l], gap2, q_span_pre, rs2, qs2, m;
			if ((int32_t)a[as1 + j].y - qe1 > max_ext || (int32_t)a[as1 + j].x - re1 > max_ext) break;
			gap2 = ((int32_t)a[as1 + j].y - (int32_t)a[as1 + j - 1].y) - (int32_t)(a[as1 + j].x - a[as1 + j - 1].x);
			q_span_pre = a[as1 + j - 1].y >> 32 & 0xff;
			rs2 = (int32_t)a[as1 + j - 1].x + q_span_pre;
			qs2 = (int32_t)a[as1 + j - 1].y + q_span_pre;
			m = rs2 - re1 < qs2 - qe1? rs2 - re1 : qs2 - qe1;
			gap2 = gap2 > 0? gap2 : -gap2;
			if (m > gap1 + gap2) break;
			re1 = (int32_t)a[as1 + j].x;
			qe1 = (int32_t)a[as1 + j].y;
			gap1 = gap2;
		}
		if (l > k + 1) {
			int j, end = K[l - 1];
			for (j = K[k]; j < end; ++j)
				a[as1 + j].y |= MG_SEED_IGNORE;
			a[as1 + end].y |= MG_SEED_FIXED;
		}
		k = l;
	}
	kfree(km, K);
}

static double print_time(double t0, int stage, const char *qname)
{
	double t;
	t = realtime();
	fprintf(stderr, "Q%d\t%s\t%.3f\n", stage, qname, t - t0);
	return t;
}

void mg_map_frag(const mg_idx_t *gi, int n_segs, const int *qlens, const char **seqs, mg_gchains_t **gcs, mg_tbuf_t *b, const mg_mapopt_t *opt, const char *qname)
{
	int i, l, rep_len, qlen_sum, n_lc, n_gc, n_mini_pos;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MG_M_SPLICE), is_sr = !!(opt->flag & MG_M_SR);
	uint32_t hash;
	int64_t n_a;
	uint64_t *u;
	int32_t *mini_pos;
	mg128_t *a;
	mg128_v mv = {0,0,0};
	mg_lchain_t *lc;
	char *seq_cat;
	km_stat_t kmst;
	float tmp, chn_pen_gap, chn_pen_skip;
	double t = 0.0;

	for (i = 0, qlen_sum = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], gcs[i] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MG_MAX_SEG) return;
	if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen) return;

	hash  = qname? kh_hash_str(qname) : 0;
	hash ^= kh_hash_uint32(qlen_sum) + kh_hash_uint32(opt->seed);
	hash  = kh_hash_uint32(hash);

	collect_minimizers(b->km, opt, gi, n_segs, qlens, seqs, &mv);
	if (opt->flag & MG_M_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->occ_max1, gi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	else a = collect_seed_hits(b->km, opt, opt->occ_max1, gi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);

	if (mg_dbg_flag & MG_DBG_SEED) {
		fprintf(stderr, "RS\t%d\n", rep_len);
		for (i = 0; i < n_a; ++i)
			fprintf(stderr, "SD\t%s\t%d\t%c\t%d\t%d\t%d\n", gi->g->seg[a[i].x>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>32&1], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
					i == 0? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));
	}

	// set max chaining gap on the query and the reference sequence
	if (is_sr)
		max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
	else max_chain_gap_qry = opt->max_gap;
	if (opt->max_gap_ref > 0) {
		max_chain_gap_ref = opt->max_gap_ref; // always honor mg_mapopt_t::max_gap_ref if set
	} else if (opt->max_frag_len > 0) {
		max_chain_gap_ref = opt->max_frag_len - qlen_sum;
		if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
	} else max_chain_gap_ref = opt->max_gap;

	tmp = expf(-opt->div * gi->k);
	chn_pen_gap = opt->chn_pen_gap * tmp;
	chn_pen_skip = opt->chn_pen_skip * tmp;

	if (mg_dbg_flag & MG_DBG_QNAME) t = realtime();
	if (n_a == 0) {
		if (a) kfree(b->km, a);
		a = 0, n_lc = 0, u = 0;
	} else {
		if (opt->flag & MG_M_RMQ) {
			a = mg_lchain_rmq(opt->max_gap, opt->max_gap_pre, opt->bw, opt->max_lc_skip, opt->rmq_size_cap, opt->min_lc_cnt, opt->min_lc_score,
							  chn_pen_gap, chn_pen_skip, n_a, a, &n_lc, &u, b->km);
		} else {
			a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_lc_skip, opt->max_lc_iter, opt->min_lc_cnt, opt->min_lc_score,
							 chn_pen_gap, chn_pen_skip, is_splice, n_segs, n_a, a, &n_lc, &u, b->km);
		}
	}
	if (mg_dbg_flag & MG_DBG_QNAME) t = print_time(t, 1, qname);

	if (opt->bw_long > opt->bw && (opt->flag & (MG_M_SPLICE|MG_M_SR)) == 0 && n_segs == 1 && n_lc > 1) { // re-chain/long-join for long sequences
		int32_t st = (int32_t)a[0].y, en = (int32_t)a[(int32_t)u[0] - 1].y;
		if (qlen_sum - (en - st) > opt->rmq_rescue_size || qlen_sum - (en - st) > qlen_sum * opt->rmq_rescue_ratio) {
			int32_t i;
			for (i = 0, n_a = 0; i < n_lc; ++i) n_a += (int32_t)u[i];
			kfree(b->km, u);
			radix_sort_128x(a, a + n_a);
			a = mg_lchain_rmq(opt->max_gap, opt->max_gap_pre, opt->bw_long, opt->max_lc_skip, opt->rmq_size_cap, opt->min_lc_cnt, opt->min_lc_score,
							  chn_pen_gap, chn_pen_skip, n_a, a, &n_lc, &u, b->km);
		}
	}

	b->frag_gap = max_chain_gap_ref;
	kfree(b->km, mv.a);

	if (n_lc) {
		lc = mg_lchain_gen(b->km, hash, qlen_sum, n_lc, u, a);
		if (n_lc > 1) {
			int32_t n_lc_new = 0;
			for (i = 0; i < n_lc; ++i) {
				mg_lchain_t *p = &lc[i];
				int32_t cnt = p->cnt, off = p->off;
				mm_fix_bad_ends(a, opt->lc_max_occ, opt->lc_max_trim, &off, &cnt);
				mm_fix_bad_ends_alt(a, p->score, opt->bw, 100, &off, &cnt);
				mm_filter_bad_seeds(b->km, off, cnt, a, 10, 40, opt->max_gap>>1, 10);
				mm_filter_bad_seeds_alt(b->km, off, cnt, a, 30, opt->max_gap>>1);
				//printf("X\t%d\t%d\t%d\t%d\t%d\t%d\n", p->qs, p->qe, p->off, p->cnt, off, cnt);
				p->off = off, p->cnt = cnt;
				if (cnt >= opt->min_lc_cnt) {
					int32_t q_span = a[p->off].y>>32 & 0xff;
					p->rs = (int32_t)a[p->off].x + 1 - q_span;
					p->qs = (int32_t)a[p->off].y + 1 - q_span;
					p->re = (int32_t)a[p->off + p->cnt - 1].x + 1;
					p->qe = (int32_t)a[p->off + p->cnt - 1].y + 1;
					lc[n_lc_new++] = *p;
				}
			}
			n_lc = n_lc_new;
		}
		for (i = 0; i < n_lc; ++i)
			mg_update_anchors(lc[i].cnt, &a[lc[i].off], n_mini_pos, mini_pos);
	} else lc = 0;
	kfree(b->km, mini_pos);
	kfree(b->km, u);
	if (mg_dbg_flag & MG_DBG_QNAME) t = print_time(t, 2, qname);

	if (mg_dbg_flag & MG_DBG_LCHAIN)
		mg_print_lchain(stdout, gi, n_lc, lc, a, qname);

	KMALLOC(b->km, seq_cat, qlen_sum);
	for (i = l = 0; i < n_segs; ++i) {
		strncpy(&seq_cat[l], seqs[i], qlens[i]);
		l += qlens[i];
	}
	n_gc = mg_gchain1_dp(b->km, gi->g, &n_lc, lc, qlen_sum, opt->bw_long, opt->bw_long, opt->bw_long, opt->max_gc_skip, opt->ref_bonus,
						 chn_pen_gap, chn_pen_skip, opt->mask_level, a, &u);
	if (mg_dbg_flag & MG_DBG_QNAME) t = print_time(t, 3, qname);
	gcs[0] = mg_gchain_gen(0, b->km, gi->g, gi->es, n_gc, u, lc, a, hash, opt->min_gc_cnt, opt->min_gc_score, opt->gdp_max_ed, n_segs, seq_cat);
	if (mg_dbg_flag & MG_DBG_QNAME) t = print_time(t, 4, qname);
	gcs[0]->rep_len = rep_len;
	kfree(b->km, a);
	kfree(b->km, lc);
	kfree(b->km, u);

	mg_gchain_set_parent(b->km, opt->mask_level, gcs[0]->n_gc, gcs[0]->gc, opt->sub_diff, 0);
	mg_gchain_flt_sub(opt->pri_ratio, gi->k * 2, opt->best_n, gcs[0]->n_gc, gcs[0]->gc);
	mg_gchain_drop_flt(b->km, gcs[0]);
	mg_gchain_set_mapq(b->km, gcs[0], qlen_sum, mv.n, opt->min_gc_score);
	if ((opt->flag&MG_M_CIGAR) && n_segs == 1) {
		mg_gchain_cigar(b->km, gi->g, gi->es, seq_cat, gcs[0], qname);
		mg_gchain_gen_ds(b->km, gi->g, gi->es, seq_cat, gcs[0]);
	}
	kfree(b->km, seq_cat);
	if (mg_dbg_flag & MG_DBG_QNAME) t = print_time(t, 5, qname);

	if (b->km) {
		km_stat(b->km, &kmst);
		if (mg_dbg_flag & MG_DBG_QNAME)
			fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
		if (kmst.n_blocks != kmst.n_cores) {
			fprintf(stderr, "[E::%s] memory leak at %s\n", __func__, qname);
			abort();
		}
		if (kmst.largest > 1U<<28 || (opt->cap_kalloc > 0 && kmst.capacity > opt->cap_kalloc)) {
			km_destroy(b->km);
			b->km = km_init();
		}
	}
}

mg_gchains_t *mg_map(const mg_idx_t *gi, int qlen, const char *seq, mg_tbuf_t *b, const mg_mapopt_t *opt, const char *qname)
{
	mg_gchains_t *gcs;
	mg_map_frag(gi, 1, &qlen, &seq, &gcs, b, opt, qname);
	return gcs;
}
