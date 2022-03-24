#include <math.h>
#include <assert.h>
#include <string.h>
#include "mgpriv.h"
#include "kalloc.h"

// reorder gcs->a[] and gcs->lc[] such that they are in the same order as gcs->gc[]
void mg_gchain_restore_order(void *km, mg_gchains_t *gcs)
{
	int32_t i, n_a, n_lc;
	mg_llchain_t *lc;
	mg128_t *a;
	KMALLOC(km, lc, gcs->n_lc);
	KMALLOC(km, a, gcs->n_a);
	for (i = 0, n_a = n_lc = 0; i < gcs->n_gc; ++i) {
		mg_gchain_t *gc = &gcs->gc[i];
		assert(gc->cnt > 0);
		memcpy(&lc[n_lc], &gcs->lc[gc->off], gc->cnt * sizeof(mg_llchain_t));
		memcpy(&a[n_a], &gcs->a[gcs->lc[gc->off].off], gc->n_anchor * sizeof(mg128_t));
		n_lc += gc->cnt, n_a += gc->n_anchor;
	}
	memcpy(gcs->lc, lc, gcs->n_lc * sizeof(mg_llchain_t));
	memcpy(gcs->a, a, gcs->n_a * sizeof(mg128_t));
	kfree(km, lc); kfree(km, a);
	for (i = 0, n_lc = 0; i < gcs->n_gc; ++i) {
		mg_gchain_t *gc = &gcs->gc[i];
		gc->off = n_lc;
		n_lc += gc->cnt;
	}
	for (i = 0, n_a = 0; i < gcs->n_lc; ++i) {
		mg_llchain_t *lc = &gcs->lc[i];
		lc->off = n_a;
		n_a += lc->cnt;
	}
}

// recompute gcs->gc[].{off,n_anchor} and gcs->lc[].off, ASSUMING they are properly ordered (see mg_gchain_restore_order)
void mg_gchain_restore_offset(mg_gchains_t *gcs)
{
	int32_t i, j, n_a, n_lc;
	for (i = 0, n_a = n_lc = 0; i < gcs->n_gc; ++i) {
		mg_gchain_t *gc = &gcs->gc[i];
		gc->off = n_lc;
		for (j = 0, gc->n_anchor = 0; j < gc->cnt; ++j) {
			mg_llchain_t *lc = &gcs->lc[n_lc + j];
			lc->off = n_a;
			n_a += lc->cnt;
			gc->n_anchor += lc->cnt;
		}
		n_lc += gc->cnt;
	}
	assert(n_lc == gcs->n_lc && n_a == gcs->n_a);
}

// sort chains by score
void mg_gchain_sort_by_score(void *km, mg_gchains_t *gcs)
{
	mg128_t *z;
	mg_gchain_t *gc;
	int32_t i;
	KMALLOC(km, z, gcs->n_gc);
	KMALLOC(km, gc, gcs->n_gc);
	for (i = 0; i < gcs->n_gc; ++i)
		z[i].x = (uint64_t)gcs->gc[i].score << 32 | gcs->gc[i].hash, z[i].y = i;
	radix_sort_128x(z, z + gcs->n_gc);
	for (i = gcs->n_gc - 1; i >= 0; --i)
		gc[gcs->n_gc - 1 - i] = gcs->gc[z[i].y];
	memcpy(gcs->gc, gc, gcs->n_gc * sizeof(mg_gchain_t));
	kfree(km, z); kfree(km, gc);
	mg_gchain_restore_order(km, gcs); // this put gcs in the proper order
}

// set r[].{id,parent,subsc}, ASSUMING r[] is sorted by score
void mg_gchain_set_parent(void *km, float mask_level, int n, mg_gchain_t *r, int sub_diff, int hard_mask_level)
{
	int i, j, k, *w;
	uint64_t *cov;
	if (n <= 0) return;
	for (i = 0; i < n; ++i) r[i].id = i;
	cov = (uint64_t*)kmalloc(km, n * sizeof(uint64_t));
	w = (int*)kmalloc(km, n * sizeof(int));
	w[0] = 0, r[0].parent = 0;
	for (i = 1, k = 1; i < n; ++i) {
		mg_gchain_t *ri = &r[i];
		int si = ri->qs, ei = ri->qe, n_cov = 0, uncov_len = 0;
		if (hard_mask_level) goto skip_uncov;
		for (j = 0; j < k; ++j) { // traverse existing primary hits to find overlapping hits
			mg_gchain_t *rp = &r[w[j]];
			int sj = rp->qs, ej = rp->qe;
			if (ej <= si || sj >= ei) continue;
			if (sj < si) sj = si;
			if (ej > ei) ej = ei;
			cov[n_cov++] = (uint64_t)sj<<32 | ej;
		}
		if (n_cov == 0) {
			goto set_parent_test; // no overlapping primary hits; then i is a new primary hit
		} else if (n_cov > 0) { // there are overlapping primary hits; find the length not covered by existing primary hits
			int j, x = si;
			radix_sort_gfa64(cov, cov + n_cov);
			for (j = 0; j < n_cov; ++j) {
				if ((int)(cov[j]>>32) > x) uncov_len += (cov[j]>>32) - x;
				x = (int32_t)cov[j] > x? (int32_t)cov[j] : x;
			}
			if (ei > x) uncov_len += ei - x;
		}
skip_uncov:
		for (j = 0; j < k; ++j) { // traverse existing primary hits again
			mg_gchain_t *rp = &r[w[j]];
			int sj = rp->qs, ej = rp->qe, min, max, ol;
			if (ej <= si || sj >= ei) continue; // no overlap
			min = ej - sj < ei - si? ej - sj : ei - si;
			max = ej - sj > ei - si? ej - sj : ei - si;
			ol = si < sj? (ei < sj? 0 : ei < ej? ei - sj : ej - sj) : (ej < si? 0 : ej < ei? ej - si : ei - si); // overlap length; TODO: this can be simplified
			if ((float)ol / min - (float)uncov_len / max > mask_level) {
				int cnt_sub = 0;
				ri->parent = rp->parent;
				rp->subsc = rp->subsc > ri->score? rp->subsc : ri->score;
				if (ri->cnt >= rp->cnt) cnt_sub = 1;
				if (cnt_sub) ++rp->n_sub;
				break;
			}
		}
set_parent_test:
		if (j == k) w[k++] = i, ri->parent = i, ri->n_sub = 0;
	}
	kfree(km, cov);
	kfree(km, w);
}

// set r[].flt, i.e. mark weak suboptimal chains as filtered
int mg_gchain_flt_sub(float pri_ratio, int min_diff, int best_n, int n, mg_gchain_t *r)
{
	if (pri_ratio > 0.0f && n > 0) {
		int i, k, n_2nd = 0;
		for (i = k = 0; i < n; ++i) {
			int p = r[i].parent;
			if (p == i) { // primary
				r[i].flt = 0, ++k;
			} else if ((r[i].score >= r[p].score * pri_ratio || r[i].score + min_diff >= r[p].score) && n_2nd < best_n) {
				if (!(r[i].qs == r[p].qs && r[i].qe == r[p].qe && r[i].ps == r[p].ps && r[i].pe == r[p].pe)) // not identical hits; TODO: check path as well
					r[i].flt = 0, ++n_2nd, ++k;
				else r[i].flt = 1;
			} else r[i].flt = 1;
		}
		return k;
	}
	return n;
}

// hard drop filtered chains, ASSUMING gcs is properly ordered
void mg_gchain_drop_flt(void *km, mg_gchains_t *gcs)
{
	int32_t i, n_gc, n_lc, n_a, n_lc0, n_a0, *o2n;
	if (gcs->n_gc == 0) return;
	KMALLOC(km, o2n, gcs->n_gc);
	for (i = 0, n_gc = 0; i < gcs->n_gc; ++i) {
		mg_gchain_t *r = &gcs->gc[i];
		o2n[i] = -1;
		if (r->flt || r->cnt == 0) {
			kfree(gcs->km, r->p);
			continue;
		}
		o2n[i] = n_gc++;
	}
	n_gc = n_lc = n_a = 0;
	n_lc0 = n_a0 = 0;
	for (i = 0; i < gcs->n_gc; ++i) {
		mg_gchain_t *r = &gcs->gc[i];
		if (o2n[i] >= 0) {
			memmove(&gcs->a[n_a], &gcs->a[n_a0], r->n_anchor * sizeof(mg128_t));
			memmove(&gcs->lc[n_lc], &gcs->lc[n_lc0], r->cnt * sizeof(mg_llchain_t));
			gcs->gc[n_gc] = *r;
			gcs->gc[n_gc].id = n_gc;
			gcs->gc[n_gc].parent = o2n[gcs->gc[n_gc].parent];
			++n_gc, n_lc += r->cnt, n_a += r->n_anchor;
		}
		n_lc0 += r->cnt, n_a0 += r->n_anchor;
	}
	assert(n_lc0 == gcs->n_lc && n_a0 == gcs->n_a);
	kfree(km, o2n);
	gcs->n_gc = n_gc, gcs->n_lc = n_lc, gcs->n_a = n_a;
	if (n_a != n_a0) {
		KREALLOC(gcs->km, gcs->a, gcs->n_a);
		KREALLOC(gcs->km, gcs->lc, gcs->n_lc);
		KREALLOC(gcs->km, gcs->gc, gcs->n_gc);
	}
	mg_gchain_restore_offset(gcs);
}

// estimate mapping quality
void mg_gchain_set_mapq(void *km, mg_gchains_t *gcs, int qlen, int max_mini, int min_gc_score)
{
	static const float q_coef = 40.0f;
	int64_t sum_sc = 0;
	float uniq_ratio, r_sc, r_cnt;
	int i, t_sc, t_cnt;
	if (gcs == 0 || gcs->n_gc == 0) return;
	t_sc = qlen < 100? qlen : 100;
	t_cnt = max_mini < 10? max_mini : 10;
	if (t_cnt < 5) t_cnt = 5;
	r_sc = 1.0 / t_sc;
	r_cnt = 1.0 / t_cnt;
	for (i = 0; i < gcs->n_gc; ++i)
		if (gcs->gc[i].parent == gcs->gc[i].id)
			sum_sc += gcs->gc[i].score;
	uniq_ratio = (float)sum_sc / (sum_sc + gcs->rep_len);
	for (i = 0; i < gcs->n_gc; ++i) {
		mg_gchain_t *r = &gcs->gc[i];
		if (r->parent == r->id) {
			int mapq, subsc;
			float pen_s1 = (r->score > t_sc? 1.0f : r->score * r_sc) * uniq_ratio;
			float x, pen_cm = r->n_anchor > t_cnt? 1.0f : r->n_anchor * r_cnt;
			pen_cm = pen_s1 < pen_cm? pen_s1 : pen_cm;
			subsc = r->subsc > min_gc_score? r->subsc : min_gc_score;
			x = (float)subsc / r->score;
			mapq = (int)(pen_cm * q_coef * (1.0f - x) * logf(r->score));
			mapq -= (int)(4.343f * logf(r->n_sub + 1) + .499f);
			mapq = mapq > 0? mapq : 0;
			if (r->score > subsc && mapq == 0) mapq = 1;
			r->mapq = mapq < 60? mapq : 60;
		} else r->mapq = 0;
	}
}
