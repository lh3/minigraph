#include <assert.h>
#include "mgpriv.h"
#include "bseq.h"
#include "mss.h"

/**********************
 * Graph augmentation *
 **********************/

void mg_ggsimple(void *km, const mg_ggopt_t *opt, gfa_t *g, int32_t n_seq, const mg_bseq1_t *seq, mg_gchains_t *const* gcs)
{
	int32_t t, i, j, *scnt, *soff, *qcnt, *qoff, max_acnt, *sc, m_ovlp = 0, *ovlp = 0, n_ins, m_ins;
	int32_t l_pseq, m_pseq;
	int64_t sum_acnt, sum_alen;
	uint64_t *meta;
	mg_intv_t *sintv, *qintv;
	double a_dens;
	gfa_ins_t *ins;
	char *pseq;

	// count the number of intervals on each segment
	KCALLOC(km, scnt, g->n_seg);
	KCALLOC(km, qcnt, n_seq);
	for (t = 0, max_acnt = 0; t < n_seq; ++t) {
		const mg_gchains_t *gt = gcs[t];
		for (i = 0; i < gt->n_gc; ++i) {
			const mg_gchain_t *gc = &gt->gc[i];
			if (gc->id != gc->parent) continue;
			if (gc->blen < opt->min_depth_len || gc->mapq < opt->min_mapq) continue;
			if (gc->n_anchor > max_acnt) max_acnt = gc->n_anchor;
			++qcnt[t];
			for (j = 0; j < gc->cnt; ++j)
				++scnt[gt->lc[gc->off + j].v>>1];
		}
	}
	if (max_acnt == 0) { // no gchain
		kfree(km, scnt); kfree(km, qcnt);
		return;
	}

	// compute soff[] and qoff[]
	KMALLOC(km, soff, g->n_seg + 1);
	KMALLOC(km, qoff, n_seq + 1);
	for (soff[0] = 0, i = 1; i <= g->n_seg; ++i)
		soff[i] = soff[i - 1] + scnt[i - 1];
	for (qoff[0] = 0, i = 1; i <= n_seq; ++i)
		qoff[i] = qoff[i - 1] + qcnt[i - 1];

	// populate the interval list
	memset(scnt, 0, 4 * g->n_seg);
	memset(qcnt, 0, 4 * n_seq);
	KMALLOC(km, sintv, soff[g->n_seg]);
	KMALLOC(km, qintv, qoff[n_seq]);
	sum_acnt = sum_alen = 0;
	for (t = 0; t < n_seq; ++t) {
		const mg_gchains_t *gt = gcs[t];
		for (i = 0; i < gt->n_gc; ++i) {
			const mg_gchain_t *gc = &gt->gc[i];
			mg_intv_t *p;
			if (gc->id != gc->parent) continue;
			if (gc->blen < opt->min_depth_len || gc->mapq < opt->min_mapq) continue;
			p = &qintv[qoff[t] + qcnt[t]];
			++qcnt[t];
			p->st = gc->qs, p->en = gc->qe, p->rev = 0, p->far = -1, p->i = -1;
			for (j = 0; j < gc->cnt; ++j) {
				const mg_llchain_t *lc = &gt->lc[gc->off + j];
				int32_t rs, re, tmp;
				if (lc->cnt > 0) { // compute start and end on the forward strand on the segment
					const mg128_t *q = &gt->a[lc->off];
					rs = (int32_t)q->x + 1 - (int32_t)(q->y>>32&0xff);
					q = &gt->a[lc->off + lc->cnt - 1];
					re = (int32_t)q->x;
					assert(rs >= 0 && re > rs && re < g->seg[lc->v>>1].len);
					sum_alen += re - rs, sum_acnt += (q->x>>32) - (gt->a[lc->off].x>>32) + 1;
					if (lc->v&1)
						tmp = rs, rs = g->seg[lc->v>>1].len - re, re = g->seg[lc->v>>1].len - tmp;
				} else rs = 0, re = g->seg[lc->v>>1].len;
				// save the interval
				p = &sintv[soff[lc->v>>1] + scnt[lc->v>>1]];
				++scnt[lc->v>>1];
				p->st = rs, p->en = re, p->rev = lc->v&1, p->far = -1, p->i = -1;
			}
		}
	}
	a_dens = (double)sum_acnt / sum_alen;

	// sort and index intervals
	for (i = 0; i < g->n_seg; ++i) {
		assert(soff[i+1] - soff[i] == scnt[i]);
		mg_intv_index(soff[i+1] - soff[i], &sintv[soff[i]]);
	}
	kfree(km, scnt);
	for (i = 0; i < n_seq; ++i) {
		assert(qoff[i+1] - qoff[i] == qcnt[i]);
		mg_intv_index(qoff[i+1] - qoff[i], &qintv[qoff[i]]);
	}
	kfree(km, qcnt);

	// extract poorly regions
	m_pseq = l_pseq = 0, pseq = 0;
	m_ins = n_ins = 0, ins = 0;
	KMALLOC(km, sc, max_acnt);
	KMALLOC(km, meta, max_acnt);
	for (t = 0; t < n_seq; ++t) {
		const mg_gchains_t *gt = gcs[t];
		for (i = 0; i < gt->n_gc; ++i) {
			const mg_gchain_t *gc = &gt->gc[i];
			int32_t off_a, off_l, n_ss, far_q;
			msseg_t *ss;
			if (gc->id != gc->parent) continue;
			if (gc->blen < opt->min_map_len || gc->mapq < opt->min_mapq) continue;
			assert(gc->cnt > 0);

			// fill sc[]. This part achieves a similar goal to the one in mg_gchain_extra(). It makes more assumptions, but is logically simpler.
			off_l = gc->off;
			off_a = gt->lc[off_l].off + 1;
			far_q = 0;
			for (j = 1; j < gc->n_anchor; ++j, ++off_a) {
				const mg128_t *q = &gt->a[off_a - 1], *p = &gt->a[off_a];
				const mg_llchain_t *lc = &gt->lc[off_l];
				int32_t s, off_l0 = off_l, pd, qd = (int32_t)p->y - (int32_t)q->y, c = (int32_t)(p->x>>32) - (int32_t)(q->x>>32) - 1;
				if ((int32_t)q->y > far_q) far_q = (int32_t)q->y; // far_q keeps the rightmost query position seen so far
				if (off_a == lc->off + lc->cnt) { // we are at the end of the current lchain
					pd = g->seg[lc->v>>1].len - (int32_t)q->x - 1;
					for (++off_l; off_l < gc->off + gc->cnt && gt->lc[off_l].cnt == 0; ++off_l)
						pd += g->seg[gt->lc[off_l].v>>1].len;
					assert(off_l < gc->off + gc->cnt);
					pd += (int32_t)p->x + 1;
				} else pd = (int32_t)p->x - (int32_t)q->x;
				if (pd == qd && c == 0) s = -opt->match_pen;
				else if ((int32_t)p->y < far_q) s = 1; // query overlap
				else if (pd > qd) s = (int32_t)(c + (pd - qd) * a_dens + .499);
				else s = c;
				sc[j - 1] = s;
				meta[j-1] = (uint64_t)pd<<32 | off_l0;
			}

			// get regions to insert
			ss = mss_find_all(0, gc->n_anchor - 1, sc, 10, 0, &n_ss);
			off_a = gt->lc[gc->off].off;
			for (j = 0; j < n_ss; ++j) {
				const mg128_t *p, *q;
				int32_t st, en, ls, le, span, pd, k, n_ovlp;
				gfa_ins_t I;

				// find the initial positions
				if (ss[j].st <= 1 || ss[j].en >= gc->n_anchor - 1) continue; // not at the ends
				st = ss[j].st - 1, en = ss[j].en;
				q = &gt->a[off_a + st];
				p = &gt->a[off_a + en];
				span = p->y>>32&0xff;
				I.ctg = t;
				ls = (int32_t)meta[st], le = (int32_t)meta[en]; // first and last lchain; CLOSED interval
				assert(ls <= le);
				I.v[0] = gt->lc[ls].v;
				I.v[1] = gt->lc[le].v;
				I.voff[0] = (int32_t)q->x + 1;
				I.voff[1] = (int32_t)p->x + 1 - span;
				I.coff[0] = (int32_t)q->y + 1;
				I.coff[1] = (int32_t)p->y + 1 - span;
				assert(I.voff[0] <= g->seg[I.v[0]>>1].len);
				assert(I.voff[1] <= g->seg[I.v[1]>>1].len);
				for (k = st, pd = 0; k < en; ++k) pd += meta[k]>>32;
				pd -= span;

				// adjust for overlapping poistions
				if (I.v[0] == I.v[1] && I.voff[0] > I.voff[1]) {
					assert(I.voff[0] - I.voff[1] <= span);
					I.coff[1] += I.voff[0] - I.voff[1];
					pd += I.voff[0] - I.voff[1];
					I.voff[1] = I.voff[0];
				}
				if (I.coff[0] > I.coff[1]) {
					int32_t d = I.coff[0] - I.coff[1];
					int32_t l1 = g->seg[I.v[1]>>1].len;
					if (d > span || d > I.voff[0] + (l1 - I.voff[1])) {
						if (mg_verbose >= 2)
							fprintf(stderr, "[W::%s] unexpected insert [%c%s:%d,%c%s:%d|%d] <=> %s:[%d,%d|%d]\n", __func__, "><"[I.v[0]&1], g->seg[I.v[0]>>1].name, I.voff[0], "><"[I.v[1]&1], g->seg[I.v[1]>>1].name, I.voff[1], pd, seq[t].name, I.coff[0], I.coff[1], I.coff[1] - I.coff[0]);
						continue; // such overlap can't be properly resolved
					}
					if (I.voff[1] + d <= l1) {
						I.voff[1] += d, pd += d;
						I.coff[1] = I.coff[0];
					} else {
						int32_t x = l1 - I.voff[1];
						d -= x;
						I.voff[1] += x, I.coff[1] += x, pd += x;
						assert(I.voff[0] >= d);
						I.voff[0] -= d, I.coff[0] -= d, pd += d;
					}
				}
				assert(I.voff[0] <= g->seg[I.v[0]>>1].len);
				assert(I.voff[1] <= g->seg[I.v[1]>>1].len);
				pd -= gfa_ins_adj(g, opt->ggs_shrink_pen, &I, seq[t].seq); // "3" is not used for now

				// filtering
				if (I.coff[1] - I.coff[0] < opt->min_var_len && pd < opt->min_var_len)
					continue;
				for (k = I.coff[0]; k < I.coff[1]; ++k) { // test ambiguous bases
					int c = seq[t].seq[k];
					if (c == 'n' || c == 'N') break;
				}
				if (k != I.coff[1]) continue; // no ambiguous bases on the insert
				n_ovlp = mg_intv_overlap(km, qoff[t+1] - qoff[t], &qintv[qoff[t]], I.coff[0], I.coff[1], &ovlp, &m_ovlp); // test overlapping on the query
				if (n_ovlp == 0) fprintf(stderr, "[W::%s] query interval %s:%d-%d is not covered\n", __func__, seq[t].name, I.coff[0], I.coff[1]);
				if (n_ovlp != 1) continue;
				for (k = ls; k <= le; ++k) { // find other mappings overlapping with the insert on the graph
					uint32_t v = gt->lc[k].v, len = g->seg[v>>1].len;
					int32_t s = 0, e = len, tmp;
					if (k == ls) s = (int32_t)gt->a[off_a+st].x + 1 - (int32_t)(gt->a[off_a+st].y>>32&0xff);
					if (k == le) e = (int32_t)gt->a[off_a+en].x + 1;
					if (v&1) tmp = s, s = len - e, e = len - tmp;
					n_ovlp = mg_intv_overlap(km, soff[(v>>1)+1] - soff[v>>1], &sintv[soff[v>>1]], s, e, &ovlp, &m_ovlp);
					if (n_ovlp == 0) fprintf(stderr, "[W::%s] graph interval %s:%d-%d is not covered\n", __func__, g->seg[v>>1].name, s, e); // this should be an assert()
					if (n_ovlp != 1) break;
				}
				if (k <= le) continue;
				if (pd - (I.coff[1] - I.coff[0]) < opt->min_var_len && (I.coff[1] - I.coff[0]) - pd < opt->min_var_len) {
					int32_t qd = I.coff[1] - I.coff[0], mlen;
					l_pseq = mg_path2seq(km, g, gt, ls, le, I.voff, &pseq, &m_pseq);
					mlen = mg_fastcmp(km, l_pseq, pseq, qd, &seq[t].seq[I.coff[0]], opt->ggs_fc_kmer, opt->ggs_fc_max_occ);
					if (mlen > (qd > pd? qd : pd) * opt->ggs_max_mlen) continue;
				}
				if (mg_dbg_flag & MG_DBG_INSERT) {
					l_pseq = mg_path2seq(km, g, gt, ls, le, I.voff, &pseq, &m_pseq);
					fprintf(stderr, "IN\t[%c%s:%d,%c%s:%d|%d] <=> %s:[%d,%d|%d]\n", "><"[I.v[0]&1], g->seg[I.v[0]>>1].name, I.voff[0], "><"[I.v[1]&1], g->seg[I.v[1]>>1].name, I.voff[1], pd, seq[t].name, I.coff[0], I.coff[1], I.coff[1] - I.coff[0]);
					fprintf(stderr, "IP\t%s\nIQ\t", pseq);
					fwrite(&seq[t].seq[I.coff[0]], 1, I.coff[1] - I.coff[0], stderr);
					fprintf(stderr, "\nIS\t%d==%d\t%d\n", pd, l_pseq, mg_fastcmp(km, l_pseq, pseq, I.coff[1] - I.coff[0], &seq[t].seq[I.coff[0]], opt->ggs_fc_kmer, opt->ggs_fc_max_occ));
				}
				if (n_ins == m_ins) KEXPAND(km, ins, m_ins);
				ins[n_ins++] = I;
			}
			kfree(0, ss);
		}
	}
	kfree(km, pseq);
	kfree(km, ovlp);
	kfree(km, sc);
	kfree(km, meta);
	kfree(km, soff); kfree(km, qoff);
	kfree(km, sintv); kfree(km, qintv);

	if (n_ins > 0) {
		char **names, **seqs;
		KMALLOC(km, names, n_seq);
		KMALLOC(km, seqs, n_seq);
		for (i = 0; i < n_seq; ++i)
			names[i] = seq[i].name, seqs[i] = seq[i].seq;
		n_ins = gfa_ins_filter(g, n_ins, ins);
		gfa_augment(g, n_ins, ins, n_seq, (const char*const*)names, (const char*const*)seqs);
		kfree(km, ins);
		kfree(km, names);
		kfree(km, seqs);
	}
	if (mg_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] inserted %d events\n", __func__,
				realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), n_ins);
}
