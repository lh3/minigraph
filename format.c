#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "kalloc.h"
#include "mgpriv.h"

static inline void str_enlarge(void *km, kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup32(s->m);
		s->s = Krealloc(km, char, s->s, s->m);
	}
}

static inline void str_copy(void *km, kstring_t *s, const char *st, const char *en)
{
	str_enlarge(km, s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

void mg_str_write(void *km, kstring_t *s, int32_t len, char *str)
{
	str_copy(km, s, str, str + len);
}

void mg_str_reserve(void *km, kstring_t *s, int32_t len)
{
	str_enlarge(km, s, len);
}

void mg_sprintf_core(void *km, kstring_t *s, const char *fmt, va_list ap)
{
	char buf[16]; // for integer to string conversion
	const char *p, *q;
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) str_copy(km, s, q, p);
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(km, s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				str_enlarge(km, s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				str_copy(km, s, r, r + strlen(r));
			} else if (*p == 'c') {
				str_enlarge(km, s, 1);
				s->s[s->l++] = va_arg(ap, int);
			} else abort();
			q = p + 1;
		}
	}
	if (p > q) str_copy(km, s, q, p);
	s->s[s->l] = 0;
}

void mg_sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	mg_sprintf_core(0, s, fmt, ap);
	va_end(ap);
}

void mg_sprintf_km(void *km, kstring_t *s, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	mg_sprintf_core(km, s, fmt, ap);
	va_end(ap);
}

void mg_print_lchain(FILE *fp, const mg_idx_t *gi, int n_lc, const mg_lchain_t *lc, const mg128_t *a, const char *qname)
{
	kstring_t str = {0,0,0};
	int i, j;
	for (i = 0; i < n_lc; ++i) {
		const mg_lchain_t *p = &lc[i];
		int mlen, blen, span = a[p->off].y>>32&0xff;
		mlen = blen = span;
		for (j = 1; j < p->cnt; ++j) {
			int ql = (int32_t)a[p->off + j].y - (int32_t)a[p->off + j - 1].y;
			int pl = (int32_t)a[p->off + j].x - (int32_t)a[p->off + j - 1].x;
			blen += pl > ql? pl : ql;
			mlen += pl > span && ql > span? span : pl < ql? pl : ql;
		}
		str.l = 0;
		mg_sprintf_lite(&str, "LC\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", qname, p->qs, p->qe, "+-"[p->v&1], gi->g->seg[p->v>>1].name, gi->g->seg[p->v>>1].len,
						p->rs, p->re, p->score, mlen, blen, p->cnt);
		for (j = 0; j < p->cnt; ++j)
			mg_sprintf_lite(&str, "%d,", (int32_t)a[p->off + j].y);
		mg_sprintf_lite(&str, "\t");
		for (j = 0; j < p->cnt; ++j)
			mg_sprintf_lite(&str, "%d,", (int32_t)a[p->off + j].x);
		mg_sprintf_lite(&str, "\t");
		for (j = 0; j < p->cnt; ++j)
			mg_sprintf_lite(&str, "%d,", (int32_t)(a[p->off + j].y>>MG_SEED_OCC_SHIFT));
		mg_sprintf_lite(&str, "\n");
		fwrite(str.s, 1, str.l, fp);
	}
	free(str.s);
}

void mg_write_gaf(kstring_t *s, const gfa_t *g, const mg_gchains_t *gs, int32_t n_seg, const int32_t *qlens, const char *qname, uint64_t flag, void *km)
{
	int32_t i, j, qlen, rev_sign = 0;
	s->l = 0;
	for (i = 0, qlen = 0; i < n_seg; ++i) qlen += qlens[i];
	if ((gs == 0 || gs->n_gc == 0) && (flag&MG_M_SHOW_UNMAP)) {
		mg_sprintf_lite(s, "%s", qname);
		if ((flag&MG_M_FRAG_MERGE) && n_seg == 2 && s->l > 2 && s->s[s->l-1] == '1' && s->s[s->l-2] == '/') s->l -= 2;
		mg_sprintf_lite(s, "\t%d\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\n", qlen);
		return;
	}
	if (gs == 0) return;
	for (i = 0; i < gs->n_gc; ++i) {
		const mg_gchain_t *p = &gs->gc[i];
		int32_t sign_pos, compact;
		if (p->id != p->parent && !(flag&MG_M_PRINT_2ND)) continue;
		if (p->cnt == 0) continue;
		mg_sprintf_lite(s, "%s", qname);
		if ((flag&MG_M_FRAG_MERGE) && n_seg == 2 && s->l > 2 && s->s[s->l-1] == '1' && s->s[s->l-2] == '/') s->l -= 2;
		mg_sprintf_lite(s, "\t%d\t%d\t%d\t+\t", qlen, p->qs, p->qe);
		assert(p->cnt > 0);
		sign_pos = s->l - 2;
		if (flag & MG_M_VERTEX_COOR) {
			compact = 0;
			for (j = 0; j < p->cnt; ++j) {
				const mg_llchain_t *q = &gs->lc[p->off + j];
				mg_sprintf_lite(s, "%c%s", "><"[q->v&1], g->seg[q->v>>1].name);
			}
		} else {
			int32_t last_pnid = -1, st = -1, en = -1, rev = -1;
			compact = flag&MG_M_NO_COMP_PATH? 0 : 1;
			for (j = 0; j < p->cnt; ++j) {
				const mg_llchain_t *q;
				const gfa_seg_t *t;
				assert(p->off + j < gs->n_lc);
				q = &gs->lc[p->off + j];
				t = &g->seg[q->v>>1];
				if (t->snid < 0) { // no stable ID; write the vertex coordinate
					compact = 0;
					if (last_pnid >= 0) mg_sprintf_lite(s, "%c%s:%d-%d", "><"[rev], g->sseq[last_pnid].name, st, en);
					last_pnid = -1, st = -1, en = -1, rev = -1;
					mg_sprintf_lite(s, "%c%s", "><"[q->v&1], g->seg[q->v>>1].name);
				} else {
					int cont = 0;
					if (last_pnid >= 0 && t->snid == last_pnid && (q->v&1) == rev) { // same stable sequence and same strand
						if (!(q->v&1)) { // forward strand
							if (t->soff == en)
								en = t->soff + t->len, cont = 1;
						} else { // reverse strand
							if (t->soff + t->len == st)
								st = t->soff, cont = 1;
						}
					}
					if (cont == 0) {
						if (last_pnid >= 0) compact = 0;
						if (last_pnid >= 0) mg_sprintf_lite(s, "%c%s:%d-%d", "><"[rev], g->sseq[last_pnid].name, st, en);
						last_pnid = t->snid, rev = q->v&1, st = t->soff, en = st + t->len;
					}
				}
			}
			if (last_pnid >= 0) {
				if (g->sseq[last_pnid].rank != 0 || g->sseq[last_pnid].min != 0)
					compact = 0;
				if (!compact) mg_sprintf_lite(s, "%c%s:%d-%d", "><"[rev], g->sseq[last_pnid].name, st, en);
			} else compact = 0;
		}
		if (compact) {
			int32_t rev = gs->lc[p->off].v&1;
			const gfa_seg_t *t = &g->seg[gs->lc[rev? p->off + p->cnt - 1 : p->off].v>>1];
			const gfa_sseq_t *ps = &g->sseq[t->snid];
			mg_sprintf_lite(s, "%s\t%d\t", ps->name, ps->max);
			if (rev) {
				rev_sign = 1;
				s->s[sign_pos] = '-';
				mg_sprintf_lite(s, "%d\t%d", t->soff + (p->plen - p->pe), t->soff + (p->plen - p->ps));
			} else {
				mg_sprintf_lite(s, "%d\t%d", t->soff + p->ps, t->soff + p->pe);
			}
		} else mg_sprintf_lite(s, "\t%d\t%d\t%d", p->plen, p->ps, p->pe);
		if (p->p) mg_sprintf_lite(s, "\t%d\t%d\t%d", p->p->mlen, p->p->blen, p->mapq);
		else mg_sprintf_lite(s, "\t%d\t%d\t%d", p->mlen, p->blen, p->mapq);
		mg_sprintf_lite(s, "\ttp:A:%c", p->id == p->parent? 'P' : 'S');
		if (p->p) mg_sprintf_lite(s, "\tNM:i:%d", p->p->blen - p->p->mlen);
		mg_sprintf_lite(s, "\tcm:i:%d\ts1:i:%d\ts2:i:%d", p->n_anchor, p->score, p->subsc);
		if (p->div >= 0.0f && p->div <= 1.0f) {
			char buf[16];
			if (p->div == 0.0f) buf[0] = '0', buf[1] = 0;
			else snprintf(buf, 16, "%.4f", p->div);
			mg_sprintf_lite(s, "\tdv:f:%s", buf);
		}
		if (n_seg > 1) {
			mg_sprintf_lite(s, "\tql:B:i");
			for (j = 0; j < n_seg; ++j) mg_sprintf_lite(s, ",%d", qlens[j]);
		}
		if (p->p) {
			mg_sprintf_lite(s, "\tcg:Z:");
			if (rev_sign)
				for (j = p->p->n_cigar - 1; j >= 0; --j)
					mg_sprintf_lite(s, "%d%c", (int32_t)(p->p->cigar[j]>>4), "MIDNSHP=XB"[p->p->cigar[j]&0xf]);
			else
				for (j = 0; j < p->p->n_cigar; ++j)
					mg_sprintf_lite(s, "%d%c", (int32_t)(p->p->cigar[j]>>4), "MIDNSHP=XB"[p->p->cigar[j]&0xf]);
		}
		if (p->ds.ds) {
			mg_sprintf_lite(s, "\tds:Z:");
			if (rev_sign) {
				const char *ds = p->ds.ds;
				int32_t i, j;
				for (i = p->ds.n_off - 1; i >= 0; --i) {
					int32_t off = p->ds.off[i], en;
					mg_sprintf_lite(s, "%c", ds[off]); // print the operator
					en = i < p->ds.n_off - 1? p->ds.off[i+1] : p->ds.len;
					if (ds[off] == ':') {
						for (j = off + 1; j < en; ++j)
							mg_sprintf_lite(s, "%c", ds[j]);
					} else if (ds[off] == '*') {
						for (j = off + 1; j < en; ++j)
							mg_sprintf_lite(s, "%c", gfa_comp_table[(uint8_t)ds[j]]);
					} else {
						for (j = en - 1; j >= off + 1; --j) {
							if (ds[j] == '[') mg_sprintf_lite(s, "]");
							else if (ds[j] == ']') mg_sprintf_lite(s, "[");
							else mg_sprintf_lite(s, "%c", gfa_comp_table[(uint8_t)ds[j]]);
						}
					}
				}
			} else {
				mg_sprintf_lite(s, "%s", p->ds.ds);
			}
		}
		mg_sprintf_lite(s, "\n");
		if ((mg_dbg_flag & MG_DBG_LCHAIN) || (flag & MG_M_WRITE_LCHAIN)) {
			char buf[16];
			for (j = 0; j < p->cnt; ++j) {
				const mg_llchain_t *lc = &gs->lc[p->off + j];
				mg_sprintf_lite(s, "*\t%c%s\t%d\t%d", "><"[lc->v&1], g->seg[lc->v>>1].name, g->seg[lc->v>>1].len, lc->cnt);
				if (lc->cnt > 0) {
					double div;
					int32_t q_span = (int32_t)(gs->a[lc->off].y>>32&0xff);
					int32_t n = (int32_t)(gs->a[lc->off + lc->cnt - 1].x>>32) - (int32_t)(gs->a[lc->off].x>>32) + 1;
					div = n == lc->cnt? 0.0 : (n > lc->cnt? log((double)n / lc->cnt) : log((double)lc->cnt / n)) / q_span;
					if (div == 0.0) buf[0] = '0', buf[1] = 0;
					else snprintf(buf, 16, "%.4f", div);
					mg_sprintf_lite(s, "\t%s", buf);
					mg_sprintf_lite(s, "\t%d\t%d", (int32_t)gs->a[lc->off].x + 1 - q_span, (int32_t)gs->a[lc->off + lc->cnt - 1].x + 1);
					mg_sprintf_lite(s, "\t%d\t%d", (int32_t)gs->a[lc->off].y + 1 - q_span, (int32_t)gs->a[lc->off + lc->cnt - 1].y + 1);
					if (flag & MG_M_WRITE_MZ) {
						int32_t i, last;
						last = (int32_t)gs->a[lc->off].x + 1 - q_span;
						mg_sprintf_lite(s, "\t%d\t", q_span);
						for (i = 1; i < lc->cnt; ++i) {
							int32_t x = (int32_t)gs->a[lc->off + i].x + 1 - q_span;
							if (i > 1) mg_sprintf_lite(s, ",");
							mg_sprintf_lite(s, "%d", x - last);
							last = x;
						}
						last = (int32_t)gs->a[lc->off].y + 1 - q_span;
						mg_sprintf_lite(s, "\t");
						for (i = 1; i < lc->cnt; ++i) {
							int32_t x = (int32_t)gs->a[lc->off + i].y + 1 - q_span;
							if (i > 1) mg_sprintf_lite(s, ",");
							mg_sprintf_lite(s, "%d", x - last);
							last = x;
						}
					}
				}
				mg_sprintf_lite(s, "\n");
			}
		}
	}
}
