#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "kalloc.h"
#include "mgpriv.h"

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

static void mg_sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	char buf[16]; // for integer to string conversion
	const char *p, *q;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) str_copy(s, q, p);
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				str_copy(s, r, r + strlen(r));
			} else if (*p == 'c') {
				str_enlarge(s, 1);
				s->s[s->l++] = va_arg(ap, int);
			} else abort();
			q = p + 1;
		}
	}
	if (p > q) str_copy(s, q, p);
	va_end(ap);
	s->s[s->l] = 0;
}

void mg_print_lchain(FILE *fp, const mg_idx_t *gi, int n_lc, const mg_lchain_t *lc, const mg128_t *a, const char *qname)
{
	kstring_t str = {0,0,0};
	int i, j;
	for (i = 0; i < n_lc; ++i) {
		const mg_lchain_t *p = &lc[i];
		str.l = 0;
		mg_sprintf_lite(&str, "LC\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t", qname, p->qs, p->qe, "+-"[p->v&1], gi->g->seg[p->v>>1].name, gi->g->seg[p->v>>1].len, p->rs, p->re, p->cnt);
		for (j = 0; j < p->cnt; ++j)
			mg_sprintf_lite(&str, "%d,", (int32_t)a[p->off + j].y);
		mg_sprintf_lite(&str, "\t");
		for (j = 0; j < p->cnt; ++j)
			mg_sprintf_lite(&str, "%d,", (int32_t)a[p->off + j].x);
		mg_sprintf_lite(&str, "\n");
		fwrite(str.s, 1, str.l, fp);
	}
	free(str.s);
}

void mg_write_paf(kstring_t *s, const gfa_t *g, const mg_gchains_t *gs, int32_t qlen, const char *qname, uint64_t flag, void *km)
{
	int32_t i, j;
	s->l = 0;
	if ((gs == 0 || gs->n_gc == 0) && (flag&MG_M_SHOW_UNMAP)) {
		mg_sprintf_lite(s, "%s\t%d\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\n", qname, qlen);
		return;
	}
	if (gs == 0) return;
	for (i = 0; i < gs->n_gc; ++i) {
		const mg_gchain_t *p = &gs->gc[i];
		int32_t sign_pos, compact;
		if (p->id != p->parent && !(flag&MG_M_PRINT_2ND)) continue;
		if (p->cnt == 0) continue;
		mg_sprintf_lite(s, "%s\t%d\t%d\t%d\t+\t", qname, qlen, p->qs, p->qe);
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
				if (t->pnid < 0) { // no stable ID; write the vertex coordinate
					compact = 0;
					if (last_pnid >= 0) mg_sprintf_lite(s, "%c%s:%d-%d", "><"[rev], g->pseq[last_pnid].name, st, en);
					last_pnid = -1, st = -1, en = -1, rev = -1;
					mg_sprintf_lite(s, "%c%s", "><"[q->v&1], g->seg[q->v>>1].name);
				} else {
					int cont = 0;
					if (last_pnid >= 0 && t->pnid == last_pnid && (q->v&1) == rev) { // same stable sequence and same strand
						if (!(q->v&1)) { // forward strand
							if (t->ppos == en)
								en = t->ppos + t->len, cont = 1;
						} else { // reverse strand
							if (t->ppos + t->len == st)
								st = t->ppos, cont = 1;
						}
					}
					if (cont == 0) {
						if (last_pnid >= 0) compact = 0;
						if (last_pnid >= 0) mg_sprintf_lite(s, "%c%s:%d-%d", "><"[rev], g->pseq[last_pnid].name, st, en);
						last_pnid = t->pnid, rev = q->v&1, st = t->ppos, en = st + t->len;
					}
				}
			}
			if (last_pnid >= 0) {
				if (g->pseq[last_pnid].rank != 0 || g->pseq[last_pnid].min != 0)
					compact = 0;
				if (!compact) mg_sprintf_lite(s, "%c%s:%d-%d", "><"[rev], g->pseq[last_pnid].name, st, en);
			} else compact = 0;
		}
		if (compact) {
			const gfa_seg_t *t = &g->seg[gs->lc[p->off].v>>1];
			const gfa_pseq_t *ps = &g->pseq[t->pnid];
			if (gs->lc[p->off].v&1) s->s[sign_pos] = '-';
			mg_sprintf_lite(s, "%s\t%d\t%d\t%d", ps->name, ps->max, t->ppos + p->ps, t->ppos + p->pe);
		} else mg_sprintf_lite(s, "\t%d\t%d\t%d", p->plen, p->ps, p->pe);
		mg_sprintf_lite(s, "\t%d\t%d\t%d", p->mlen, p->blen, p->mapq);
		mg_sprintf_lite(s, "\ttp:A:%c\tcm:i:%d\ts1:i:%d\ts2:i:%d", p->id == p->parent? 'P' : 'S', p->n_anchor, p->score, p->subsc);
		if (p->div >= 0.0f && p->div <= 1.0f) {
			char buf[16];
			if (p->div == 0.0f) buf[0] = '0', buf[1] = 0;
			else snprintf(buf, 16, "%.4f", p->div);
			mg_sprintf_lite(s, "\tdv:f:%s", buf);
		}
		mg_sprintf_lite(s, "\n");
	}
}
