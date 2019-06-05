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
			mg_sprintf_lite(&str, "%d,", (int32_t)a[p->as + j].y);
		mg_sprintf_lite(&str, "\t");
		for (j = 0; j < p->cnt; ++j)
			mg_sprintf_lite(&str, "%d,", (int32_t)a[p->as + j].x);
		mg_sprintf_lite(&str, "\n");
		fwrite(str.s, 1, str.l, fp);
	}
	free(str.s);
}

void mg_write_paf(kstring_t *s, const gfa_t *g, const mg_gchains_t *gs, int32_t qlen, const char *qname, void *km)
{
	int32_t i, j;
	s->l = 0;
	for (i = 0; i < gs->n_g; ++i) {
		const mg_gchain_t *p = &gs->g[i];
		mg_sprintf_lite(s, "%s\t%d\t%d\t%d\t+\t", qname, qlen, p->qs, p->qe);
		for (j = 0; j < p->cnt; ++j) {
			const mg_llchain_t *q = &gs->l[p->ls + j];
			mg_sprintf_lite(s, "%c%s", "><"[q->v&1], g->seg[q->v>>1].name);
		}
		mg_sprintf_lite(s, "\t%d\t%d\t%d\n", p->path_len, 0, 0);
	}
}

void mg_print_paf(FILE *fp, const gfa_t *g, const mg_gchains_t *gs, int32_t qlen, const char *qname, void *km)
{
	kstring_t str = {0,0,0};
	mg_write_paf(&str, g, gs, qlen, qname, km);
	fputs(str.s, fp);
	free(str.s);
}
