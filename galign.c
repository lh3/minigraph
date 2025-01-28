#include <assert.h>
#include <string.h>
#include "mgpriv.h"
#include "kalloc.h"
#include "miniwfa.h"

/******************
 * Generate cigar *
 ******************/

static void append_cigar1(void *km, mg64_v *c, int32_t op, int32_t len)
{
	if (c->n > 0 && (c->a[c->n - 1]&0xf) == op) {
		c->a[c->n - 1] += (uint64_t)len<<4;
	} else {
		if (c->n == c->m) {
			c->m += (c->m>>1) + 16;
			KREALLOC(km, c->a, c->m);
		}
		c->a[c->n++] = (uint64_t)len<<4 | op;
	}
}

static void append_cigar(void *km, mg64_v *c, int32_t n_cigar, const uint32_t *cigar)
{
	int32_t k;
	if (n_cigar == 0) return;
	append_cigar1(km, c, cigar[0]&0xf, cigar[0]>>4);
	if (c->n + n_cigar - 1 > c->m) {
		c->m = c->n + n_cigar - 1;
		kroundup32(c->m);
		KREALLOC(km, c->a, c->m);
	}
	for (k = 0; k < n_cigar - 1; ++k)
		c->a[c->n + k] = cigar[1 + k];
	c->n += n_cigar - 1;
}

void mg_gchain_cigar(void *km, const gfa_t *g, const gfa_edseq_t *es, const char *qseq, mg_gchains_t *gt, const char *qname) // qname for debugging only
{
	int32_t i, l_seq = 0, m_seq = 0;
	char *seq = 0;
	void *km2;
	mg64_v cigar = {0,0,0};
	km2 = km_init2(km, 0);
	for (i = 0; i < gt->n_gc; ++i) {
		mg_gchain_t *gc = &gt->gc[i];
		int32_t l0 = gc->off;
		int32_t off_a0 = gt->lc[l0].off;
		int32_t j, j0 = 0, k, l;
		cigar.n = 0;
		append_cigar1(km, &cigar, 7, gt->a[off_a0].y>>32&0xff);
		for (j = 1; j < gc->n_anchor; ++j) {
			const mg128_t *q, *p = &gt->a[off_a0 + j];
			if ((p->y & MG_SEED_IGNORE) && j != gc->n_anchor - 1) continue;
			q = &gt->a[off_a0 + j0];
			// find the lchain that contains the anchor
			for (l = l0; l < gc->off + gc->cnt; ++l) {
				mg_llchain_t *r = &gt->lc[l];
				if (off_a0 + j >= r->off && off_a0 + j < r->off + r->cnt)
					break;
			}
			assert(l < gc->off + gc->cnt);
			assert((int32_t)q->x < g->seg[gt->lc[l0].v>>1].len);
			// calculate the target sequence length
			if (l == l0) {
				l_seq = (int32_t)p->x - (int32_t)q->x;
			} else {
				l_seq = g->seg[gt->lc[l0].v>>1].len - (int32_t)q->x - 1;
				for (k = l0 + 1; k < l; ++k)
					l_seq += es[gt->lc[k].v].len;
				l_seq += (int32_t)p->x + 1;
			}
			if (l_seq + 1 > m_seq) {
				m_seq = l_seq + 1;
				kroundup32(m_seq);
				KREALLOC(km, seq, m_seq);
			}
			// get the target sequence
			if (l == l0) { // on the same vertex
				memcpy(seq, &es[gt->lc[l0].v].seq[(int32_t)q->x + 1], l_seq);
			} else {
				uint32_t v = gt->lc[l0].v;
				l_seq = g->seg[v>>1].len - (int32_t)q->x - 1;
				memcpy(seq, &es[v].seq[(int32_t)q->x + 1], l_seq);
				for (k = l0 + 1; k < l; ++k) {
					v = gt->lc[k].v;
					memcpy(&seq[l_seq], es[v].seq, es[v].len);
					l_seq += es[v].len;
				}
				memcpy(&seq[l_seq], es[gt->lc[l].v].seq, (int32_t)p->x + 1);
				l_seq += (int32_t)p->x + 1;
			}
			{
				int32_t qlen = (int32_t)p->y - (int32_t)q->y;
				const char *qs = &qseq[(int32_t)q->y + 1];
				assert(l_seq > 0 || qlen > 0);
				if (l_seq == 0) append_cigar1(km, &cigar, 1, qlen);
				else if (qlen == 0) append_cigar1(km, &cigar, 2, l_seq);
				else if (l_seq == qlen && qlen <= (q->y>>32&0xff)) append_cigar1(km, &cigar, 7, qlen);
				else {
					mwf_opt_t opt;
					mwf_rst_t rst;
					mwf_opt_init(&opt);
					opt.flag |= MWF_F_CIGAR;
					mwf_wfa_auto(km2, &opt, l_seq, seq, qlen, qs, &rst);
					append_cigar(km, &cigar, rst.n_cigar, rst.cigar);
					kfree(km2, rst.cigar);
					if ((mg_dbg_flag&MG_DBG_MINIWFA) && l_seq > 5000 && qlen > 5000 && rst.s >= 10000)
						fprintf(stderr, "WL\t%s\t%d\t%d\t%d\t%d\t%d\n", qname, i, (int32_t)q->y + 1, (int32_t)p->y - (int32_t)q->y, l_seq, rst.s);
					if (rst.s >= 10000 && l_seq > 5000 && qlen > 5000) {
						km_destroy(km2);
						km2 = km_init2(km, 0);
					}
					if ((mg_dbg_flag&MG_DBG_MWF_SEQ) && l_seq > 5000 && qlen > 5000 && rst.s >= 10000) {
						char *str;
						str = Kmalloc(km, char, qlen + l_seq + strlen(qname) + 100);
						k = sprintf(str, "WL\t%s\t%d\t%d\t%d\nWT\t%.*s\nWQ\t%.*s\n", qname, i, (int32_t)q->y + 1, rst.s, l_seq, seq, qlen, qs);
						fwrite(str, 1, k, stderr);
						kfree(km, str);
					}
				}
			}
			j0 = j, l0 = l;
		}
		// save the CIGAR to gt->gc[i]
		gc->p = (mg_cigar_t*)kcalloc(gt->km, 1, cigar.n * 8 + sizeof(mg_cigar_t));
		gc->p->ss = (int32_t)gt->a[off_a0].x + 1 - (int32_t)(gt->a[off_a0].y>>32&0xff);
		gc->p->ee = (int32_t)gt->a[off_a0 + gc->n_anchor - 1].x + 1;
		gc->p->n_cigar = cigar.n;
		memcpy(gc->p->cigar, cigar.a, cigar.n * 8);
		for (j = 0, l = 0; j < gc->p->n_cigar; ++j) {
			int32_t op = gc->p->cigar[j]&0xf, len = gc->p->cigar[j]>>4;
			if (op == 7) gc->p->mlen += len, gc->p->blen += len;
			else gc->p->blen += len;
			if (op != 1) gc->p->aplen += len;
			if (op != 2) l += len;
		}
		memset(&gc->ds, 0, sizeof(gc->ds));
		assert(l == gc->qe - gc->qs && gc->p->aplen == gc->pe - gc->ps);
	}
	km_destroy(km2);
	kfree(km, seq);
	kfree(km, cigar.a);
}

/***********************
 * Generate the ds tag *
 ***********************/

#define mg_get_nucl(s, i) (seq_nt4_table[(uint8_t)(s)[(i)]]) // get the base in the "nt4" encoding

static void write_indel(void *km, kstring_t *str, int64_t len, const char *seq, int64_t ll, int64_t lr) // write an indel to ds
{
	int64_t i;
	if (ll + lr >= len) {
		mg_sprintf_km(km, str, "[");
		for (i = 0; i < len; ++i)
			mg_sprintf_km(km, str, "%c", "acgtn"[mg_get_nucl(seq, i)]);
		mg_sprintf_km(km, str, "]");
	} else {
		int64_t k = 0;
		if (ll > 0) {
			mg_sprintf_km(km, str, "[");
			for (i = 0; i < ll; ++i)
				mg_sprintf_km(km, str, "%c", "acgtn"[mg_get_nucl(seq, k+i)]);
			mg_sprintf_km(km, str, "]");
			k += ll;
		}
		for (i = 0; i < len - lr - ll; ++i)
			mg_sprintf_km(km, str, "%c", "acgtn"[mg_get_nucl(seq, k+i)]);
		k += len - lr - ll;
		if (lr > 0) {
			mg_sprintf_km(km, str, "[");
			for (i = 0; i < lr; ++i)
				mg_sprintf_km(km, str, "%c", "acgtn"[mg_get_nucl(seq, k+i)]);
			mg_sprintf_km(km, str, "]");
		}
	}
}

void mg_gchain_gen_ds(void *km, const gfa_t *g, const gfa_edseq_t *es, const char *qseq, mg_gchains_t *gt)
{
	int32_t i, m_off = 0, n_off = 0, *off = 0;
	void *km2;
	kstring_t str = {0,0,0}, seq = {0,0,0};
	km2 = km_init2(km, 0);
	for (i = 0; i < gt->n_gc; ++i) {
		mg_gchain_t *gc = &gt->gc[i];
		int32_t j;
		int64_t x, y, ds_len;
		str.l = seq.l = n_off = 0;
		if (gc->p->aplen > seq.m) {
			seq.s = Krealloc(km2, char, seq.s, gc->p->aplen);
			seq.m = gc->p->aplen;
		}
		for (j = 0, seq.l = 0; j < gc->cnt; ++j) { // extract the aligned sequence in the graph
			int32_t k = gc->off + j;
			uint32_t v = gt->lc[k].v;
			int32_t slen = es[v].len;
			int32_t st = j > 0? 0 : gc->p->ss;
			int32_t en = j < gc->cnt - 1? slen : gc->p->ee;
			assert(seq.l + (en - st) <= gc->p->aplen);
			memcpy(&seq.s[seq.l], &es[v].seq[st], en - st);
			seq.l += en - st;
		}
		assert(seq.l == gc->p->aplen);
		for (j = 0, x = 0, y = gc->qs, ds_len = 0, n_off = 0; j < gc->p->n_cigar; ++j) { // estimate the approximate length of ds
			int64_t op = gc->p->cigar[j]&0xf, len = gc->p->cigar[j]>>4, z;
			if (op == 0 || op == 7 || op == 8) { // alignment match
				int32_t l = 0;
				++n_off;
				for (z = 0; z < len; ++z) {
					if (mg_get_nucl(seq.s, x+z) != mg_get_nucl(qseq, y+z))
						ds_len += 3, ds_len += 6, n_off += 2, l = 0;
					else ++l;
				}
				ds_len += 6;
				x += len, y += len;
			} else if (op == 1) { // insertion
				ds_len += len + 1, ++n_off, y += len;
			} else if (op == 2) { // deletion
				ds_len += len + 1, ++n_off, x += len;
			}
		}
		if (n_off > m_off) {
			m_off = n_off + (n_off>>1) + 16;
			off = Krealloc(km2, int32_t, off, m_off);
		}
		mg_str_reserve(km2, &str, ds_len);
		for (j = 0, x = 0, y = gc->qs, n_off = 0; j < gc->p->n_cigar; ++j) { // write ds
			int64_t op = gc->p->cigar[j]&0xf, len = gc->p->cigar[j]>>4;
			assert(n_off < m_off);
			if (op == 0 || op == 7 || op == 8) { // alignment match
				int64_t z;
				int32_t l = 0;
				for (z = 0; z < len; ++z) {
					uint8_t cx = mg_get_nucl(seq.s, x+z);
					uint8_t cy = mg_get_nucl(qseq, y+z);
					if (cx != cy) {
						if (l > 0) {
							off[n_off++] = str.l;
							mg_sprintf_km(km2, &str, ":%d", l);
						}
						off[n_off++] = str.l;
						mg_sprintf_km(km2, &str, "*%c%c", "acgtn"[cx], "acgtn"[cy]);
						l = 0;
					} else ++l;
				}
				if (l > 0) {
					off[n_off++] = str.l;
					mg_sprintf_km(km2, &str, ":%d", l);
				}
				x += len, y += len;
			} else if (op == 1) { // insertion
				int64_t z, ll, lr;
				for (z = 1; z <= len; ++z)
					if (y - z < gc->qs || qseq[y + len - z] != qseq[y - z])
						break;
				lr = z - 1;
				for (z = 0; z < len; ++z)
					if (y + len + z >= gc->qe || qseq[y + len + z] != qseq[y + z])
						break;
				ll = z;
				off[n_off++] = str.l;
				mg_sprintf_km(km2, &str, "+");
				write_indel(km2, &str, len, &qseq[y], ll, lr);
				y += len;
			} else if (op == 2) { // deletion
				int64_t z, ll, lr;
				for (z = 1; z <= len; ++z)
					if (x - z < 0 || seq.s[x + len - z] != seq.s[x - z])
						break;
				lr = z - 1;
				for (z = 0; z < len; ++z)
					if (x + len + z >= gc->p->aplen || seq.s[x + z] != seq.s[x + len + z])
						break;
				ll = z;
				off[n_off++] = str.l;
				mg_sprintf_km(km2, &str, "-");
				write_indel(km2, &str, len, &seq.s[x], ll, lr);
				x += len;
			}
		}
		gc->ds.len = str.l;
		gc->ds.ds = Kcalloc(gt->km, char, str.l + 1);
		memcpy(gc->ds.ds, str.s, str.l);
		gc->ds.n_off = n_off;
		gc->ds.off = Kcalloc(gt->km, int32_t, n_off);
		memcpy(gc->ds.off, off, n_off * sizeof(int32_t));
	}
	km_destroy(km2); // this frees both str.s and seq.s
}
