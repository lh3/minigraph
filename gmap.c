#include <stdlib.h>
#include <assert.h>
#include "kthread.h"
#include "kalloc.h"
#include "bseq.h"
#include "sys.h"
#include "mgpriv.h"
#include "gfa-priv.h"

typedef struct {
	int64_t mini_batch_size;
	int n_processed, n_threads, n_fp;
	const mg_mapopt_t *opt;
	mg_bseq_file_t **fp;
	const mg_idx_t *gi;
	kstring_t str;
	double *c_seg, *c_link;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int n_seq, n_frag;
	mg_bseq1_t *seq;
	int *seg_off, *n_seg;
	mg_gchains_t **gcs;
	mg_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *s = (step_t*)_data;
	int qlens[MG_MAX_SEG], j, off = s->seg_off[i], pe_ori = s->p->opt->pe_ori;
	const char *qseqs[MG_MAX_SEG];
	mg_tbuf_t *b = s->buf[tid];
	assert(s->n_seg[i] <= MG_MAX_SEG);
	if (mg_dbg_flag & MG_DBG_QNAME)
		fprintf(stderr, "QR\t%s\t%d\t%d\n", s->seq[off].name, tid, s->seq[off].l_seq);
	for (j = 0; j < s->n_seg[i]; ++j) {
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1))))
			mg_revcomp_bseq(&s->seq[off + j]);
		qlens[j] = s->seq[off + j].l_seq;
		qseqs[j] = s->seq[off + j].seq;
	}
	if (s->p->opt->flag & MG_M_INDEPEND_SEG) {
		for (j = 0; j < s->n_seg[i]; ++j)
			mg_map_frag(s->p->gi, 1, &qlens[j], &qseqs[j], &s->gcs[off+j], b, s->p->opt, s->seq[off+j].name);
	} else {
		mg_map_frag(s->p->gi, s->n_seg[i], qlens, qseqs, &s->gcs[off], b, s->p->opt, s->seq[off].name);
	}
#if 0 // for paired-end reads
	for (j = 0; j < s->n_seg[i]; ++j) // flip the query strand and coordinate to the original read strand
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1)))) {
			int k, t;
			mg_revcomp_bseq(&s->seq[off + j]);
			for (k = 0; k < s->n_reg[off + j]; ++k) {
				mg_lchain_t *r = &s->reg[off + j][k];
				t = r->qs;
				r->qs = qlens[j] - r->qe;
				r->qe = qlens[j] - t;
				r->v ^= 1;
			}
		}
#endif
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i, j, k;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
		int with_qual = !(p->opt->flag & MG_M_NO_QUAL);
		int with_comment = !!(p->opt->flag & MG_M_COPY_COMMENT);
		int frag_mode = (p->n_fp > 1 || !!(p->opt->flag & MG_M_FRAG_MODE));
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
		if (p->n_fp > 1) s->seq = mg_bseq_read_frag(p->n_fp, p->fp, p->mini_batch_size, with_qual, with_comment, &s->n_seq);
		else s->seq = mg_bseq_read(p->fp[0], p->mini_batch_size, with_qual, with_comment, frag_mode, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i)
				mg_toupper(s->seq[i].l_seq, s->seq[i].seq);
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].rid = p->n_processed++;
			s->buf = (mg_tbuf_t**)calloc(p->n_threads, sizeof(mg_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = mg_tbuf_init();
			s->seg_off = (int*)calloc(2 * s->n_seq, sizeof(int));
			s->n_seg = s->seg_off + s->n_seq; // n_seg, rep_len and frag_gap are allocated together with seg_off
			KCALLOC(0, s->gcs, s->n_seq);
			for (i = 1, j = 0; i <= s->n_seq; ++i)
				if (i == s->n_seq || !frag_mode || !mg_qname_same(s->seq[i-1].name, s->seq[i].name)) {
					s->n_seg[s->n_frag] = i - j;
					s->seg_off[s->n_frag++] = j;
					j = i;
				}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_frag);
		return in;
    } else if (step == 2) { // step 2: output
		void *km = 0;
        step_t *s = (step_t*)in;
		for (i = 0; i < p->n_threads; ++i) mg_tbuf_destroy(s->buf[i]);
		free(s->buf);
		if (!(mg_dbg_flag & MG_DBG_NO_KALLOC)) km = km_init();
		for (k = 0; k < s->n_frag; ++k) {
			int seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
			if ((p->opt->flag & MG_M_FRAG_MODE) && (p->opt->flag & MG_M_FRAG_MERGE)) {
				mg_bseq1_t *t = &s->seq[seg_st];
				int32_t *qlens;
				KMALLOC(km, qlens, seg_en - seg_st); // TODO: if this is an issue (quite unlikely), preallocate
				for (i = seg_st; i < seg_en; ++i)
					qlens[i - seg_st] = s->seq[i].l_seq;
				if (p->opt->flag & MG_M_CAL_COV)
					mg_cov_map(p->gi->g, s->gcs[seg_st], p->opt->min_cov_mapq, p->opt->min_cov_blen, p->c_seg, p->c_link, t->name);
				else mg_write_gaf(&p->str, p->gi->g, s->gcs[seg_st], seg_en - seg_st, qlens, t->name, p->opt->flag, km);
				kfree(km, qlens);
				if (p->str.l) mg_err_fputs(p->str.s, stdout);
			} else {
				for (i = seg_st; i < seg_en; ++i) {
					mg_bseq1_t *t = &s->seq[i];
					if (p->opt->flag & MG_M_CAL_COV)
						mg_cov_map(p->gi->g, s->gcs[i], p->opt->min_cov_mapq, p->opt->min_cov_blen, p->c_seg, p->c_link, t->name);
					else mg_write_gaf(&p->str, p->gi->g, s->gcs[i], 1, &t->l_seq, t->name, p->opt->flag, km);
					if (p->str.l) mg_err_fputs(p->str.s, stdout);
				}
			}
			for (i = seg_st; i < seg_en; ++i) {
				mg_gchain_free(s->gcs[i]);
				free(s->seq[i].seq); free(s->seq[i].name);
				if (s->seq[i].qual) free(s->seq[i].qual);
				if (s->seq[i].comment) free(s->seq[i].comment);
			}
		}
		free(s->gcs); free(s->seg_off); free(s->seq); // n_seg, rep_len and frag_gap were allocated with seg_off; no memory leak here
		if (km) km_destroy(km);
		if (mg_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), s->n_seq);
		free(s);
	}
    return 0;
}

static mg_bseq_file_t **open_bseqs(int n, const char **fn)
{
	mg_bseq_file_t **fp;
	int i, j;
	fp = (mg_bseq_file_t**)calloc(n, sizeof(mg_bseq_file_t*));
	for (i = 0; i < n; ++i) {
		if ((fp[i] = mg_bseq_open(fn[i])) == 0) {
			if (mg_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s'\n", fn[i]);
			for (j = 0; j < i; ++j)
				mg_bseq_close(fp[j]);
			free(fp);
			return 0;
		}
	}
	return fp;
}

int mg_map_file_frag(const mg_idx_t *idx, int n_segs, const char **fn, const mg_mapopt_t *opt, int n_threads, double *c_seg, double *c_link)
{
	int i, pl_threads;
	pipeline_t pl;
	if (n_segs < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp = open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.gi = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	pl.mini_batch_size = opt->mini_batch_size;
	pl.c_seg = c_seg, pl.c_link = c_link;
	pl_threads = n_threads == 1? 1 : (opt->flag&MG_M_2_IO_THREADS)? 3 : 2;
	kt_pipeline(pl_threads, worker_pipeline, &pl, 3);

	free(pl.str.s);
	for (i = 0; i < pl.n_fp; ++i)
		mg_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 0;
}

int mg_map_files(gfa_t *g, int n_fn, const char **fn, const mg_idxopt_t *ipt, const mg_mapopt_t *opt0, int n_threads)
{
	mg_mapopt_t opt = *opt0;
	mg_idx_t *gi;
	int i, ret = 0;
	double *cov_seg = 0, *cov_link = 0;
	if ((gi = mg_index(g, ipt, n_threads, &opt)) == 0) return -1;
	if (opt.flag & MG_M_CAL_COV) {
		KCALLOC(0, cov_seg,  g->n_seg);
		KCALLOC(0, cov_link, g->n_arc);
	}
	if (opt.flag & MG_M_FRAG_MODE) {
		ret = mg_map_file_frag(gi, n_fn, fn, &opt, n_threads, cov_seg, cov_link);
	} else {
		for (i = 0; i < n_fn; ++i) {
			ret = mg_map_file_frag(gi, 1, &fn[i], &opt, n_threads, cov_seg, cov_link);
			if (ret != 0) break;
		}
	}
	if (opt.flag & MG_M_CAL_COV) {
		gfa_aux_update_cv(g, "dc", cov_seg, cov_link);
		free(cov_seg); free(cov_link);
	}
	mg_idx_destroy(gi);
	return ret;
}
