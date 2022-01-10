#include <string.h>
#include "mgpriv.h"
#include "sys.h"

void mg_idxopt_init(mg_idxopt_t *io)
{
	memset(io, 0, sizeof(mg_idxopt_t));
	io->k = 15;
	io->w = 10;
	io->bucket_bits = 14;
}

void mg_mapopt_init(mg_mapopt_t *mo)
{
	memset(mo, 0, sizeof(mg_mapopt_t));
	mo->seed = 11;
	mo->occ_max1 = 50, mo->occ_max1_cap = 250;
	mo->occ_max1_frac = 2e-4f;
	mo->occ_weight = 20;
	mo->max_gap = 5000;
	mo->max_gap_ref = -1;
	mo->max_gap_pre = 1000;
	mo->max_lc_skip = 25, mo->max_gc_skip = 25;
	mo->max_lc_iter = 5000;
	mo->bw = 500, mo->bw_long = 20000;
	mo->rmq_size_cap = 100000;
	mo->rmq_rescue_size = 1000;
	mo->rmq_rescue_ratio = 0.1f;
	mo->mini_batch_size = 500000000;
	mo->div = 0.1f;
	mo->chn_pen_gap = 1.0f, mo->chn_pen_skip = 0.05f;
	mo->min_lc_cnt = 2, mo->min_lc_score = 30;
	mo->min_gc_cnt = 3, mo->min_gc_score = 50;
	mo->mask_level = 0.5f;
	mo->sub_diff = 6;
	mo->best_n = 5;
	mo->pri_ratio = 0.8f;
	mo->ref_bonus = 0;
	mo->pe_ori = 0; // FF
	mo->max_gc_seq_ext = 5;
	mo->min_cov_mapq = 20;
	mo->min_cov_blen = 1000;
	mo->cap_kalloc = 1000000000;
}

void mg_ggopt_gen_sc_mat(int8_t mat[25], int32_t a, int32_t b)
{
	int i, j, m = 5;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = -1;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = -1;
}

void mg_ggopt_init(mg_ggopt_t *go)
{
	memset(go, 0, sizeof(mg_ggopt_t));
	go->algo = MG_G_NONE;
	go->flag |= MG_G_NO_QOVLP;
	go->min_map_len = 100000;
	go->min_depth_len = 10000;
	go->min_mapq = 5;
	go->min_var_len = 50;
	go->match_pen = 10;
	// for ggs
	go->ggs_shrink_pen = 9;
	go->ggs_min_end_cnt = 10;
	go->ggs_min_end_frac = 0.1f;
	mg_ggopt_gen_sc_mat(go->scmat, 2, -4);
	go->gapo = 4, go->gapo2 = 24;
	go->gape = 2, go->gape2 = 1;
	go->ggs_max_iden = 0.80f;
	go->ggs_min_inv_iden = 0.95f;
}

int mg_opt_set(const char *preset, mg_idxopt_t *io, mg_mapopt_t *mo, mg_ggopt_t *go)
{
	if (preset == 0) {
		mg_idxopt_init(io);
		mg_mapopt_init(mo);
		mg_ggopt_init(go);
	} else if (strcmp(preset, "lr") == 0) { // this is the default
	} else if (strcmp(preset, "asm") == 0 || strcmp(preset, "ggs") == 0) {
		io->k = 19, io->w = 10;
		mo->flag |= MG_M_RMQ;
		mo->bw = 1000, mo->bw_long = 100000;
		mo->max_gap = 10000, mo->max_gap_pre = 1000;
		mo->min_lc_cnt = 3, mo->min_lc_score = 40;
		mo->min_gc_cnt = 5, mo->min_gc_score = 1000;
		mo->min_cov_mapq = 5;
		mo->min_cov_blen = 100000;
		mo->max_lc_skip = mo->max_gc_skip = 50;
		mo->div = 0.01f;
		mo->mini_batch_size = 4000000000LL;
		if (strcmp(preset, "ggs") == 0)
			go->algo = MG_G_GGSIMPLE, mo->best_n = 0;
	} else if (strcmp(preset, "se") == 0 || strcmp(preset, "sr") == 0) {
		io->k = 21, io->w = 10;
		mo->flag |= MG_M_SR | MG_M_HEAP_SORT | MG_M_2_IO_THREADS;
		mo->occ_max1 = 1000;
		mo->occ_max1_cap = 2500;
		mo->max_gap = 100;
		mo->bw = mo->bw_long = 100;
		mo->max_frag_len = 800;
		mo->pri_ratio = 0.5f;
		mo->min_lc_cnt = 2, mo->min_lc_score = 25;
		mo->min_gc_cnt = 3, mo->min_gc_score = 40;
		mo->mini_batch_size = 50000000;
		mo->min_cov_blen = 50;
		mo->chn_pen_gap = 0.2f;
		mo->ref_bonus = 1;
		if (strcmp(preset, "sr") == 0) {
			mo->flag |= MG_M_FRAG_MODE | MG_M_FRAG_MERGE;
			mo->pe_ori = 0<<1|1; // FR
		}
	} else return -1;
	return 0;
}

int mg_opt_check(const mg_idxopt_t *io, const mg_mapopt_t *mo, const mg_ggopt_t *go)
{
	if ((mo->flag & MG_M_FRAG_MODE) && !(mo->flag & MG_M_FRAG_MERGE)) {
		if (mg_verbose >= 1)
			fprintf(stderr, "[ERROR]\033[1;31m the fragment-without-merge mode is not implemented\033[0m\n");
		return -1;
	}
	return 0;
}

void mg_opt_update(const mg_idx_t *gi, mg_mapopt_t *mo, mg_ggopt_t *go)
{
	float f[2];
	int32_t q[2];
	f[0] = 0.05f, f[1] = mo->occ_max1_frac;
	mg_idx_cal_quantile(gi, 2, f, q);
	if (q[0] + 1 > mo->occ_weight) mo->occ_weight = q[0] + 1;
	if (q[1] + 1 > mo->occ_max1)   mo->occ_max1   = q[1] + 1;
	if (mo->occ_max1 > mo->occ_max1_cap) mo->occ_max1 = mo->occ_max1_cap;
	if (mo->bw_long < mo->bw) mo->bw_long = mo->bw;
	if (mg_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] occ_weight=%d, occ_max1=%d; 95 percentile: %d\n", __func__,
				realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), mo->occ_weight, mo->occ_max1, q[0]);
}
