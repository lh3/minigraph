#include <string.h>
#include "mgpriv.h"

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
	mo->mid_occ = 100;
	mo->mid_occ_frac = 2e-4f;
	mo->max_gap = 5000;
	mo->max_gap_ref = -1;
	mo->max_chain_skip = 25;
	mo->bw = 2000;
	mo->mini_batch_size = 500000000;
	mo->min_lc_cnt = 2, mo->min_lc_score = 30;
	mo->min_gc_cnt = 3, mo->min_gc_score = 50;
	mo->mask_level = 0.5f;
	mo->sub_diff = 6;
	mo->best_n = 5;
	mo->pri_ratio = 0.8f;
}

void mg_ggopt_init(mg_ggopt_t *go)
{
	memset(go, 0, sizeof(mg_ggopt_t));
	go->algo = MG_G_NONE;
	go->flag |= MG_G_NO_QOVLP;
	go->min_map_len = 50000;
	go->min_depth_len = 10000;
	go->min_mapq = 5;
	go->min_var_len = 250;
	go->match_pen = 10;
	// for ggs
	go->ggs_shrink_pen = 9;
	go->ggs_fc_kmer = 9, go->ggs_fc_max_occ = 10;
	go->ggs_min_end_cnt = 10;
	go->ggs_min_end_frac = 0.1f;
	go->ggs_max_kiden = 0.8f;
}

int mg_opt_set(const char *preset, mg_idxopt_t *io, mg_mapopt_t *mo, mg_ggopt_t *go)
{
	if (preset == 0) {
		mg_idxopt_init(io);
		mg_mapopt_init(mo);
		mg_ggopt_init(go);
	} else if (strcmp(preset, "lr") == 0) {
		io->k = 15, io->w = 10;
		mo->bw = 2000, mo->max_gap = 5000;
	} else if (strcmp(preset, "asm20") == 0) {
		io->k = 19, io->w = 10;
		mo->bw = 10000, mo->max_gap = 10000;
		mo->min_lc_cnt = 3, mo->min_lc_score = 40;
		mo->min_gc_cnt = 5, mo->min_gc_score = 1000;
	} else if (strcmp(preset, "ggs") == 0 || strcmp(preset, "ggsimple") == 0) {
		io->k = 19, io->w = 10;
		go->algo = MG_G_GGSIMPLE;
		mo->best_n = 0;
		mo->max_gap = mo->bw = 10000;
		mo->min_lc_cnt = 3, mo->min_lc_score = 40;
		mo->min_gc_cnt = 5, mo->min_gc_score = 1000;
	} else if (strcmp(preset, "se") == 0) {
		io->k = 21, io->w = 10;
		mo->flag |= MG_M_HEAP_SORT;
		mo->mid_occ = 1000;
		mo->max_gap = 100, mo->bw = 100;
		mo->pri_ratio = 0.5f;
		mo->min_lc_cnt = 2, mo->min_lc_score = 25;
		mo->min_gc_cnt = 3, mo->min_gc_score = 40;
		mo->mini_batch_size = 50000000;
	} else return -1;
	return 0;
}

int mg_opt_check(const mg_idxopt_t *io, const mg_mapopt_t *mo, const mg_ggopt_t *go)
{
	return 0;
}

void mg_opt_update(const mg_idx_t *gi, mg_mapopt_t *mo, mg_ggopt_t *go)
{
	int32_t mid_occ;
	mid_occ = mg_idx_cal_max_occ(gi, mo->mid_occ_frac);
	if (mid_occ > mo->mid_occ)
		mo->mid_occ = mid_occ;
	if (mg_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] mid_occ = %d\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), mo->mid_occ);
}
