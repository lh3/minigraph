#include <string.h>
#include "minigraph.h"

void mg_idxopt_init(mg_idxopt_t *io)
{
	memset(io, 0, sizeof(mg_idxopt_t));
	io->k = 15;
	io->w = 10;
	io->bucket_bits = 14;
}

void mg_mapopt_init(mg_mapopt_t *mo)
{
	memset(mo, 0, sizeof(mg_idxopt_t));
	mo->seed = 11;
	mo->mid_occ = 200;
	mo->max_gap = 5000;
	mo->max_gap_ref = -1;
	mo->max_chain_skip = 25;
	mo->bw = 500;
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
	go->min_map_len = 50000;
}

int mg_opt_set(const char *preset, mg_idxopt_t *io, mg_mapopt_t *mo, mg_ggopt_t *go)
{
	if (preset == 0) {
		mg_idxopt_init(io);
		mg_mapopt_init(mo);
		mg_ggopt_init(go);
	} else if (strcmp(preset, "ggsimple") == 0) {
		go->algo = MG_G_GGSIMPLE;
		mo->best_n = 0;
	} else return -1;
	return 0;
}

int mg_opt_check(const mg_idxopt_t *io, const mg_mapopt_t *mo, const mg_ggopt_t *go)
{
	return 0;
}
