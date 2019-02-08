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
	mo->bw = 500;
}

int mg_opt_set(const char *preset, mg_idxopt_t *io, mg_mapopt_t *mo)
{
	if (preset == 0) {
		mg_idxopt_init(io);
		mg_mapopt_init(mo);
	} else return -1;
	return 0;
}

int mg_opt_check(const mg_idxopt_t *io, const mg_mapopt_t *mo)
{
	return 0;
}
