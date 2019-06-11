#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bseq.h"
#include "mgpriv.h"
#include "ketopt.h"

#ifdef __linux__
#include <sys/resource.h>
#include <sys/time.h>
void liftrlimit()
{
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
}
#else
void liftrlimit() {}
#endif

static ko_longopt_t long_options[] = {
	{ "print-gfa",    ko_no_argument,       301 },
	{ "no-kalloc",    ko_no_argument,       302 },
	{ "vcoor",        ko_no_argument,       303 },
	{ 0, 0, 0 }
};

static inline int64_t mg_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm') x *= 1e6;
	else if (*p == 'K' || *p == 'k') x *= 1e3;
	return (int64_t)(x + .499);
}

int main(int argc, char *argv[])
{
	const char *opt_str = "x:k:w:t:r:m:n:g:K:o:p:N:Pq:d:l:";
	ketopt_t o = KETOPT_INIT;
	mg_mapopt_t opt;
	mg_idxopt_t ipt;
	mg_ggopt_t gpt;
	int i, c, n_threads = 4, print_gfa = 0;
//	char *rg = 0;
	char *s;
	FILE *fp_help = stderr;
	mg_idx_t *gi;

	mg_verbose = 3;
	liftrlimit();
	mg_realtime0 = realtime();
	mg_opt_set(0, &ipt, &opt, &gpt);

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) { // test command line options and apply option -x/preset first
		if (c == 'x') {
			if (mg_opt_set(o.arg, &ipt, &opt, &gpt) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", o.arg);
				return 1;
			}
		} else if (c == ':') {
			fprintf(stderr, "[ERROR] missing option argument\n");
			return 1;
		} else if (c == '?') {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[o.i - 1]);
			return 1;
		}
	}
	o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'w') ipt.w = atoi(o.arg);
		else if (c == 'k') ipt.k = atoi(o.arg);
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'r') opt.bw = mg_parse_num(o.arg);
		else if (c == 'g') opt.max_gap = mg_parse_num(o.arg);
		else if (c == 'K') opt.mini_batch_size = mg_parse_num(o.arg);
		else if (c == 'p') opt.pri_ratio = atof(o.arg);
		else if (c == 'N') opt.best_n = mg_parse_num(o.arg);
		else if (c == 'P') opt.flag |= MG_F_ALL_CHAINS;
		else if (c == 'l') gpt.min_map_len = mg_parse_num(o.arg);
		else if (c == 'd') gpt.min_depth_len = mg_parse_num(o.arg);
		else if (c == 'q') gpt.min_mapq = atoi(o.arg);
		else if (c == 301) print_gfa = 1;
		else if (c == 302) mg_dbg_flag |= MG_DBG_NO_KALLOC;
		else if (c == 303) opt.flag |= MG_M_VERTEX_COOR;
		else if (c == 'n') {
			opt.min_gc_cnt = strtol(o.arg, &s, 10);
			if (*s == ',') opt.min_lc_cnt = strtol(s + 1, &s, 10);
		} else if (c == 'm') {
			opt.min_gc_score = strtol(o.arg, &s, 10);
			if (*s == ',') opt.min_lc_score = strtol(s + 1, &s, 10);
		} else if (c == 'o') {
			if (strcmp(o.arg, "-") != 0) {
				if (freopen(o.arg, "wb", stdout) == NULL) {
					fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m\n", o.arg);
					exit(1);
				}
			}
		}
	}
	if (mg_opt_check(&ipt, &opt, &gpt) < 0)
		return 1;

	if (argc == o.ind || fp_help == stdout) {
		fprintf(fp_help, "Usage: minigraph [options] <target.gfa> <query.fa> [...]\n");
		fprintf(fp_help, "Options:\n");
		fprintf(fp_help, "  Indexing:\n");
		fprintf(fp_help, "    -k INT       k-mer size (no larger than 28) [%d]\n", ipt.k);
		fprintf(fp_help, "    -w INT       minizer window size [%d]\n", ipt.w);
		fprintf(fp_help, "  Mapping:\n");
		fprintf(fp_help, "    -g NUM       stop chain enlongation if there are no minimizers in INT-bp [%d]\n", opt.max_gap);
		fprintf(fp_help, "    -r NUM       bandwidth used in chaining and DP-based alignment [%d]\n", opt.bw);
		fprintf(fp_help, "    -n INT[,INT] minimal number of minimizers on a graph/linear chain [%d,%d]\n", opt.min_gc_cnt, opt.min_lc_cnt);
		fprintf(fp_help, "    -m INT[,INT] minimal graph/linear chaining score [%d,%d]\n", opt.min_gc_score, opt.min_lc_score);
		fprintf(fp_help, "    -p FLOAT     min secondary-to-primary score ratio [%g]\n", opt.pri_ratio);
		fprintf(fp_help, "    -N INT       retain at most INT secondary alignments [%d]\n", opt.best_n);
		fprintf(fp_help, "  Graph generation:\n");
		fprintf(fp_help, "    -q INT       min mapping quality [%d]\n", gpt.min_mapq);
		fprintf(fp_help, "    -l NUM       min alignment length [%d]\n", gpt.min_map_len);
		fprintf(fp_help, "    -d NUM       min alignment length for depth calculation [%d]\n", gpt.min_depth_len);
		fprintf(fp_help, "  Input/output:\n");
		fprintf(fp_help, "    -t INT       number of threads [%d]\n", n_threads);
		fprintf(fp_help, "    -o FILE      output alignments to FILE [stdout]\n");
		fprintf(fp_help, "    -K NUM       minibatch size for mapping [500M]\n");
		fprintf(fp_help, "  Preset:\n");
		fprintf(fp_help, "    -x STR       preset []\n");
		fprintf(fp_help, "                 - ggsimple: simple algorithm for graph generation\n");
		return fp_help == stdout? 0 : 1;
	}

	gi = mg_index_file(argv[o.ind], ipt.k, ipt.w, ipt.bucket_bits, n_threads);
	if (gi == 0) {
		fprintf(stderr, "[ERROR] failed to load the graph from file '%s'\n", argv[o.ind]);
		return 1;
	}
	if (print_gfa) {
		gfa_print(gi->g, stdout, 1);
		goto free_gfa;
	}

#if 0
	int sid1 = gfa_name2id(gi->g, "MTh0");
	int sid2 = gfa_name2id(gi->g, "MTh13516");
	int sid3 = gfa_name2id(gi->g, "MTo8961");
	int32_t n_pathv;
	gfa_path_dst_t dst[3];
	gfa_pathv_t *path;
	if (sid1 < 0 || sid2 < 0) abort();
	dst[0].v = sid2<<1|0, dst[0].target_dist = 13516;
	dst[1].v = sid3<<1|0, dst[1].target_dist = 10000;
	path = gfa_shortest_k(0, gi->g, sid1<<1|0, 2, dst, 20000, 7, &n_pathv);
	gfa_sub_print_path(stderr, gi->g, n_pathv, path);
	free(path);
#endif

	if (gpt.algo == MG_G_NONE) {
		for (i = o.ind + 1; i < argc; ++i)
			mg_map_file(gi, argv[i], &opt, n_threads);
	} else {
		for (i = o.ind + 1; i < argc; ++i)
			mg_ggen(gi, argv[i], &opt, &gpt, n_threads);
	}

free_gfa:
	mg_idx_destroy(gi);

	if (fflush(stdout) == EOF) {
		fprintf(stderr, "[ERROR] failed to write the results\n");
		exit(EXIT_FAILURE);
	}

	if (mg_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MG_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - mg_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
