#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mgpriv.h"
#include "gfa-priv.h"
#include "sys.h"
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
	{ "version",      ko_no_argument,       300 },
	{ "vc",           ko_no_argument,       301 },
	{ "secondary",    ko_required_argument, 302 },
	{ "ins-qovlp",    ko_required_argument, 303 },
	{ "heap-sort",    ko_required_argument, 304 },
	{ "show-unmap",   ko_required_argument, 305 },
	{ "ggen",         ko_optional_argument, 306 },
	{ "gg-vkmer",     ko_required_argument, 307 },
	{ "gg-vkiden",    ko_required_argument, 308 },
	{ "gg-min-end-cnt",  ko_required_argument, 309 },
	{ "gg-min-end-frac", ko_required_argument, 310 },
	{ "no-comp-path", ko_no_argument,       312 },
	{ "gg-match-pen", ko_required_argument, 313 },
	{ "frag",         ko_no_argument,       314 },
	{ "cov",          ko_no_argument,       315 },
	{ "min-cov-blen", ko_required_argument, 316 },
	{ "min-cov-mapq", ko_required_argument, 317 },
	{ "gap-pen",      ko_required_argument, 318 },
	{ "no-kalloc",    ko_no_argument,       401 },
	{ "dbg-qname",    ko_no_argument,       402 },
	{ "dbg-lchain",   ko_no_argument,       403 },
	{ "dbg-insert",   ko_no_argument,       404 },
	{ "dbg-shortk",   ko_no_argument,       405 },
	{ "dbg-gc1",      ko_no_argument,       406 },
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

static inline void yes_or_no(uint64_t *flag_, int f, int long_idx, const char *arg, int yes_to_set)
{
	uint64_t flag = *flag_;
	if (yes_to_set) {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) flag |= f;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) flag &= ~f;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	} else {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) flag &= ~f;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) flag |= f;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	}
	*flag_ = flag;
}

int main(int argc, char *argv[])
{
	const char *opt_str = "x:k:w:t:r:m:n:g:K:o:p:N:Pq:d:l:f:U:M:F:j:L:";
	ketopt_t o = KETOPT_INIT;
	mg_mapopt_t opt;
	mg_idxopt_t ipt;
	mg_ggopt_t gpt;
	int i, c, ret, n_threads = 4;
	char *s;
	FILE *fp_help = stderr;
	gfa_t *g;

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
		else if (c == 'f') opt.mid_occ_frac = atof(o.arg);
		else if (c == 'U') opt.mid_occ = atoi(o.arg);
		else if (c == 'r') opt.bw = mg_parse_num(o.arg);
		else if (c == 'g') opt.max_gap = mg_parse_num(o.arg);
		else if (c == 'F') opt.max_frag_len = mg_parse_num(o.arg);
		else if (c == 'K') opt.mini_batch_size = mg_parse_num(o.arg);
		else if (c == 'p') opt.pri_ratio = atof(o.arg);
		else if (c == 'N') opt.best_n = mg_parse_num(o.arg);
		else if (c == 'P') opt.flag |= MG_M_ALL_CHAINS;
		else if (c == 'M') opt.mask_level = atof(o.arg);
		else if (c == 'j') opt.div = atof(o.arg);
		else if (c == 'l') gpt.min_map_len = mg_parse_num(o.arg);
		else if (c == 'd') gpt.min_depth_len = mg_parse_num(o.arg);
		else if (c == 'q') gpt.min_mapq = atoi(o.arg);
		else if (c == 'L') gpt.min_var_len = atoi(o.arg);
		else if (c == 301) opt.flag |= MG_M_VERTEX_COOR;      // --vc
		else if (c == 307) gpt.ggs_fc_kmer = atoi(o.arg);     // --gg-vkmer
		else if (c == 308) gpt.ggs_max_kiden = atof(o.arg);   // --gg-vkiden
		else if (c == 309) gpt.ggs_min_end_cnt = atoi(o.arg);  // --gg-min-end-cnt
		else if (c == 310) gpt.ggs_min_end_frac = atof(o.arg); // --gg-min-end-frac
		else if (c == 312) opt.flag |= MG_M_NO_COMP_PATH;     // --no-comp-path
		else if (c == 313) gpt.match_pen = atoi(o.arg);       // --gg-match-pen
		else if (c == 314) opt.flag |= MG_M_FRAG_MODE | MG_M_FRAG_MERGE;       // --frag
		else if (c == 315) opt.flag |= MG_M_CAL_COV, gpt.flag |= MG_G_CAL_COV; // --cov
		else if (c == 316) opt.min_cov_blen = mg_parse_num(o.arg);             // --min-cov-blen
		else if (c == 317) opt.min_cov_mapq = atoi(o.arg);                     // --min-cov-mapq
		else if (c == 318) opt.chn_pen_gap = atof(o.arg);     // --gap-pen
		else if (c == 401) mg_dbg_flag |= MG_DBG_NO_KALLOC;   // --no-kalloc
		else if (c == 402) mg_dbg_flag |= MG_DBG_QNAME;       // --dbg-qname
		else if (c == 403) mg_dbg_flag |= MG_DBG_LCHAIN;      // --dbg-lchain
		else if (c == 404) mg_dbg_flag |= MG_DBG_INSERT;      // --dbg-insert
		else if (c == 405) mg_dbg_flag |= MG_DBG_SHORTK;      // --dbg-shortk
		else if (c == 406) mg_dbg_flag |= MG_DBG_GC1;         // --dbg-gc1
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
		} else if (c == 306) { // --ggen
			if (o.arg) {
				if (strcmp(o.arg, "none") == 0) gpt.algo = MG_G_NONE;
				else if (strcmp(o.arg, "simple") == 0) gpt.algo = MG_G_GGSIMPLE;
				else {
					fprintf(stderr, "ERROR: unknown graph generation algorithm \"%s\"\n", o.arg);
					return 1;
				}
			} else gpt.algo = MG_G_GGSIMPLE;
		} else if (c == 302) { // --secondary
			yes_or_no(&opt.flag, MG_M_PRINT_2ND, o.longidx, o.arg, 1);
		} else if (c == 303) { // --ins-qovlp
			yes_or_no(&gpt.flag, MG_G_NO_QOVLP, o.longidx, o.arg, 1);
		} else if (c == 304) { // --heap-sort
			yes_or_no(&opt.flag, MG_M_HEAP_SORT, o.longidx, o.arg, 1);
		} else if (c == 305) { // --show-unmap
			yes_or_no(&opt.flag, MG_M_SHOW_UNMAP, o.longidx, o.arg, 1);
		} else if (c == 300) { // --version
			puts(MG_VERSION);
			return 0;
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
		fprintf(fp_help, "    -f FLOAT     ignore top FLOAT fraction of repetitive minimizers [%g]\n", opt.mid_occ_frac);
		fprintf(fp_help, "    -U INT       ignore minimizers with occurrences above INT [%d]\n", opt.mid_occ);
		fprintf(fp_help, "    -j FLOAT     expected sequence divergence [%g]\n", opt.div);
		fprintf(fp_help, "    -g NUM       stop chain enlongation if there are no minimizers in INT-bp [%d]\n", opt.max_gap);
		fprintf(fp_help, "    -F NUM       max fragment length (effective with -xsr or in the fragment mode) [%d]\n", opt.max_frag_len);
		fprintf(fp_help, "    -r NUM       bandwidth used in chaining [%d]\n", opt.bw);
		fprintf(fp_help, "    -n INT[,INT] minimal number of minimizers on a graph/linear chain [%d,%d]\n", opt.min_gc_cnt, opt.min_lc_cnt);
		fprintf(fp_help, "    -m INT[,INT] minimal graph/linear chaining score [%d,%d]\n", opt.min_gc_score, opt.min_lc_score);
		fprintf(fp_help, "    -p FLOAT     min secondary-to-primary score ratio [%g]\n", opt.pri_ratio);
		fprintf(fp_help, "    -N INT       retain at most INT secondary mappings [%d]\n", opt.best_n);
		fprintf(fp_help, "    --cov        remap and generate segment and link use frequencies; output GFA\n");
		fprintf(fp_help, "  Graph generation:\n");
		fprintf(fp_help, "    --ggen       perform incremental graph generation\n");
		fprintf(fp_help, "    -q INT       min mapping quality [%d]\n", gpt.min_mapq);
		fprintf(fp_help, "    -l NUM       min alignment length [%d]\n", gpt.min_map_len);
		fprintf(fp_help, "    -d NUM       min alignment length for depth calculation [%d]\n", gpt.min_depth_len);
		fprintf(fp_help, "    -L INT       min variant length [%d]\n", gpt.min_var_len);
		fprintf(fp_help, "  Input/output:\n");
		fprintf(fp_help, "    -t INT       number of threads [%d]\n", n_threads);
		fprintf(fp_help, "    -o FILE      output mappings to FILE [stdout]\n");
		fprintf(fp_help, "    -K NUM       minibatch size for mapping [500M]\n");
		fprintf(fp_help, "  Preset:\n");
		fprintf(fp_help, "    -x STR       preset []\n");
		fprintf(fp_help, "                 - lr: noisy long read mapping (the default)\n");
		fprintf(fp_help, "                 - asm20: asm-to-ref mapping, for ~5%% sequence divergence\n");
		fprintf(fp_help, "                 - sr: short reads\n");
		fprintf(fp_help, "                 - ggs: simple algorithm for graph generation\n");
		return fp_help == stdout? 0 : 1;
	}

	g = gfa_read(argv[o.ind]);
	if (g == 0) {
		fprintf(stderr, "[ERROR] failed to load the graph from file '%s'\n", argv[o.ind]);
		return 1;
	} else if (mg_verbose >= 3) {
		fprintf(stderr, "[M::%s::%.3f*%.2f] loaded the graph from \"%s\"\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), argv[o.ind]);
	}

	if (gpt.algo == MG_G_NONE) {
		ret = mg_map_files(g, argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &ipt, &opt, n_threads);
	} else {
		ret = mg_ggen(g, argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &ipt, &opt, &gpt, n_threads);
	}

	if (gpt.algo != MG_G_NONE || (opt.flag & MG_M_CAL_COV))
		gfa_print(g, stdout, 0);
	gfa_destroy(g);

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
	return !!ret;
}
