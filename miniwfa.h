/*
  The MIT License

  Copyright (c) 2022-     Dana-Farber Cancer Institute

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#ifndef MINIWFA_H
#define MINIWFA_H

#include <stdint.h>

#define MWF_F_CIGAR    0x1
#define MWF_F_DEBUG    0x100

typedef struct {
	int32_t flag;
	int32_t x, o1, e1, o2, e2;
	int32_t step;
} mwf_opt_t;

typedef struct {
	int32_t s;       // score
	int32_t n_cigar; // number of CIGAR operators
	uint32_t *cigar; // CIGAR in the htslib packing: len<<4|op
} mwf_rst_t;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Set default parameters
 *
 * @param opt     (out) options
 */
void mwf_opt_init(mwf_opt_t *opt);

/**
 * Align two sequences with WFA
 *
 * @param km      kalloc handler. Set to NULL to use malloc.
 * @param opt     parameters
 * @param tl      target sequence length
 * @param ts      target sequence
 * @param ql      query sequence length
 * @param qs      query sequence
 * @param r       (out) results
 */
void mwf_wfa(void *km, const mwf_opt_t *opt, int32_t tl, const char *ts, int32_t ql, const char *qs, mwf_rst_t *r);

// These functions are in "mwf-dbg.c". For debugging only.
int32_t mwf_cigar2score(const mwf_opt_t *opt, int32_t n_cigar, const uint32_t *cigar, int32_t *tl, int32_t *ql);
void mwf_assert_cigar(const mwf_opt_t *opt, int32_t n_cigar, const uint32_t *cigar, int32_t tl0, int32_t ql0, int32_t s0);

#ifdef __cplusplus
}
#endif

#endif
