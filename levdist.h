/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2023 Chenxi Zhou <chnx.zhou@gmail.com>                          *
 *                                                                               *
 * Permission is hereby granted, free of charge, to any person obtaining a copy  *
 * of this software and associated documentation files (the "Software"), to deal *
 * in the Software without restriction, including without limitation the rights  *
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     *
 * copies of the Software, and to permit persons to whom the Software is         *
 * furnished to do so, subject to the following conditions:                      *
 *                                                                               *
 * The above copyright notice and this permission notice shall be included in    *
 * all copies or substantial portions of the Software.                           *
 *                                                                               *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, *
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE *
 * SOFTWARE.                                                                     *
 *********************************************************************************/

/********************************** Revision History *****************************
 *                                                                               *
 * 24/07/23 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#ifndef LEV_DIST_H
#define LEV_DIST_H

#include <stdlib.h>
#include <stdint.h>

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/*
 * Diagonal and a couple of common routines
 */
typedef struct {
    int32_t d, k, p; // diagonal (query pos minus target pos), target pos, traceback variable
} wf_diag1_t;

typedef struct {
    size_t n, m;
    wf_diag1_t *a;
} wf_diag_t;

/*
 * Traceback arrays
 */
typedef struct {
    int32_t n, d0;
    uint64_t *a;
} wf_tb1_t;

typedef struct {
    size_t m, n;
    wf_tb1_t *a;
} wf_tb_t;

/*
 * CIGAR operations
 */
typedef struct {
    int32_t m, n;
    uint32_t *cigar;
} wf_cigar_t;

typedef struct {
    char *ts, *qs;
    int32_t tl, ql;
    wf_diag_t *wf_diag;
    wf_tb_t *wf_tb;
    int32_t is_ext, bw;
    int32_t score, t_end, q_end, n_cigar;
    uint32_t *cigar;
} wf_config_t;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Find the edit distance between two sequences
 *
 * @param tl         target sequence length
 * @param ts         target sequence
 * @param ql         query sequence length
 * @param qs         query sequence
 * @param is_ext     extension alignment if true (stop when reaching the end of either query or target)
 * @param bw         bandwidth
 * @param step       step size; only when n_cigar is not NULL
 * @param score      (out) edit distance
 * @param t_endl     (out) length of target in the alignment
 * @param q_endl     (out) length of query in the alignment
 * @param n_cigar    (in/out) number of cigar operations; NULL if don't need CIGAR
 *
 * @return CIGAR in the htslib packing
 */
uint32_t *wf_ed(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t is_ext, int32_t bw, int32_t *score, 
        int32_t *t_endl, int32_t *q_endl, int32_t *n_cigar);
void wf_ed_core(wf_config_t *conf);
void wf_print_cigar(uint32_t *cigar, int32_t n_cigar, FILE *fo);
void wf_print_alignment(char *ts, int tl, char *qs, int ql, uint32_t *cigar, int32_t n_cigar, int lwd, FILE *fo);
void wf_config_destroy(wf_config_t *conf, int clean_seq);

#ifdef __cplusplus
}
#endif

#endif

