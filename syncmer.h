/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2022 Chenxi Zhou <chnx.zhou@gmail.com>                          *
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
 * 11/08/22 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#ifndef __SYNCMER_H__
#define __SYNCMER_H__

#include <stdlib.h>
#include <stdint.h>

#include "sstream.h"

extern unsigned char seq_nt4_table[256];
extern char char_nt4_table[4];

#define uint128_t __uint128_t
#define MAX_RD_NUM 0x7FFFFFFFFFULL
#define MAX_RD_LEN 0xFFFFFFULL
#define MAX_RD_SCM 0xFFFFFFULL

typedef struct {
    uint64_t sid; // seq id
    char *sname; // seq name
    uint32_t hoco_l; // hoco seq length
    // homopolymer compressed seq
    // every byte packs 4 bases (00->A, 01->C, 10->G, 11->T)
    // ambiguous bases are converted to 'A'
    uint8_t *hoco_s;
    // homopolymer run length (rl)
    // the real rl of a position with rl=255 is stored in ho_l_rl
    uint8_t *ho_rl;
    // large homopolymer run length (>254)
    // mostly would be NULL
    uint32_t *ho_l_rl;
    // positions of ambiguous bases [hoco space]
    // this first number stores the number of ambiguous bases
    // mostly would be NULL
    uint32_t *n_nucl;
    uint32_t n; // number s/kmer
    uint32_t *m_pos; // s/kmer positions: pos << 1 | rev [hoco space]
    uint64_t *s_mer; // smer << 1 | o/c
    uint128_t *k_mer_h; // kmer hashes
} sr_t;

typedef struct {size_t n, m; sr_t *a;} sr_v;

typedef struct {
    uint128_t h; // kmer hash
    uint64_t s; // smer << 1 | o/c
    // s/kmer positions: sid:39 | pos:24 | rev:1 [hoco space]
    // pos is the syncmer index on reads
    // refer to sr_t m_pos for the real physical position
    uint64_t *m_pos;
    /*** this is expensive and we do not really need this
    char *seq; // consensus sequence
    uint32_t len; // sequence length
    **/
    uint32_t k_cov:31, del:1; // kmer coverage, deleted
} syncmer_t;

typedef struct {
    uint64_t syncmer_n;
    double syncmer_per_read, syncmer_avg_dist, smer_avg_cnt, kmer_avg_cnt;
    int smer_unique, smer_singleton, smer_peak_hom, smer_peak_het;
    int kmer_unique, kmer_singleton, kmer_peak_hom, kmer_peak_het;
} sr_stat_t;

#ifdef __cplusplus
extern "C" {
#endif

void sr_read(sstream_t *s_stream, sr_v *sr, int k, int w, int n_threads);
void print_syncmer_on_seq(sr_t *sr, uint32_t n, int k, int w, FILE *fo);
void print_all_syncmers_on_seq(sr_t *sr, int k, int w, FILE *fo);
void print_seq(sr_t *sr, FILE *fo);
void sr_stat(sr_v *sr, sr_stat_t *stats, int w, FILE *fo, int more);
void sr_destroy(sr_t *sr);
void sr_v_destroy(sr_v *sr_v);
int sr_validate(sr_v *sr);

#ifdef __cplusplus
}
#endif

#endif

