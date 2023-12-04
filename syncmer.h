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

#include "kstring.h"
#include "sstream.h"

extern const unsigned char seq_nt4_table[256];
extern const char seq_nt4_comp_table[128];
extern const char char_nt4_table[4];

#define uint128_t __uint128_t
#define MAX_RD_NUM 0xFFFFFFFFULL
#define MAX_RD_LEN 0x7FFFFFFFULL
#define MAX_RD_SCM 0x7FFFFFFFULL

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
    uint64_t *k_mer; // kmer id << 1 | ec error corrected
} sr_t;

typedef struct {
    uint64_t syncmer_n;
    double syncmer_per_read, syncmer_avg_dist, smer_avg_cnt, kmer_avg_cnt;
    int smer_unique, smer_singleton, smer_peak_hom, smer_peak_het;
    int kmer_unique, kmer_singleton, kmer_peak_hom, kmer_peak_het;
} sr_stat_t;

typedef struct {
    size_t n, m;
    sr_t *a;
    int k, s; // kmer and smer size
    sr_stat_t *stats;
} sr_db_t;

typedef struct {
    uint64_t h, s; // kmer hash and smer
    // this is expensive and we do not really need this
    // char *seq; // consensus sequence
    // uint32_t len; // sequence length
    uint32_t cov:31, del:1; // kmer coverage, deleted
    // s/kmer positions: sid:32 | pos:31 | rev:1 [hoco space]
    // pos is the syncmer index on reads
    // refer to sr_t m_pos for the real physical position
    uint64_t *m_pos;
} syncmer_t;

// syncmer database
typedef struct {
    size_t n, m;
    syncmer_t *a; // sorted syncmer
    // different copies of the same syncmer
    // first and last 8 bits indicate number of copies before [excluded]
    // and after [included] the current position
    // any number larger than 254 [included] will be 255
    // the total number of copies is (uint32_t) (c>>8) + (uint8_t) c
    // the start and end position is [p - (c>>8), p + (uint8_t) c - 1]
    uint16_t *c;
    // haplotype/reference infomration
    // e.g., the haplotype assignment
    // e.g., the reference position the syncmer mapped to
    // 31 bits chr + 32 bits pos + 1 bit strand
    uint64_t *h;
} syncmer_db_t;

#ifdef __cplusplus
extern "C" {
#endif

void sr_read(sstream_t *s_stream, sr_db_t *sr_db, size_t m_data, int n_threads);
syncmer_db_t *collect_syncmer_from_reads(sr_db_t *sr_db);
int syncmer_link_coverage_analysis(sr_db_t *sr_db, syncmer_db_t *scm_db, uint32_t min_k_cov, 
        uint32_t min_n_seq, uint32_t min_pt, double min_f, double **_beta, 
        double **_bse, double **_r2, int verbose);
void print_syncmer_on_seq(sr_t *sr, uint32_t n, int k, int w, FILE *fo);
void print_all_syncmers_on_seq(sr_t *sr, int k, int w, FILE *fo);
void print_aligned_syncmers_on_seq(sr_t *sr, int w, uint32_t beg, uint32_t end, FILE *fo);
void print_hoco_seq(sr_t *sr, FILE *fo);
void get_hoco_seq(sr_t *sr, kstring_t *s);
void get_kmer_seq(uint8_t *hoco_s, uint32_t pos, int l, uint32_t rev, uint8_t *kmer_s);
void get_kmer_dna_seq(uint8_t *hoco_s, uint32_t pos, int l, uint32_t rev, char *dna_seq);
void sr_destroy(sr_t *sr);
void sr_db_init(sr_db_t *sr_db, int k, int s);
void sr_db_clean(sr_db_t *sr_db);
void sr_db_destroy(sr_db_t *sr_db);
int sr_db_validate(sr_db_t *sr_db);
void sr_db_stat(sr_db_t *sr_db, FILE *fo, int more);
void syncmer_db_init(syncmer_db_t *scm_db);
void syncmer_db_clean(syncmer_db_t *scm_db);
void syncmer_db_destroy(syncmer_db_t *scm_db);

#ifdef __cplusplus
}
#endif

#endif

