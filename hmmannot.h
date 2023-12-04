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
 * 04/08/22 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#ifndef __HMMANNOT_H__
#define __HMMANNOT_H__

#include <stdint.h>
#include <stdio.h>

static const int athaliana_pltd_g71n = 71;
static const char *const athaliana_pltd_g71[] = {
    "psbA",  "matK",  "rps16", "psbK",  "psbI", "atpA",  "atpF",  "atpH", "atpI",  "rps2",
    "rpoC2", "rpoC1", "rpoB",  "ycf6",  "psbM", "psbD",  "psbC",  "ycf9", "rps14", "psaB",
    "psaA",  "ycf3",  "rps4",  "ndhJ",  "psbG", "ndhC",  "atpE",  "atpB", "rbcL",  "accD",
    "psaI",  "ycf4",  "cemA",  "petA",  "psbJ", "psbL",  "psbF",  "psbE", "ORF31", "petG",
    "psaJ",  "rpl33", "rps18", "rpl20", "clpP", "psbB",  "psbT",  "psbN", "psbH",  "petB",
    "petD",  "rpoA",  "rps11", "rpl36", "rps8", "rpl14", "rpl16", "rps3", "rpl22", "rps19",
    "ndhF",  "rpl32", "ycf5",  "ndhD",  "psaC", "ndhE",  "ndhG",  "ndhI", "ndhA",  "ndhH",
    "rps15"
};

typedef enum organelle_type {
    OG_UNCLASSIFIED = 0,
    OG_MITO = 1,
    OG_PLTD = 2,
    OG_MINI = 3
} OG_TYPE_t;

#define MAX_BED_SCORE 1000

static const char *const OG_TYPES[] = {"unclassified", "mito", "pltd", "mini"};

// flag bit
typedef struct {
    char *gname, *sname;
    uint32_t hmmfrom, hmmto, alifrom, alito, envfrom, envto, modlen;
    double evalue, score, bias;
    uint32_t gid:29, og_type:2, strand:1, sid;
} hmm_annot_t;

typedef enum {
    ORDER_UNKNOWN  =-1,
    ORDER_UNSORTED = 0,
    ORDER_GNAME    = 1, // key order: gname
    ORDER_GID      = 2, // key order: gid
    ORDER_SNAME    = 3, // key order: sname
    ORDER_SID      = 4, // key order: sid
    ORDER_SID_OG   = 5, // key order: sid - og_type - gid - score
    ORDER_SID_CO   = 6  // key order: sid - alifrom - alito
} hmm_annot_so_t;

typedef struct {
    size_t n, m; 
    hmm_annot_t *a;
    uint32_t n_gene, n_seg;
    char **gnames; // gene names
    char **snames; // seg names
    void *h_gnames; // gene name index map
    void *h_snames; // seg name index map
    hmm_annot_so_t so;
    uint32_t n_idx;  // index size: n_gene or n_seg
    uint64_t *index; // index for the first key
} hmm_annot_db_t;

typedef struct {
    char *sname, *gname, strand;
    uint32_t alifrom, alito;
    int score;
} hmm_annot_bed6_t;

typedef struct { 
    size_t n, m; 
    hmm_annot_bed6_t *a; 
    size_t n_seg, m_seg; 
    char **snames; // seg name dictionary
} hmm_annot_bed6_db_t;

#ifdef __cplusplus
extern "C" {
#endif

void hmm_annot_db_destroy(hmm_annot_db_t *annot_db);
hmm_annot_db_t *hmm_annot_read(char *annot_file, hmm_annot_db_t *annot_db, OG_TYPE_t og_type);
void hmm_annot_db_sort(hmm_annot_db_t *annot_db, hmm_annot_so_t so);
int is_trn(hmm_annot_t *annot);
int is_rrn(hmm_annot_t *annot);
hmm_annot_t *hmm_annot_db_index_query(hmm_annot_db_t *annot_db, uint32_t id, uint32_t *n);
hmm_annot_t *hmm_annot_db_sname_query(hmm_annot_db_t *annot_db, char *sname, uint32_t *n);
hmm_annot_t *hmm_annot_db_gname_query(hmm_annot_db_t *annot_db, char *gname, uint32_t *n);
uint32_t hmm_annot_gname2id(hmm_annot_db_t *annot_db, char *gname);
uint32_t hmm_annot_sname2id(hmm_annot_db_t *annot_db, char *sname);
void hmm_annot_print(hmm_annot_t *hmm_annot, size_t n, FILE *fo);
void hmm_annot_db_print(hmm_annot_db_t *annot_db, FILE *fo);
void hmm_annot_formatted_print_index_list(hmm_annot_db_t *annot_db, uint64_t *index_list, size_t n, 
        FILE *fo, OG_TYPE_t og_type, double max_evalue, int header);
void hmm_annot_formatted_print_sname_list(hmm_annot_db_t *annot_db, char **sname_list, size_t n,
        FILE *fo, OG_TYPE_t og_type, double max_evalue, int header);
void hmm_annot_print_bed6(hmm_annot_bed6_db_t *annots, FILE *fo, int header);
void hmm_annot_bed6_sname_add(hmm_annot_bed6_db_t *annots, hmm_annot_db_t *annot_db, char *cname, 
        char *sname, uint32_t len, uint32_t beg, int rev,
        uint32_t offset, OG_TYPE_t og_type, double max_evalue);
hmm_annot_bed6_db_t *hmm_annot_bed6_db_init();
void hmm_annot_bed6_db_destroy(hmm_annot_bed6_db_t *annot_db);

#ifdef __cplusplus
}
#endif

#endif
