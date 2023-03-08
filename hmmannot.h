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

#define OG_UNCLASSIFIED 0
#define OG_MITO 1
#define OG_PLTD 2
#define OG_MINI 3

static const char *const OG_TYPES[] = {"unclassified", "mito", "pltd", "mini"};

// flag bit
typedef struct {
    char *gname, *sname;
    uint32_t hmmfrom, hmmto, alifrom, alito, envfrom, envto, modlen;
    double evalue, score, bias;
    uint32_t gid:27, og_type:2, trn:1, core:1, del:1;
    uint64_t sid:63, strand:1;
} hmm_annot_t;

typedef struct {size_t n, m; hmm_annot_t *a; void *h_names; char **dict; } hmm_annot_v;

#ifdef __cplusplus
extern "C" {
#endif

void hmm_annot_destroy(hmm_annot_t *hmm_annot);
void hmm_annot_v_destroy(hmm_annot_v *annot_v);
hmm_annot_v *hmm_annot_read(char *annot_file, hmm_annot_v *annot_v, uint8_t og_type);
int is_trn(hmm_annot_t *annot);
void hmm_annot_index(hmm_annot_v *annot_v);
uint32_t hmm_annot_name2id(hmm_annot_v *annot_v, char *gname);
void hmm_annot_print(hmm_annot_t *hmm_annot, size_t n, FILE *fo);

#ifdef __cplusplus
}
#endif

#endif
