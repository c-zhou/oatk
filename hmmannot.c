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
 * 17/10/22 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "kvec.h"
#include "khashl.h"
#include "misc.h"
#include "hmmannot.h"

KHASHL_MAP_INIT(KH_LOCAL, kh_str_t, kh_str, kh_cstr_t, uint32_t, kh_hash_str, kh_eq_str)

void hmm_annot_destroy(hmm_annot_t *hmm_annot)
{
    if (!hmm_annot) return;
    if (hmm_annot->gname) free(hmm_annot->gname);
    if (hmm_annot->sname) free(hmm_annot->sname);
    free(hmm_annot);
}

void hmm_annot_v_destroy(hmm_annot_v *annot_v)
{
    if (!annot_v) return;
    size_t i;
    for (i = 0; i < annot_v->n; ++i) {
        if (annot_v->a[i].gname) free(annot_v->a[i].gname);
        if (annot_v->a[i].sname) free(annot_v->a[i].sname);
    }
    free(annot_v->a);
    if (annot_v->h_names) {
        kh_str_t *h_names = (kh_str_t *) annot_v->h_names;
        khint_t n = kh_size(h_names);
        for (i = 0; i < n; ++i)
            free(annot_v->dict[i]);
        free(annot_v->dict);
        kh_str_destroy(h_names);
    }
    free(annot_v);
}

static int is_empty_line(char *line)
{
    char *c;
    c = line;
    while (isspace((unsigned char) *c))
        c++;
    if(*c == 0)
        return 1;
    return 0;
}

hmm_annot_v *hmm_annot_read(char *annot_file, hmm_annot_v *annot_v, uint8_t og_type) // og_type: OG_MITO or OG_PLTD or OG_MINI
{
    if (!annot_file) return 0;

    FILE *fp;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char gname[512], sname[512], strand[4];
    hmm_annot_t *annot;

    if (!annot_v) {
        MYMALLOC(annot_v, 1);
        annot_v->n = annot_v->m = 0;
        annot_v->a = 0;
        annot_v->h_names = 0;
        annot_v->dict = 0;
    }

    fp = fopen(annot_file, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, annot_file);
        exit(EXIT_FAILURE);
    }

    while ((read = getline(&line, &ln, fp)) != -1) {
        if (is_empty_line(line) || !strncmp(line, "#", 1))
            continue;
        
        kv_pushp(hmm_annot_t, *annot_v, &annot);
        
        sscanf(line, "%s %*s %s %*s %u %u %u %u %u %u %u %s %lf %lf %lf %*s",
                gname, sname, &annot->hmmfrom, &annot->hmmto, &annot->alifrom, 
                &annot->alito, &annot->envfrom, &annot->envto, &annot->modlen, 
                strand, &annot->evalue, &annot->score, &annot->bias);

        annot->gname = strdup(gname);
        annot->sname = strdup(sname);
        annot->strand = strand[0]=='+'? 0 : 1;
        annot->og_type = og_type;
    }
    
    if (line) free(line);
    fclose(fp);

    return annot_v;
}

void hmm_annot_index(hmm_annot_v *annot_v)
{
    // make gene dictionary and update annotation gid
    int absent;
    size_t i;
    uint32_t gid;
    kh_str_t *h_names;
    khint32_t k, n;
    char **dict;

    if (annot_v->h_names)
        kh_str_destroy((kh_str_t *) annot_v->h_names);

    h_names = kh_str_init();
    gid = 0;
    for (i = 0; i < annot_v->n; ++i) {    
        k = kh_str_put(h_names, annot_v->a[i].gname, &absent);
        if (absent) kh_val(h_names, k) = gid++;
        annot_v->a[i].gid = kh_val(h_names, k);
    }
    n = kh_size(h_names);
    MYMALLOC(dict, n);
    for (k = (khint32_t) 0; k < kh_end(h_names); ++k) {
        if (kh_exist(h_names, k))
            dict[kh_val(h_names, k)] = strdup(kh_key(h_names, k));
    }

    annot_v->h_names = h_names;
    annot_v->dict = dict;
}

uint32_t hmm_annot_name2id(hmm_annot_v *annot_v, char *gname)
{
    kh_str_t *h_names;
    khint32_t k;
    h_names = (kh_str_t *) annot_v->h_names;
    k = kh_str_get(h_names, gname);
    return k < kh_end(h_names)? kh_val(h_names, k) : UINT32_MAX;
}

void hmm_annot_print(hmm_annot_t *hmm_annot, size_t n, FILE *fo)
{
    size_t i;
    hmm_annot_t *annot;
    for (i = 0; i < n; ++i) {
        annot = &hmm_annot[i];
        fprintf(fo, "[M::%s] %s %s %u %u %u %u %u %u %u %e %.1f %.1f\n", __func__, 
                annot->gname, annot->sname, annot->hmmfrom, annot->hmmto, annot->alifrom, annot->alito, 
                annot->envfrom, annot->envto, annot->modlen, annot->evalue, annot->score, annot->bias);
    }
}

int is_trn(hmm_annot_t *annot)
{
    return !strncmp(annot->gname, "trn", 3);
}

