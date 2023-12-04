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
#include <math.h>

#include "kvec.h"
#include "khashl.h"
#include "misc.h"
#include "hmmannot.h"

KHASHL_MAP_INIT(KH_LOCAL, kh_str_t, kh_str, kh_cstr_t, uint32_t, kh_hash_str, kh_eq_str)

void hmm_annot_db_destroy(hmm_annot_db_t *annot_db)
{
    if (!annot_db) return;
    size_t i;
    kh_str_t *h_gnames = (kh_str_t *) annot_db->h_gnames;
    kh_str_destroy(h_gnames);
    for (i = 0; i < annot_db->n_gene; ++i)
        free(annot_db->gnames[i]);
    free(annot_db->gnames);
    kh_str_t *h_snames = (kh_str_t *) annot_db->h_snames;
    kh_str_destroy(h_snames);
    for (i = 0; i < annot_db->n_seg; ++i)
        free(annot_db->snames[i]);
    free(annot_db->snames);
    free(annot_db->a);
    free(annot_db->index);
    free(annot_db);
}

hmm_annot_bed6_db_t *hmm_annot_bed6_db_init()
{
    hmm_annot_bed6_db_t *annot_db;
    MYCALLOC(annot_db, 1);
    return annot_db;
}

void hmm_annot_bed6_db_destroy(hmm_annot_bed6_db_t *annot_db)
{
    size_t i;
    for (i = 0; i < annot_db->n_seg; ++i)
        free(annot_db->snames[i]);
    free(annot_db->snames);
    free(annot_db->a);
    free(annot_db);
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

static inline char **update_name_list(char *name, char **dict, size_t *_n, size_t *_m, kh_str_t *h_dict, uint32_t *index)
{
    khint32_t k;

    k = kh_str_get(h_dict, name);
    if (k < kh_end(h_dict)) {
        // name already exists
        *index = kh_val(h_dict, k);
        return dict;
    }

    // new name
    int absent;
    size_t n, m;
    n = *_n;
    m = *_m;
    if (n == m) {
        ++m;
        kroundup32(m);
        MYREALLOC(dict, m);
    }
    dict[n] = strdup(name);
    k = kh_str_put(h_dict, dict[n], &absent);
    kh_val(h_dict, k) = n;

    *_n = n + 1;
    *_m = m;
    *index = n;

    return dict;
}

hmm_annot_db_t *hmm_annot_read(char *annot_file, hmm_annot_db_t *annot_db, OG_TYPE_t og_type) // og_type: OG_MITO or OG_PLTD or OG_MINI
{
    if (!annot_file) return 0;

    FILE *fp;
    char *line = NULL;
    size_t ln = 0;
    ssize_t read;
    char gname[512], sname[512], strand[4];
    hmm_annot_t *annot;
    kh_str_t *h_gnames, *h_snames;
    char **gnames, **snames;
    size_t n_gene, m_gene, n_seg, m_seg;
    uint32_t index;

    if (annot_db) {
        n_gene = m_gene = annot_db->n_gene;
        n_seg = m_seg = annot_db->n_seg;
        gnames = annot_db->gnames;
        snames = annot_db->snames;
        h_gnames = annot_db->h_gnames;
        h_snames = annot_db->h_snames;
    } else {
        MYCALLOC(annot_db, 1);
        n_gene = n_seg = 0;
        m_gene = m_seg = 256;
        MYMALLOC(gnames, m_gene);
        MYMALLOC(snames, m_seg);
        h_gnames = kh_str_init();
        h_snames = kh_str_init();
    }

    fp = fopen(annot_file, "r");
    if (fp == NULL) {
        fprintf(stderr, "[E::%s] cannot open file %s for reading\n", __func__, annot_file);
        exit(EXIT_FAILURE);
    }

    while ((read = getline(&line, &ln, fp)) != -1) {
        if (is_empty_line(line) || !strncmp(line, "#", 1))
            continue;
        
        kv_pushp(hmm_annot_t, *annot_db, &annot);
        
        sscanf(line, "%s %*s %s %*s %u %u %u %u %u %u %u %s %lf %lf %lf %*s",
                gname, sname, &annot->hmmfrom, &annot->hmmto, &annot->alifrom, 
                &annot->alito, &annot->envfrom, &annot->envto, &annot->modlen, 
                strand, &annot->evalue, &annot->score, &annot->bias);

        annot->strand = strand[0]=='+'? 0 : 1;
        annot->og_type = (uint32_t) og_type;

        if (annot->strand) {
            // reverse strand
            SWAP(annot->alifrom, annot->alito);
            SWAP(annot->envfrom, annot->envto);
        }

        snames = update_name_list(sname, snames, &n_seg, &m_seg, h_snames, &index);
        annot->sid = index;
        annot->sname = snames[index];
        gnames = update_name_list(gname, gnames, &n_gene, &m_gene, h_gnames, &index);
        annot->gid = index;
        annot->gname = gnames[index];
    }
    
    if (line) free(line);
    fclose(fp);
    MYREALLOC(gnames, n_gene);
    MYREALLOC(snames, n_seg);
    annot_db->n_gene = n_gene;
    annot_db->n_seg = n_seg;
    annot_db->gnames = gnames;
    annot_db->snames = snames;
    annot_db->h_gnames = h_gnames;
    annot_db->h_snames = h_snames;

    return annot_db;
}

uint32_t hmm_annot_gname2id(hmm_annot_db_t *annot_db, char *gname)
{
    kh_str_t *h_gnames;
    khint32_t k;
    h_gnames = (kh_str_t *) annot_db->h_gnames;
    k = kh_str_get(h_gnames, gname);
    return k < kh_end(h_gnames)? kh_val(h_gnames, k) : UINT32_MAX;
}

uint32_t hmm_annot_sname2id(hmm_annot_db_t *annot_db, char *sname)
{
    kh_str_t *h_snames;
    khint32_t k;
    h_snames = (kh_str_t *) annot_db->h_snames;
    k = kh_str_get(h_snames, sname);
    return k < kh_end(h_snames)? kh_val(h_snames, k) : UINT32_MAX;
}

static inline void hmm_annot_print_all_fields(hmm_annot_t *annot, FILE *fo)
{
    fprintf(fo, "%s [%u] %s [%u] %u %u %u %u %u %u %u %c %e %.1f %.1f\n",
            annot->gname, annot->gid, annot->sname, (uint32_t) annot->sid, 
            annot->hmmfrom, annot->hmmto, annot->alifrom, annot->alito,
            annot->envfrom, annot->envto, annot->modlen, annot->strand? '-' : '+', 
            annot->evalue, annot->score, annot->bias);
}

void hmm_annot_print(hmm_annot_t *hmm_annot, size_t n, FILE *fo)
{
    size_t i;
    for (i = 0; i < n; ++i) hmm_annot_print_all_fields(&hmm_annot[i], fo);
}

void hmm_annot_db_print(hmm_annot_db_t *annot_db, FILE *fo)
{
    hmm_annot_print(annot_db->a, annot_db->n, fo);
}

static int hmm_annot_cmpfunc1(const void *a, const void *b)
{
    return strcmp(((hmm_annot_t *) a)->gname, ((hmm_annot_t *) b)->gname);
}

static int hmm_annot_cmpfunc2(const void *a, const void *b)
{
    uint32_t x, y;
    x = ((hmm_annot_t *) a)->gid;
    y = ((hmm_annot_t *) b)->gid;
    return (x > y) - (x < y);
}

static int hmm_annot_cmpfunc3(const void *a, const void *b)
{
    return strcmp(((hmm_annot_t *) a)->sname, ((hmm_annot_t *) b)->sname);
}

static int hmm_annot_cmpfunc4(const void *a, const void *b)
{
    uint32_t x, y;
    x = ((hmm_annot_t *) a)->sid;
    y = ((hmm_annot_t *) b)->sid;
    return (x > y) - (x < y);
}

// annotation comparison function for ORDER_SID_OG key order: sid - og_type - gid - score
static int hmm_annot_cmpfunc5(const void *a, const void *b)
{
    hmm_annot_t *x, *y;
    x = (hmm_annot_t *) a;
    y = (hmm_annot_t *) b;

    if (x->sid != y->sid)
        return (x->sid > y->sid) - (x->sid < y->sid);
    if (x->og_type != y->og_type)
        return (x->og_type > y->og_type) - (x->og_type < y->og_type);
    if (x->gid != y->gid)
        return (x->gid > y->gid) - (x->gid < y->gid);
    return (x->score < y->score) - (x->score > y->score);
}

// annotation comparison function for ORDER_SID_CO key order: sid - alifrom - alito
static int hmm_annot_cmpfunc6(const void *a, const void *b)
{
    hmm_annot_t *x, *y;
    x = (hmm_annot_t *) a;
    y = (hmm_annot_t *) b;
    if (x->sid != y->sid)
        return (x->sid > y->sid) - (x->sid < y->sid);
    if (x->alifrom != y->alifrom)
        return (x->alifrom > y->alifrom) - (x->alifrom < y->alifrom);
    return (x->alito > y->alito) - (x->alito < y->alito);
}

static uint64_t *hmm_annot_db_sort_index_gid(hmm_annot_db_t *annot_db)
{
    // index annotation for fast access to list of a given gene
    uint32_t i, j, n_gene, gid;
    hmm_annot_t *annots;
    annots = annot_db->a;
    n_gene = annot_db->n_gene;
    uint64_t *index;
    MYCALLOC(index, n_gene);
    gid = annots[0].gid;
    for (i = 1, j = 0; i < annot_db->n; ++i) {
        if (gid != annots[i].gid) {
            index[gid] = (uint64_t) j << 32 | (i - j);
            j = i;
            gid = annots[i].gid;
        }
    }
    index[gid] = (uint64_t) j << 32 | (i - j);
    return index;
}

static uint64_t *hmm_annot_db_sort_index_sid(hmm_annot_db_t *annot_db)
{
    // index annotation for fast access to gene list of seqs
    uint32_t i, j, n_seg, sid;
    hmm_annot_t *annots;
    annots = annot_db->a;
    n_seg = annot_db->n_seg;
    uint64_t *index;
    MYCALLOC(index, n_seg);
    sid = annots[0].sid;
    for (i = 1, j = 0; i < annot_db->n; ++i) {
        if (sid != annots[i].sid) {
            index[sid] = (uint64_t) j << 32 | (i - j);
            j = i;
            sid = annots[i].sid;
        }
    }
    index[sid] = (uint64_t) j << 32 | (i - j);
    return index;
}

static void hmm_annot_db_sort_index(hmm_annot_db_t *annot_db)
{
    if (annot_db->index) free(annot_db->index);
    annot_db->index = 0;
    hmm_annot_so_t so = annot_db->so;
    if (so == ORDER_UNKNOWN || so == ORDER_UNSORTED ||
            so == ORDER_GNAME || so == ORDER_SNAME)
        return;
    switch(so) {
        case ORDER_GID:
            annot_db->index = hmm_annot_db_sort_index_gid(annot_db);
            annot_db->n_idx = annot_db->n_gene;
            break;
        case ORDER_SID:
        case ORDER_SID_OG:
        case ORDER_SID_CO:
            annot_db->index = hmm_annot_db_sort_index_sid(annot_db);
            annot_db->n_idx = annot_db->n_seg;
            break;
        default:
            fprintf(stderr, "[E::%s] undefined sorting: %d\n", __func__, so);
            exit(EXIT_FAILURE);
    }
}

void hmm_annot_db_sort(hmm_annot_db_t *annot_db, hmm_annot_so_t so)
{
    if (so == annot_db->so) return;
    switch(so) {
        case ORDER_GNAME:
            qsort(annot_db->a, annot_db->n, sizeof(hmm_annot_t), hmm_annot_cmpfunc1);
            break;
        case ORDER_GID:
            qsort(annot_db->a, annot_db->n, sizeof(hmm_annot_t), hmm_annot_cmpfunc2);
            break;
        case ORDER_SNAME:
            qsort(annot_db->a, annot_db->n, sizeof(hmm_annot_t), hmm_annot_cmpfunc3);
            break;
        case ORDER_SID:
            qsort(annot_db->a, annot_db->n, sizeof(hmm_annot_t), hmm_annot_cmpfunc4);
            break;
        case ORDER_SID_OG:
            qsort(annot_db->a, annot_db->n, sizeof(hmm_annot_t), hmm_annot_cmpfunc5);
            break;
        case ORDER_SID_CO:
            qsort(annot_db->a, annot_db->n, sizeof(hmm_annot_t), hmm_annot_cmpfunc6);
            break;
        default:
            fprintf(stderr, "[E::%s] undefined sorting: %d\n", __func__, so);
            exit(EXIT_FAILURE);
    }
    annot_db->so = so;
    hmm_annot_db_sort_index(annot_db);
}

// TODO check sort order
hmm_annot_t *hmm_annot_db_index_query(hmm_annot_db_t *annot_db, uint32_t id, uint32_t *n)
{
    hmm_annot_t *q = 0;
    *n = 0;
    if (id < annot_db->n_idx) {
        q = &annot_db->a[annot_db->index[id] >> 32];
        *n = (uint32_t) annot_db->index[id];
    }
    return q;
}

// TODO check sort order
hmm_annot_t *hmm_annot_db_sname_query(hmm_annot_db_t *annot_db, char *sname, uint32_t *n)
{
    return hmm_annot_db_index_query(annot_db, hmm_annot_sname2id(annot_db, sname), n);
}

// TODO check sort order
hmm_annot_t *hmm_annot_db_gname_query(hmm_annot_db_t *annot_db, char *gname, uint32_t *n)
{
    return hmm_annot_db_index_query(annot_db, hmm_annot_gname2id(annot_db, gname), n);
}

static void print_hmm_annot_bed6_header(FILE *fo)
{
    fprintf(fo, "#seq_name align_from align_to gene_name score_capped_at_%d strand\n", MAX_BED_SCORE);
}

void hmm_annot_formatted_print_index_list(hmm_annot_db_t *annot_db, uint64_t *index_list, size_t n, FILE *fo,
        OG_TYPE_t og_type, double max_evalue, int header)
{
    uint32_t i, j, n_annot;
    int score;
    hmm_annot_t *annot;
    hmm_annot_db_sort(annot_db, ORDER_SID_CO);
    if (header)
        print_hmm_annot_bed6_header(fo);
    for (i = 0; i < n; ++i) {
        annot = hmm_annot_db_index_query(annot_db, index_list[i], &n_annot);
        for (j = 0; j < n_annot; ++j) {
            if (annot[j].og_type != (uint32_t) og_type || 
                    annot[j].evalue > max_evalue) continue;
            score = lround(annot[j].score);
            score = MIN(score, MAX_BED_SCORE);
            fprintf(fo, "%s\t%u\t%u\t%s\t%d\t%c\n", annot[j].sname, annot[j].alifrom, 
                    annot[j].alito, annot[j].gname, score, annot[j].strand? '-' : '+');
        }
    }
}

void hmm_annot_formatted_print_sname_list(hmm_annot_db_t *annot_db, char **sname_list, size_t n, FILE *fo,
        OG_TYPE_t og_type, double max_evalue, int header)
{
    uint32_t i, j, n_annot;
    int score;
    hmm_annot_t *annot;
    hmm_annot_db_sort(annot_db, ORDER_SID_CO);
    if (header)
        print_hmm_annot_bed6_header(fo);
    for (i = 0; i < n; ++i) {
        annot = hmm_annot_db_sname_query(annot_db, sname_list[i], &n_annot);
        for (j = 0; j < n_annot; ++j) {
            if (annot[j].og_type != (uint32_t) og_type ||
                    annot[j].evalue > max_evalue) continue;
            score = lround(annot[j].score);
            score = MIN(score, MAX_BED_SCORE);
            fprintf(fo, "%s\t%u\t%u\t%s\t%d\t%c\n", annot[j].sname, annot[j].alifrom,
                    annot[j].alito, annot[j].gname, score, annot[j].strand? '-' : '+');
        }
    }
}

void hmm_annot_bed6_sname_add(hmm_annot_bed6_db_t *annots, hmm_annot_db_t *annot_db, char *cname, 
        char *sname, uint32_t len, uint32_t beg, int rev, uint32_t offset, OG_TYPE_t og_type, double max_evalue)
{
    uint32_t i, n_annot, alifrom, alito, alilen, strand;
    int score;
    hmm_annot_t *annot;
    hmm_annot_db_sort(annot_db, ORDER_SID_CO);
    annot = hmm_annot_db_sname_query(annot_db, sname, &n_annot);
    for (i = 0; i < n_annot; ++i) {
        if (annot[i].og_type != (uint32_t) og_type ||
                annot[i].evalue > max_evalue) continue;
        alifrom = annot[i].alifrom;
        alito = annot[i].alito;
        if (alifrom > alito) continue; // this should never happend unless nhmmer bug
        alilen = alito - alifrom;
        strand = annot[i].strand;
        score = lround(annot[i].score);
        score = MIN(score, MAX_BED_SCORE);
        if (rev) {
            SWAP(alifrom, alito);
            alifrom = len - alifrom;
            alito = len - alito;
            strand = !strand;
        }
        alifrom = MAX(alifrom, beg);
        alito = MAX(alito, beg);
        alifrom -= beg;
        alito -= beg;
        // for annotation clipped at the ends
        // keep it only when at least [half] the sequence retained
        if ((alito - alifrom) < (double) alilen * .5)
            continue; 
        alifrom += offset;
        alito += offset;
            
        hmm_annot_bed6_t bed6 = {cname, annot[i].gname, strand? '-' : '+', alifrom, alito, score};
        kv_push(hmm_annot_bed6_t, *annots, bed6);
    }
}

static int annot_bed6_cmpfunc(const void *a, const void *b)
{
    hmm_annot_bed6_t *x, *y;
    x = (hmm_annot_bed6_t *) a;
    y = (hmm_annot_bed6_t *) b;
    int c;
    c = strcmp(x->sname, y->sname);
    if (c) return c;
    c = (x->alifrom > y->alifrom) - (x->alifrom < y->alifrom);
    if (c) return c;
    return (x->alito > y->alito) - (x->alito < y->alito);
}

void hmm_annot_print_bed6(hmm_annot_bed6_db_t *annots, FILE *fo, int header)
{
    if (!annots->n) return;

    if (header)
        print_hmm_annot_bed6_header(fo);
    size_t i;
    hmm_annot_bed6_t *bed6;
    
    qsort(annots->a, annots->n, sizeof(hmm_annot_bed6_t), annot_bed6_cmpfunc);
    
    for (i = 0; i < annots->n; ++i) {
        bed6 = &annots->a[i];
        fprintf(fo, "%s\t%u\t%u\t%s\t%d\t%c\n", bed6->sname, bed6->alifrom,
                bed6->alito, bed6->gname, bed6->score, bed6->strand);
    }
}

int is_trn(hmm_annot_t *annot)
{
    return !strncmp(annot->gname, "trn", 3);
}

int is_rrn(hmm_annot_t *annot)
{
    return !strncmp(annot->gname, "rrn", 3);
}
