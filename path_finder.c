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
 * 03/08/22 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "ketopt.h"
#include "kvec.h"
#include "khashl.h"

#include "path.h"
#include "graph.h"
#include "hmmannot.h"

#define PATHFINDER_VERSION "0.1"

int VERBOSE = 0;
double realtime0;

char TAG_ARC_COV[4]; // arc coverage
char TAG_SEQ_COV[4]; // seq coverage
char TAG_SBP_COV[4]; // seq total base coverage

KHASHL_MAP_INIT(KH_LOCAL, kh_s32_t, kh_s32, khint32_t, uint32_t, kh_hash_dummy, kh_eq_generic);

typedef struct {
    uint8_t type; // og type
    double score; // annotation score (top 'n_core' genes for each seg)
    double sscore; // secondary annotation score
    uint32_t len; // total length of segs
    uint32_t nv; // number of segs
    uint32_t *v; // seg ids
    uint32_t ng; // number of genes (controlled by 'no_trn', 'max_eval')
    uint64_t *g; // sorted list of genes (best hit) gid << 32 | score
} og_component_t;

typedef struct {size_t n, m; og_component_t *a; } og_component_v;

static int annot_cmpfunc_s(const void *a, const void *b)
{
    hmm_annot_t *x, *y;
    x = (hmm_annot_t *) a;
    y = (hmm_annot_t *) b;

    return x->sid != y->sid? 
        (x->sid > y->sid) - (x->sid < y->sid) : 
        (x->score < y->score) - (x->score > y->score);
}

static int og_cmpfunc_r(const void *a, const void *b)
{
    return (((og_component_t *) a)->score < ((og_component_t *) b)->score) - 
        (((og_component_t *) a)->score > ((og_component_t *) b)->score);
}

static int u64_cmpfunc_r(const void *a, const void *b)
{
    return (*(uint64_t *) a < *(uint64_t *) b) - (*(uint64_t *) a > *(uint64_t *) b);
}

static int is_trn(hmm_annot_t *annot)
{
    return !strncmp(annot->gname, "trn", 3);
}

static void og_component_destroy(og_component_t *og_component)
{
    if (og_component->v) free(og_component->v);
    if (og_component->g) free(og_component->g);
}

static void og_component_v_destroy(og_component_v *component_v)
{
    size_t i;
    for (i = 0; i < component_v->n; ++i)
        og_component_destroy(&component_v->a[i]);
    if (component_v->a)
        free(component_v->a);
    free(component_v);
}

// if the graph size is larger than COMMON_MAX_PLTD_SIZE, the sequence is likely mito
// size will only include one copy of IR
static uint32_t COMMON_MAX_PLTD_SIZE = 200000;
static uint32_t COMMON_MIN_PLTD_SIZE =  80000;
// pltd to mito score fold threshold to mark graph as plat without considering other conditions
static double PLTD_TO_MITO_FST = 5.0;

og_component_v *annot_seq_og_type(hmm_annot_v *annot_v, asg_t *asg, int no_trn, 
        double max_eval, int n_core, int min_len, int min_score)
{
    uint32_t i, j, k, n_seg, n_gene, sid, gid;
    hmm_annot_t *annots, *annot;
    double m_score, p_score, *mito_score, *pltd_score;
    int m_n, p_n, len;
    uint8_t *visited, og_t;
    og_component_t *component;
    og_component_v *component_v;
    asmg_t *g;

    g = asg->asmg;
    n_gene = annot_v->n;
    annots = annot_v->a;
    n_seg = asg->n_seg;

    if (n_gene == 0) return 0;
    
    // set seg id
    for (i = 0; i < n_gene; ++i)
        annots[i].sid = asg_name2id(asg, annots[i].sname);

    // sort annotations by sid - score
    qsort(annots, n_gene, sizeof(hmm_annot_t), annot_cmpfunc_s);

    // calculate annotation score for each seg
    sid = annots[0].sid;
    m_score = p_score = .0;
    m_n = p_n = 0;
    MYCALLOC(mito_score, n_seg);
    MYCALLOC(pltd_score, n_seg);
    for (i = 0; i < n_gene; ++i) {
        annot = &annots[i];
        if (annot->sid != sid) {
            mito_score[sid] = m_score;
            pltd_score[sid] = p_score;
            sid = annot->sid;
            m_score = p_score = .0;
            m_n = p_n = 0;
        }
        if (annot->evalue > max_eval || (no_trn && is_trn(annot)))
            continue;
        if (annot->mito && m_n++ < n_core)
            m_score += annot->score;
        if (annot->pltd && p_n++ < n_core)
            p_score += annot->score;
    }
    mito_score[sid] = m_score;
    pltd_score[sid] = p_score;

    // extract each subgraph and do classification
    // TODO: deal with graphs when mito and pltd are in the same subgraph
    MYCALLOC(component_v, 1);
    MYCALLOC(visited, n_seg);
    for (i = 0; i < n_seg; ++i) {
        if (visited[i]) continue;
        asg_subgraph(asg, &i, 1, 0);
        m_score = p_score = .0;
        m_n = len = 0;
        for (j = 0; j < n_seg; ++j) {
            if (!g->vtx[j].del) {
                m_score += mito_score[j];
                p_score += pltd_score[j];
                ++m_n;
                len += g->vtx[j].len;
                visited[j] = 1;
            }
        }
        og_t = 0;
        if (len < min_len && m_score < min_score && p_score < min_score)
            // mark as non-organelle sequence
            og_t = OG_NONE;
        else
            // 1, mito; 2, pltd; 3, unclassified
            og_t = m_score == p_score? OG_UNCLASSIFIED : (m_score > p_score? OG_MITO : OG_PLTD);

        if (og_t == OG_MITO || og_t == OG_PLTD) {
            // make seg list
            uint32_t *comp_v;
            MYMALLOC(comp_v, m_n);
            k = 0;
            for (j = 0; j < n_seg; ++j)
                if (!g->vtx[j].del)
                    comp_v[k++] = j;
            // make gene list
            // get all genes in the component
            kvec_t(uint64_t) comp_g;
            kv_init(comp_g);
            for (j = 0; j < n_gene; ++j) {
                annot = &annots[j];
                if (g->vtx[annot->sid].del || 
                        annot->evalue > max_eval || 
                        (no_trn && is_trn(annot)))
                    continue;
                // gid:30 | OG_TYPE:2 | score:32
                kv_push(uint64_t, comp_g, 
                        ((uint64_t) hmm_annot_name2id(annot_v, annot->gname) << 2 | 
                         (annot->mito? OG_MITO : OG_PLTD)) << 32 | (uint32_t) annot->score);
            }
            // sort by gid and then score - descending order
            qsort(comp_g.a, comp_g.n, sizeof(uint64_t), u64_cmpfunc_r);
            // keep only the best hit for each gene
            k = 0;
            gid = UINT32_MAX;
            for (j = 0; j < comp_g.n; ++j) {
                if (comp_g.a[j] >> 32 != gid) {
                    comp_g.a[k++] = comp_g.a[j];
                    gid = comp_g.a[j] >> 32;
                }
            }
            comp_g.n = comp_g.m = k;
            MYREALLOC(comp_g.a, k);

            // add og_component
            kv_pushp(og_component_t, *component_v, &component);
            component->type = og_t;
            component->score = m_score > p_score? m_score : p_score;
            component->sscore = m_score < p_score? m_score : p_score;        
            component->len = len;
            component->nv = m_n;
            component->v = comp_v;
            component->ng = comp_g.n;
            component->g = comp_g.a;
        }

        if (VERBOSE > 0)
            fprintf(stderr, "[M::%s] subgraph seeding from %s: segs, %d; size, %d; mito score, %.3f; pltd score, %.3f; classification, %d\n", 
                    __func__, asg->seg[i].name, m_n, len, m_score, p_score, og_t);
    }
    free(visited);
    free(mito_score);
    free(pltd_score);

    // TODO fix this by maybe improving the specificity of the HMM database
    // revisit to deal with mito misclassification
    // find the best pltd component
    uint32_t p_b, p_b1;
    double p_s, p_s1;
    p_b = p_b1 = UINT32_MAX;
    p_s = p_s1 = 0;
    for (i = 0; i < component_v->n; ++i) {
        component = &component_v->a[i];
        if (component->type != OG_PLTD)
            continue;
        if (component->score > p_s && component->len >= COMMON_MIN_PLTD_SIZE) {
            if (component->len <= COMMON_MAX_PLTD_SIZE) {
                p_b = i;
                p_s = component->score;
            }
            p_b1 = i;
            p_s1 = component->score;
        }
    }
    if (p_b == UINT32_MAX) p_b = p_b1;
    if (p_b != UINT32_MAX) {
        // PLTD is likely there
        for (i = 0; i < component_v->n; ++i) {
            if (i == p_b) continue;
            component = &component_v->a[i];
            if (component->type != OG_PLTD)
                continue;
            if (component->score < component->sscore * PLTD_TO_MITO_FST || 
                    component->len > COMMON_MAX_PLTD_SIZE) {
                component->type = OG_MITO;
                double tmp = component->score;
                component->score = component->sscore;
                component->sscore = tmp;
            }
        }
    }

    qsort(component_v->a, component_v->n, sizeof(og_component_t), og_cmpfunc_r);

    return component_v;
}

static void print_og_classification_summary(asg_t *asg, hmm_annot_v *annot_v, og_component_v *og_components, FILE *fo)
{
    uint32_t i, j;
    og_component_t *component;
    for (i = 0; i < og_components->n; ++i) {
        component = &og_components->a[i];
        fprintf(fo, "[M::%s] OG component %u \n", __func__, i);
        fprintf(fo, "[M::%s] OG component %u og_type: %s\n", __func__, i, component->type==OG_MITO? "mito" : "pltd");
        fprintf(fo, "[M::%s] OG component %u og_score: %.1f\n", __func__, i, component->score);
        fprintf(fo, "[M::%s] OG component %u og_sscore: %.1f\n", __func__, i, component->sscore);
        fprintf(fo, "[M::%s] OG component %u og_len: %u\n", __func__, i, component->len);
        fprintf(fo, "[M::%s] OG component %u og_nv: %u\n", __func__, i, component->nv);
        fprintf(fo, "[M::%s] OG component %u og_v: %s", __func__, i, asg->seg[component->v[0]].name);
        for (j = 1; j < component->nv; ++j)
            fprintf(fo, " %s", asg->seg[component->v[j]].name);
        fprintf(fo, "\n");
        fprintf(fo, "[M::%s] OG component %u og_ng: %u\n", __func__, i, component->ng);
        for (j = 0; j < component->ng; ++j)
            fprintf(fo, "[M::%s] OG component %u og_g: %s %u\n", __func__, i,
                    annot_v->dict[component->g[j]>>34], (uint32_t) component->g[j]);
    }
}

static int clean_graph_by_sequence_coverage(asg_t *asg, double min_cf, int max_copy)
{
    uint32_t i, j, n_seg, nv;
    uint8_t *visited;
    double avg_cov;
    kvec_t(uint32_t) rm_v;
    asmg_t *g;

    g = asg->asmg;
    n_seg = asg->n_seg;
    kv_init(rm_v);

    MYCALLOC(visited, n_seg);
    for (i = 0; i < n_seg; ++i) {
        if (visited[i]) continue;
        asg_subgraph(asg, &i, 1, 0);
        avg_cov = estimate_sequence_copy_number_from_coverage(asg, 0, max_copy);
        for (j = 0; j < n_seg; ++j) {
            if (g->vtx[j].del) continue;
            if (avg_cov && asg->seg[j].cov / avg_cov < min_cf)
                kv_push(uint32_t, rm_v, j);
            visited[j] = 1;
        }
    }
    
    nv = rm_v.n;
    // restore all vtx and arc
    // except for arcs linked to vtx to clean
    for (i = 0; i < g->n_arc; ++i)
        g->arc[i].del = 0;
    for (i = 0; i < nv; ++i)
        asmg_vtx_del(g, rm_v.a[i], 1);
    for (i = 0; i < n_seg; ++i)
        g->vtx[i].del = 0;
    asmg_finalize(g, 1);

    kv_destroy(rm_v);
    free(visited);
    
    if (VERBOSE > 1) {
        if (VERBOSE > 2) {
            fprintf(stderr, "[M::%s] graph after cleaning\n", __func__);
            asg_print(asg, stderr, 1);
        }
        fprintf(stderr, "[M::%s] number sequence cleaned: %u\n", __func__, nv);
        fprintf(stderr, "[M::%s] graph stats after cleaning\n", __func__);
        asg_stat(asg, stderr);
    }

    return nv;
}

static double sequence_covered_by_path(asg_t *asg, path_t *path, uint32_t len)
{
    uint32_t i, n, l, *v;
    uint8_t *flag;
    MYCALLOC(flag, asg->n_seg);
    v = path->v;
    l = 0;
    for (i = 0, n = path->nv; i < n; ++i) {
        if (!flag[v[i]>>1]) {
            l += asg->seg[v[i]>>1].len;
            flag[v[i]>>1] = 1;
        }
    }
    free(flag);
    return (double) l / len;
}

static void parse_organelle_component(asg_t *asg, hmm_annot_v *annot_v, og_component_v *og_components,
        int max_copy, int min_ex_g, uint32_t max_d_len, double seq_cf, int do_clean, 
        char *out_pref, int out_opt, uint8_t og_type)
{
    assert(og_type == OG_MITO || og_type == OG_PLTD);

    FILE *out_stats, *out_seq, *out_utg, *out_gfa;
    char *m_str;
    MYMALLOC(m_str, strlen(out_pref) + 64);
    sprintf(m_str, "%s.%s.stats", out_pref, og_type==OG_MITO? "mito" : "pltd");
    out_stats = fopen(m_str, "w");
    if (!out_stats) {
        fprintf(stderr, "[E::%s] failed to open file %s to write\n", __func__, m_str);
        exit(EXIT_FAILURE);
    }
    sprintf(m_str, "%s.%s.fasta", out_pref, og_type==OG_MITO? "mito" : "pltd");
    out_seq = fopen(m_str, "w");
    if (!out_seq) {
        fprintf(stderr, "[E::%s] failed to open file %s to write\n", __func__, m_str);
        exit(EXIT_FAILURE);
    }
    sprintf(m_str, "%s.%s.unassembled.fasta", out_pref, og_type==OG_MITO? "mito" : "pltd");
    out_utg = fopen(m_str, "w");
    if (!out_utg) {
        fprintf(stderr, "[E::%s] failed to open file %s to write\n", __func__, m_str);
        exit(EXIT_FAILURE);
    }
    sprintf(m_str, "%s.%s.gfa", out_pref, og_type==OG_MITO? "mito" : "pltd");
    out_gfa = fopen(m_str, "w");
    if (!out_gfa) {
        fprintf(stderr, "[E::%s] failed to open file %s to write\n", __func__, m_str);
        exit(EXIT_FAILURE);
    }
    free(m_str);

    uint32_t i, j, opt_circ;
    int c, ex_g, absent;
    uint64_t n_seg;
    og_component_t *component;
    kh_s32_t *h_genes;
    khint32_t k;
    kvec_t(uint32_t) sub_v;

    n_seg = asg->n_seg;
    h_genes = kh_s32_init();
    kv_init(sub_v);
    c = 0;
    opt_circ = 0;
    for (i = 0; i < og_components->n; ++i) {
        component = &og_components->a[i];
        if (component->type != og_type)
            continue;
        if (VERBOSE > 0)
            fprintf(stderr, "[M::%s] processing subgraph seeding from %s: "
                    "type, %s; score, %.1f; sscore, %.1f; len, %u; nv, %u; ng, %u\n", 
                    __func__, asg->seg[component->v[0]].name, component->type==OG_MITO? "mito" : "pltd", 
                    component->score, component->sscore, component->len, component->nv, component->ng);
        
        // find number of extra genes would be added
        ex_g = 0;
        for (j = 0; j < component->ng; ++j) {
            if ((component->g[j]>>32 & 0x3) != og_type)
                continue;
            k = kh_s32_get(h_genes, component->g[j]>>32);
            if (k == kh_end(h_genes) || kh_val(h_genes, k) < (uint32_t) component->g[j])
                ++ex_g;
        }
        
        // skip component with no enough extra genes
        if (ex_g < min_ex_g) {
            if (VERBOSE > 0)
                fprintf(stderr, "[M::%s] subgraph seeding from %s SKIPPED due to insufficient gene gain (%d)\n",
                        __func__, asg->seg[component->v[0]].name, ex_g);
            continue;
        }

        // processing the subgraph
        asg_subgraph(asg, component->v, 1, 0);
        
        if (VERBOSE > 1) {
            if (VERBOSE > 2) {
                fprintf(stderr, "[M::%s] subgraph created\n", __func__);
                asg_print(asg, stderr, 1);
            }    
            fprintf(stderr, "[M::%s] subgraph stats\n", __func__);
            asg_stat(asg, stderr);
        }
        
        if (do_clean) {
            asmg_drop_tip(asg->asmg, 10, 10000, 0, VERBOSE);
            if (VERBOSE > 1) {
                if (VERBOSE > 2) {
                    fprintf(stderr, "[M::%s] subgraph tip-dropped\n", __func__);
                    asg_print(asg, stderr, 1);
                }
                fprintf(stderr, "[M::%s] tip-dropped subgraph stats\n", __func__);
                asg_stat(asg, stderr);
            }
            
            asmg_pop_bubble(asg->asmg, 100000, 0, 0, 10000, 0, VERBOSE);
            if (VERBOSE > 1) {
                if (VERBOSE > 2) {
                    fprintf(stderr, "[M::%s] subgraph bubble-popped\n", __func__);
                    asg_print(asg, stderr, 1);
                }
                fprintf(stderr, "[M::%s] bubble-popped subgraph stats\n", __func__);
                asg_stat(asg, stderr);
            }
            // TODO filter low coverage arcs
        }
        
        if (asmg_vtx_n1(asg->asmg) == 0)
            continue;
        
        int *copy_number;
        MYMALLOC(copy_number, n_seg);
        estimate_sequence_copy_number_from_coverage(asg, copy_number, max_copy);
        // make a copy of the graph
        asg_t *asg_copy = asg_make_copy(asg);
        kh_u32_t *seg_dups = sequence_duplication_by_copy_number(asg_copy, copy_number);
        path_v paths;
        kv_init(paths);
        graph_path_finder(asg_copy, seg_dups, &paths);
        kh_u32_destroy(seg_dups);
        asg_destroy(asg_copy);

        if (paths.n == 0) {
            // not possible to solve the graph
            // likely because the graph is too complicated
            if (VERBOSE > 0)
                fprintf(stderr, "[M::%s] subgraph seeding from %s is unresolvable, output unitigs as unassembled\n",
                        __func__, asg->seg[component->v[0]].name);
            // FIXME sequences removed by graph clean are not included
            asg_print_fa(asg, out_utg, 60);
            kv_push(uint32_t, sub_v, component->v[0]);
        } else {
            uint32_t b, is_circ;
            double f;
            if (og_type == OG_PLTD) // only for pltd
                for (j = 0; j < paths.n; ++j)
                    path_rotate(asg, &paths.a[j], annot_v, 2);
            path_sort(&paths);
            // TODO improve the selection criteria
            b = select_best_seq(asg, &paths, 0, out_opt, seq_cf, 0, og_type == OG_PLTD);
            f = sequence_covered_by_path(asg, &paths.a[b], component->len);
            is_circ = paths.a[b].circ;
            if (VERBOSE > 0)
                fprintf(stderr, "[M::%s] best path after first pass: type, %s; coverage, %.3f\n", 
                        __func__, is_circ? "circular" : "linear", f);
            if (!is_circ || f < seq_cf) {
                // adjust the sequence copy number and redo path finding
                int updated = adjust_sequence_copy_number_by_graph_layout(asg, copy_number, max_copy);
                if (updated) {
                    // make another copy of the graph
                    asg_t *asg_copy1 = asg_make_copy(asg);
                    // do path finding
                    kh_u32_t *seg_dups1 = sequence_duplication_by_copy_number(asg_copy1, copy_number);
                    path_v paths1;
                    kv_init(paths1);
                    graph_path_finder(asg_copy1, seg_dups1, &paths1);
                    kh_u32_destroy(seg_dups1);
                    asg_destroy(asg_copy1);

                    uint32_t b1, is_circ1;
                    double f1;
                    if (og_type == OG_PLTD) // only for pltd
                        for (j = 0; j < paths1.n; ++j)
                            path_rotate(asg, &paths1.a[j], annot_v, 2);
                    path_sort(&paths1);
                    // TODO improve the selection criteria
                    b1 = select_best_seq(asg, &paths1, 0, out_opt, seq_cf, 0, og_type == OG_PLTD);
                    f1 = sequence_covered_by_path(asg, &paths1.a[b1], component->len);
                    is_circ1 = paths1.a[b1].circ;

                    // compare new path to old to make final choice
                    if ((is_circ1 && f1 >= seq_cf) || ((is_circ1 || !is_circ) && f1 > f)) {
                        // choose the new one
                        b = b1;
                        f = f1;
                        is_circ = is_circ1;
                        for (j = 0; j < paths.n; ++j)
                            path_destroy(&paths.a[j]);
                        kv_destroy(paths);
                        paths = paths1;
                    } else {
                        for (j = 0; j < paths1.n; ++j)
                            path_destroy(&paths1.a[j]);
                        kv_destroy(paths1);
                    }
                    if (VERBOSE > 0)
                        fprintf(stderr, "[M::%s] best path after second pass: type, %s; coverage, %.3f\n",
                                __func__, is_circ? "circular" : "linear", f);

                }
            }

            if (is_circ || !opt_circ || component->len > max_d_len) {
                if (!opt_circ) opt_circ = is_circ;
                
                ++c;
                kv_push(uint32_t, sub_v, component->v[0]);

                // update gene table
                for (j = 0; j < component->ng; ++j) {
                    if ((component->g[j]>>32 & 0x3) != og_type)
                        continue;
                    k = kh_s32_get(h_genes, component->g[j]>>32);
                    if (k == kh_end(h_genes) || kh_val(h_genes, k) < (uint32_t) component->g[j]) {
                        k = kh_s32_put(h_genes, component->g[j]>>32, &absent);
                        kh_val(h_genes, k) = (uint32_t) component->g[j];
                    }
                }

                path_stats(asg, &paths, out_stats);
                print_seq(asg, &paths.a[b], out_seq, c, 0, 60);
            
                if (VERBOSE > 0)
                    fprintf(stderr, "[M::%s] processing subgraph seeding from %s DONE, %d better genes gained\n",
                            __func__, asg->seg[component->v[0]].name, ex_g);
            }
        }

        for (j = 0; j < paths.n; ++j)
            path_destroy(&paths.a[j]);
        kv_destroy(paths);

        free(copy_number);
    }

    // extract and output organelle gfa components
    if (sub_v.n > 0) {
        // in subgraph all vtx and arcs deleted during expansion are back
        asg_subgraph(asg, sub_v.a, sub_v.n, 0);
        asg_print(asg, out_gfa, 0);
    }
    kv_destroy(sub_v);

    fclose(out_stats);
    fclose(out_seq);
    fclose(out_utg);
    fclose(out_gfa);
}

static ko_longopt_t long_options[] = {
    { "longest",        ko_no_argument,       301 },
    { "circular",       ko_no_argument,       302 },
    { "all",            ko_no_argument,       303 },
    { "edge-c-tag",     ko_required_argument, 304 },
    { "kmer-c-tag",     ko_required_argument, 305 },
    { "seq-c-tag",      ko_required_argument, 306 },
    { "include-trn",    ko_no_argument,       307 },
    { "max-eval",       ko_required_argument, 308 },
    { "min-s-length",   ko_required_argument, 309 },
    { "max-d-length",   ko_required_argument, 310 },
    { "no-graph-clean", ko_no_argument,       311 },
    { "mito-annot",     ko_required_argument, 'm' },
    { "pltd-annot",     ko_required_argument, 'p' },
    { "core-gene",      ko_required_argument, 'g' },
    { "min-score",      ko_required_argument, 's' },
    { "min-b-gene",     ko_required_argument, 'a' },
    { "min-seq-cov",    ko_required_argument, 'q' },
    { "max-copy",       ko_required_argument, 'c' },
    { "version",        ko_no_argument,       'V' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    const char *opt_str = "a:m:p:f:g:q:s:c:o:Vv:h";
    ketopt_t opt = KETOPT_INIT;
    int c, out_s, out_c, max_copy, ret = 0;
    asg_t *asg;
    hmm_annot_v *annot_v;
    og_component_v *og_components;
    FILE *fp_help;
    char *out_pref, *mito_annot, *pltd_annot, *ec_tag, *kc_tag, *sc_tag;
    int no_trn, n_core, min_len, min_ex_g, max_d_len, do_graph_clean;
    double max_eval, min_score, min_cf, seq_cf;

    sys_init();

    out_s = -1;
    out_c = 0;
    out_pref = 0;
    max_copy = 10;
    fp_help = stderr;
    mito_annot = 0;
    pltd_annot = 0;
    ec_tag = 0;
    kc_tag = 0;
    sc_tag = 0;
    annot_v = 0;
    no_trn = 1;
    do_graph_clean = 1;
    n_core = 10;
    max_eval = 1e-12;
    min_len = 10000;
    min_score = 100;
    min_ex_g = 3;
    max_d_len = 50000;
    seq_cf = .95;
    min_cf = .20;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >=0 ) {
        if (c == 'm') mito_annot = opt.arg;
        else if (c == 'p') pltd_annot = opt.arg;
        else if (c == 'f') seq_cf = atof(opt.arg);
        else if (c == 'g') n_core = atoi(opt.arg);
        else if (c == 's') min_score = atof(opt.arg);
        else if (c == 'c') max_copy = atoi(opt.arg);
        else if (c == 'a') min_ex_g = atoi(opt.arg);
        else if (c == 'q') min_cf = atof(opt.arg);
        else if (c == 'o') out_pref = opt.arg;
        else if (c == 301) out_s = 0, ++out_c;
        else if (c == 302) out_s = 1, ++out_c;
        else if (c == 303) out_s = 2, ++out_c;
        else if (c == 304) ec_tag = opt.arg;
        else if (c == 305) kc_tag = opt.arg;
        else if (c == 306) sc_tag = opt.arg;
        else if (c == 307) no_trn = 0;
        else if (c == 308) max_eval = atof(opt.arg);
        else if (c == 309) min_len = atoi(opt.arg);
        else if (c == 310) max_d_len = atoi(opt.arg);
        else if (c == 311) do_graph_clean = 0;
        else if (c == 'v') VERBOSE = atoi(opt.arg);
        else if (c == 'h') fp_help = stdout;
        else if (c == 'V') {
            puts(PATHFINDER_VERSION);
            return 0;
        }
        else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
        else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (argc == opt.ind || fp_help == stdout) {
        fprintf(fp_help, "Usage: path_finder [options] <file>[.gfa[.gz]] [path_str]\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "  Input/Output:\n");
        fprintf(fp_help, "    -m FILE              mitochondria core gene annotation file [NULL]\n");
        fprintf(fp_help, "    -p FILE              plastid core gene annotation file [NULL]\n");
        fprintf(fp_help, "    -o FILE              prefix of output files [oatk.asm]\n");
        fprintf(fp_help, "    -f FLOAT             prefer circular path to longest if >= FLOAT sequence covered [%.3f]\n", seq_cf);
        fprintf(fp_help, "    --longest            output only the longest path [default]\n");
        fprintf(fp_help, "    --circular           output only the longest circular path\n");
        fprintf(fp_help, "    --all                output all best paths\n");
        fprintf(fp_help, "    --edge-c-tag STR     edge coverage tag in the GFA file [EC:i] \n");
        fprintf(fp_help, "    --kmer-c-tag STR     kmer coverage tag in the GFA file [KC:i] \n");
        fprintf(fp_help, "    --seq-c-tag STR      sequence coverage tag in the GFA file [SC:f]\n");
        fprintf(fp_help, "    -v INT               verbose level [%d]\n", VERBOSE);
        fprintf(fp_help, "    --version            show version number\n");
        fprintf(fp_help, "  Classification:\n");
        fprintf(fp_help, "    -g INT               number of top core gene annotations to consider [%d]\n", n_core);
        fprintf(fp_help, "    -s FLOAT             minimum annotation score of a subgraph [%.1f]\n", min_score);
        fprintf(fp_help, "    -a INT               minimum number of addtional core genes to include a sequence [%d]\n", min_ex_g);
        fprintf(fp_help, "    --include-trn        include TRN type genes\n");
        fprintf(fp_help, "    --max-eval  FLOAT    maximum E-value of a core gene [%.3e]\n", max_eval);
        fprintf(fp_help, "    --min-s-length INT   minimum length of a singleton sequence to keep [%d]\n", min_len);
        fprintf(fp_help, "    --max-d-length INT   maximum length of a singleton sequence to delete [%d]\n", max_d_len);
        fprintf(fp_help, "  Path-finding:\n");
        fprintf(fp_help, "    -c INT               maximum copy number to consider [%d]\n", max_copy);
        fprintf(fp_help, "    -q FLOAT             minimum coverage of a sequence compared to the subgraph average [%.3f]\n", min_cf);
        fprintf(fp_help, "    --no-graph-clean     do not do assembly graph clean\n");
        return fp_help == stdout? 0 : 1;
    }

    if (argc - opt.ind < 1) {
        fprintf(stderr, "[E::%s] missing input: please specify the GFA file\n", __func__);
        return 1;
    }

    if (out_c > 1) {
        fprintf(stderr, "[E::%s] options --longest, --circular and --all are mutually exclusive\n", __func__);
        return 1;
    }
    
    if (ec_tag != 0) {
        if(is_valid_gfa_tag(ec_tag)) {
            memcpy(TAG_ARC_COV, ec_tag, 4);
        } else {
            fprintf(stderr, "[E::%s] invalid GFA tag (Regexp: [A-Za-z][A-Za-z0-9]:[A|i|f|Z|B]): %s\n", __func__, ec_tag);
            return 1;
        }
    }

    if (kc_tag != 0) {
        if(is_valid_gfa_tag(kc_tag)) {
            memcpy(TAG_SBP_COV, kc_tag, 4);
        } else {
            fprintf(stderr, "[E::%s] invalid GFA tag (Regexp: [A-Za-z][A-Za-z0-9]:[A|i|f|Z|B]): %s\n", __func__, kc_tag);
            return 1;
        }
    }

    if (sc_tag != 0) {
        if(is_valid_gfa_tag(sc_tag)) {
            memcpy(TAG_SEQ_COV, sc_tag, 4);
        } else {
            fprintf(stderr, "[E::%s] invalid GFA tag (Regexp: [A-Za-z][A-Za-z0-9]:[A|i|f|Z|B]): %s\n", __func__, sc_tag);
            return 1;
        }
    }

    if (out_s < 0) out_s = 0;
    
    if (!out_pref) out_pref = "oatk.asm";
    
    if (mito_annot) annot_v = hmm_annot_read(mito_annot, annot_v, OG_MITO);
    if (pltd_annot) annot_v = hmm_annot_read(pltd_annot, annot_v, OG_PLTD);
    if (annot_v) hmm_annot_index(annot_v);
    if (VERBOSE > 2) hmm_annot_print(annot_v->a, annot_v->n, stderr);

    asg = asg_read(argv[opt.ind]);
    if (asg == 0) {
        fprintf(stderr, "[E::%s] failed to read the graph: %s\n", __func__, argv[opt.ind]);
        exit(EXIT_FAILURE);
    }

    if (VERBOSE > 1) {
        if (VERBOSE > 2) {
            fprintf(stderr, "[M::%s] graph loaded\n", __func__);
            asg_print(asg, stderr, 1);
        }
        fprintf(stderr, "[M::%s] graph stats\n", __func__);
        asg_stat(asg, stderr);
    }
    
    if (do_graph_clean && min_cf > .0) clean_graph_by_sequence_coverage(asg, min_cf, max_copy);

    og_components = annot_seq_og_type(annot_v, asg, no_trn, max_eval, n_core, min_len, min_score);
    if (VERBOSE > 1) print_og_classification_summary(asg, annot_v, og_components, stderr);

    // graph will be changed with extra copies of sequences added
    if (mito_annot) parse_organelle_component(asg, annot_v, og_components, max_copy, min_ex_g, max_d_len, seq_cf, do_graph_clean, out_pref, out_s, OG_MITO);
    if (pltd_annot) parse_organelle_component(asg, annot_v, og_components, max_copy, min_ex_g, max_d_len, seq_cf, do_graph_clean, out_pref, out_s, OG_PLTD);
    
    if (ret) {
        fprintf(stderr, "[E::%s] failed to analysis the GFA file\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (fflush(stdout) == EOF) {
        perror("[E::%s] failed to write the results");
        exit(EXIT_FAILURE);
    }

    asg_destroy(asg);
    hmm_annot_v_destroy(annot_v);
    og_component_v_destroy(og_components);

    if (VERBOSE >= 0) {
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, PATHFINDER_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        int i;
        for (i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    return 0;
}
