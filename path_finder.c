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

#include "kvec.h"
#include "khashl.h"
#include "kthread.h"

#include "path.h"
#include "graph.h"
#include "syncasm.h"
#include "hmmannot.h"

#define PATHFINDER_VERSION "0.1"

#undef DEBUG_MINICIRCLE_REPEAT_UNIT

char TAG_ARC_COV[4]; // arc coverage
char TAG_SEQ_COV[4]; // seq coverage
char TAG_SBP_COV[4]; // seq total base coverage

KHASHL_MAP_INIT(KH_LOCAL, kh_s32_t, kh_s32, khint32_t, uint32_t, kh_hash_dummy, kh_eq_generic);

static void parse_organelle_component(asg_t *asg, hmm_annot_v *annot_v, og_component_v *og_components,
        int max_copy, int min_ex_g, uint32_t max_d_len, double seq_cf, int do_clean, int bubble_size, 
        int tip_size, double weak_cross, char *out_pref, int out_opt, uint8_t og_type, int VERBOSE)
{
    assert(og_type == OG_MITO || og_type == OG_PLTD || og_type == OG_MINI);

    FILE *out_stats, *out_seq, *out_utg, *out_gfa;
    char *m_str;
    MYMALLOC(m_str, strlen(out_pref) + 64);
    sprintf(m_str, "%s.%s.stats", out_pref, OG_TYPES[og_type]);
    out_stats = fopen(m_str, "w");
    if (!out_stats) {
        fprintf(stderr, "[E::%s] failed to open file %s to write\n", __func__, m_str);
        exit(EXIT_FAILURE);
    }
    sprintf(m_str, "%s.%s.fasta", out_pref, OG_TYPES[og_type]);
    out_seq = fopen(m_str, "w");
    if (!out_seq) {
        fprintf(stderr, "[E::%s] failed to open file %s to write\n", __func__, m_str);
        exit(EXIT_FAILURE);
    }
    sprintf(m_str, "%s.%s.unassembled.fasta", out_pref, OG_TYPES[og_type]);
    out_utg = fopen(m_str, "w");
    if (!out_utg) {
        fprintf(stderr, "[E::%s] failed to open file %s to write\n", __func__, m_str);
        exit(EXIT_FAILURE);
    }
    sprintf(m_str, "%s.%s.gfa", out_pref, OG_TYPES[og_type]);
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
                    __func__, asg->seg[component->v[0]].name, OG_TYPES[component->type], 
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
        asg_subgraph(asg, component->v, 1, 0, 0);
        
        if (VERBOSE > 1) {
            if (VERBOSE > 2) {
                fprintf(stderr, "[M::%s] subgraph created\n", __func__);
                asg_print(asg, stderr, 1);
            }    
            fprintf(stderr, "[M::%s] subgraph stats\n", __func__);
            asg_stat(asg, stderr);
        }

        if (do_clean) {
            asmg_drop_tip(asg->asmg, 10, tip_size, 0, VERBOSE);
            if (VERBOSE > 1) {
                if (VERBOSE > 2) {
                    fprintf(stderr, "[M::%s] subgraph tip-dropped\n", __func__);
                    asg_print(asg, stderr, 1);
                }
                fprintf(stderr, "[M::%s] tip-dropped subgraph stats\n", __func__);
                asg_stat(asg, stderr);
            }

            asmg_pop_bubble(asg->asmg, bubble_size, 0, 0, 1, 0, VERBOSE);
            if (VERBOSE > 1) {
                if (VERBOSE > 2) {
                    fprintf(stderr, "[M::%s] subgraph bubble-popped\n", __func__);
                    asg_print(asg, stderr, 1);
                }
                fprintf(stderr, "[M::%s] bubble-popped subgraph stats\n", __func__);
                asg_stat(asg, stderr);
            }

            asmg_remove_weak_crosslink(asg->asmg, weak_cross, 0, VERBOSE);
            if (VERBOSE > 1) {
                if (VERBOSE > 2) {
                    fprintf(stderr, "[M::%s] subgraph weak cross-link removed\n", __func__);
                    asg_print(asg, stderr, 1);
                }
                fprintf(stderr, "[M::%s] weak cross-link removed subgraph stats\n", __func__);
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
        asg_subgraph(asg, sub_v.a, sub_v.n, 0, 0);
        asg_print(asg, out_gfa, 0);
    }
    kv_destroy(sub_v);

    fclose(out_stats);
    fclose(out_seq);
    fclose(out_utg);
    fclose(out_gfa);
}

typedef struct {
    scg_ra_v *ra;
    uint64_t sid;
    uint64_t *mcs; // beg<<33|end<<1|rev end is inclusive
} mc_shared_t;

static void minicircle_analysis_thread(void *_data, long i, int tid) // kt_for() callback
{
    uint32_t j, k, beg, end, rev, nfrg;
    uint64_t sid, uid;
    scg_ra_t *ra;
    ra_frg_t *frgs;
    
    ra = &(((mc_shared_t *)_data)->ra->a[i]);
    sid = ((mc_shared_t *)_data)->sid;
    nfrg = ra->n;
    frgs = ra->a;

    if (nfrg < 2) {
        ((mc_shared_t *)_data)->mcs[i] = UINT64_MAX;
        return;
    }

    // find begin and end position
    beg = end = rev = UINT32_MAX;
    for (j = 0; j < nfrg; ++j) {
        uid = frgs[j].uid;
        if (uid>>1 != sid)
            continue;
        if (beg == UINT32_MAX)
            beg = j;
        else if (end == UINT32_MAX)
            end = j - 1;
        if (rev == UINT32_MAX) {
            rev = uid & 1;
        } else if (rev != (uid & 1)) {
            // inconsistent alignment orientation found
            rev = UINT32_MAX;
            break;
        }
    }

    if (beg == UINT32_MAX || end == UINT32_MAX || rev == UINT32_MAX) {
        // either no repeat unit or inconsistent orientation
        ((mc_shared_t *)_data)->mcs[i] = UINT64_MAX;
        return;
    }
    
    // check repeat unit consistency across the entire alignment path
    int valid = 1;
    if (beg > 0 || end < nfrg - 2) {
        uint32_t r = end - beg;
        if (beg > r) {
            valid = 0;
        } else {
            k = r - beg;
            if (++k > r) k = 0;
            for (j = 0; j < nfrg; ++j) {
                if (frgs[j].uid != frgs[beg+k].uid) {
                    valid = 0;
                    break;
                }
                if (++k > r) k = 0;
            }
        }
    }

    ((mc_shared_t *)_data)->mcs[i] = valid? ((uint64_t) beg<<33 | (uint64_t) end<<1 | rev) : UINT64_MAX;
}

static int extract_minicircles_with_anchor(scg_ra_v *ra, uint64_t anchor_sid, int n_threads)
{
    uint64_t i, j, nr, *mcs;
    
    nr = ra->n;
    MYMALLOC(mcs, nr);
    mc_shared_t mc_shared = {ra, anchor_sid, mcs};
    kt_for(n_threads, minicircle_analysis_thread, &mc_shared, nr);

#ifdef DEBUG_MINICIRCLE_REPEAT_UNIT
    scg_rv_print(ra, stderr);
    for (i = 0; i < nr; ++i) {
        if (mcs[i] == UINT64_MAX) continue;
        scg_ra_print(&ra->a[i], stderr);
        fprintf(stderr, "[DEBUG_MINICIRCLE_REPEAT_UNIT::%s] repeat unit: beg %lu; end %u; %c\n", __func__, 
                mcs[i]>>33, (uint32_t) (mcs[i]>>1), "+-"[mcs[i]&1]);
    }
#endif

    free(mcs);

    return 0;
}

static int parse_organelle_minicircle(asg_t *asg, hmm_annot_v *annot_v, og_component_v *og_components,
        double *seg_annot_score, scg_meta_t *scg_meta, int n_threads, int VERBOSE)
{
    if (og_components->n == 0) {
        if (VERBOSE > 0)
            fprintf(stderr, "[M::%s] no OG component found\n", __func__);
        return 1;
    }

    uint32_t i, anchor_sid;
    uint64_t n_seg;
    double max_s;
    og_component_t *component;
    asmg_t *asmg;

    n_seg = asg->n_seg;
    asmg = scg_meta->scg->utg_asmg;

    // find the anchor sequence in the first og component
    component = &og_components->a[0];
    if (component->type != OG_MINI)
        return 1;
    
    max_s = .0;
    anchor_sid = 0;
    for (i = 0; i < component->nv; ++i) {
        uint32_t sid = component->v[i];
        double s = seg_annot_score[sid*4+OG_MINI];
        if (s > max_s) {
            max_s = s;
            anchor_sid = sid;
        }
    }

    if (VERBOSE > 0)
        fprintf(stderr, "[M::%s] anchor sequence found: %s [len %u; score, %.3f]\n", __func__,
                asg->seg[anchor_sid].name, asg->seg[anchor_sid].len, max_s);

    assert(n_seg == asmg->n_vtx);
    for (i = 0; i < n_seg; ++i)
        assert(asg->seg[i].len == asmg->vtx[i].len);

    // check if circular path existed
    int path_exists;
    uint32_t step;
    uint64_t dist;
    path_exists = asmg_path_exists(asmg, anchor_sid<<1, anchor_sid<<1, 0, COMMON_MAX_MINICIRCLE_SIZE, &step, &dist);
    if (VERBOSE > 0)
        fprintf(stderr, "[M::%s] circlar path %s found between anchor sequence in the original assembly graph: r=%u, d=%lu\n",
                __func__, path_exists? "WAS" : "NOT", step, dist);
    if (!path_exists) {
        // path_exists = extend_and_assemble_from_anchor_sequence();
        if (VERBOSE > 0)
            fprintf(stderr, "[M::%s] circlar path %s found between anchor sequence in the reassembled graph\n",
                    __func__, path_exists? "WAS" : "NOT");
    }

    if (!path_exists) return 1;

    // align reads and find minicircles
    asmg_clean_consensus(scg_meta->scg->utg_asmg);
    scg_read_alignment(scg_meta->sr, scg_meta->ra, scg_meta->scg, n_threads, 0);
    extract_minicircles_with_anchor(scg_meta->ra, anchor_sid, n_threads);

    return 0;
}

int pathfinder_minicircle(char *asg_file, char *mini_annot, scg_meta_t *scg_meta, int n_core, int min_len,
        int min_ex_g, int max_d_len, int max_copy, double max_eval, double min_score, double min_cf, double seq_cf,
        int no_trn, int do_graph_clean, int bubble_size, int tip_size, double weak_cross,
        int out_s, char *out_pref, int n_threads, int VERBOSE)
{
    asg_t *asg;
    hmm_annot_v *annot_v;
    og_component_v *og_components;
    double *seg_annot_score;
    int ret = 0;

    asg = 0;
    annot_v = 0;
    og_components = 0;
    seg_annot_score = 0;

    annot_v = hmm_annot_read(mini_annot, annot_v, OG_MINI);
    if (annot_v == 0) {
        fprintf(stderr, "[E::%s] failed to read the annotation file\n", __func__);
        ret = 1;
        goto do_clean;
    }
    if (annot_v) hmm_annot_index(annot_v);
    if (VERBOSE > 2) hmm_annot_print(annot_v->a, annot_v->n, stderr);

    asg = asg_read(asg_file);
    if (asg == 0) {
        fprintf(stderr, "[E::%s] failed to read the graph: %s\n", __func__, asg_file);
        ret = 1;
        goto do_clean;
    }

    og_components = annot_seq_og_type(annot_v, asg, no_trn, max_eval, n_core, min_len, min_score, &seg_annot_score, VERBOSE);
    if (VERBOSE > 1) print_og_classification_summary(asg, annot_v, og_components, stderr);

    parse_organelle_minicircle(asg, annot_v, og_components, seg_annot_score, scg_meta, n_threads, VERBOSE);

    parse_organelle_component(asg, annot_v, og_components, max_copy, min_ex_g, max_d_len, seq_cf,
            do_graph_clean, bubble_size, tip_size, weak_cross, out_pref, out_s, OG_MINI, VERBOSE);


do_clean:
    asg_destroy(asg);
    hmm_annot_v_destroy(annot_v);
    og_component_v_destroy(og_components);
    free(seg_annot_score);

    return ret;
}

int pathfinder(char *asg_file, char *mito_annot, char *pltd_annot, int n_core, int min_len, int min_ex_g,
        int max_d_len, int max_copy, double max_eval, double min_score, double min_cf, double seq_cf,
        int no_trn, int do_graph_clean, int bubble_size, int tip_size, double weak_cross,
        int out_s, char *out_pref, int VERBOSE)
{
    asg_t *asg;
    hmm_annot_v *annot_v;
    og_component_v *og_components;
    int ret = 0;

    asg = 0;
    annot_v = 0;
    og_components = 0;

    if (mito_annot) annot_v = hmm_annot_read(mito_annot, annot_v, OG_MITO);
    if (pltd_annot) annot_v = hmm_annot_read(pltd_annot, annot_v, OG_PLTD);
    if (annot_v) hmm_annot_index(annot_v);
    if (VERBOSE > 2) hmm_annot_print(annot_v->a, annot_v->n, stderr);

    asg = asg_read(asg_file);
    if (asg == 0) {
        fprintf(stderr, "[E::%s] failed to read the graph: %s\n", __func__, asg_file);
        ret = 1;
        goto do_clean;
    }

    if (VERBOSE > 1) {
        if (VERBOSE > 2) {
            fprintf(stderr, "[M::%s] graph loaded\n", __func__);
            asg_print(asg, stderr, 1);
        }
        fprintf(stderr, "[M::%s] graph stats\n", __func__);
        asg_stat(asg, stderr);
    }

    if (do_graph_clean && min_cf > .0) clean_graph_by_sequence_coverage(asg, min_cf, max_copy, VERBOSE);

    og_components = annot_seq_og_type(annot_v, asg, no_trn, max_eval, n_core, min_len, min_score, 0, VERBOSE);
    if (VERBOSE > 1) print_og_classification_summary(asg, annot_v, og_components, stderr);

    // graph will be changed with extra copies of sequences added
    if (mito_annot)
        parse_organelle_component(asg, annot_v, og_components, max_copy, min_ex_g, max_d_len, seq_cf,
                do_graph_clean, bubble_size, tip_size, weak_cross, out_pref, out_s, OG_MITO, VERBOSE);
    if (pltd_annot)
        parse_organelle_component(asg, annot_v, og_components, max_copy, min_ex_g, max_d_len, seq_cf,
                do_graph_clean, bubble_size, tip_size, weak_cross, out_pref, out_s, OG_PLTD, VERBOSE);

do_clean:
    asg_destroy(asg);
    hmm_annot_v_destroy(annot_v);
    og_component_v_destroy(og_components);

    return ret;
}

#ifdef PATHFINDER_MAIN
#include "ketopt.h"

int VERBOSE = 0;
double realtime0;

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
    { "max-bubble",     ko_required_argument, 312 },
    { "max-tip",        ko_required_argument, 313 },
    { "weak-cross",     ko_required_argument, 314 },
    { "mito-annot",     ko_required_argument, 'm' },
    { "pltd-annot",     ko_required_argument, 'p' },
    { "core-gene",      ko_required_argument, 'g' },
    { "min-score",      ko_required_argument, 's' },
    { "min-gain",       ko_required_argument, 'a' },
    { "min-s-cov",      ko_required_argument, 'q' },
    { "max-copy",       ko_required_argument, 'c' },
    { "version",        ko_no_argument,       'V' },
    { "help",           ko_no_argument,       'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    const char *opt_str = "a:m:p:f:g:q:s:c:o:Vv:h";
    ketopt_t opt = KETOPT_INIT;
    int c, out_s, out_c, max_copy, ret = 0;
    FILE *fp_help;
    char *out_pref, *mito_annot, *pltd_annot, *ec_tag, *kc_tag, *sc_tag;
    int no_trn, n_core, min_len, min_ex_g, max_d_len, bubble_size, tip_size, do_graph_clean;
    double max_eval, min_score, min_cf, seq_cf, weak_cross;

    sys_init();

    out_s = -1;
    out_c = 0;
    out_pref = "oatk.asm";
    max_copy = 10;
    fp_help = stderr;
    mito_annot = 0;
    pltd_annot = 0;
    ec_tag = 0;
    kc_tag = 0;
    sc_tag = 0;
    no_trn = 1;
    do_graph_clean = 1;
    bubble_size = 100000;
    tip_size = 10000;
    weak_cross = 0.3;
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
        else if (c == 312) bubble_size = atoi(opt.arg);
        else if (c == 313) tip_size = atoi(opt.arg);
        else if (c == 314) weak_cross = atof(opt.arg);
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
        fprintf(fp_help, "Usage: pathfinder [options] <file>[.gfa[.gz]] [path_str]\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "  Input/Output:\n");
        fprintf(fp_help, "    -m FILE              mitochondria core gene annotation file [NULL]\n");
        fprintf(fp_help, "    -p FILE              plastid core gene annotation file [NULL]\n");
        fprintf(fp_help, "    -o FILE              prefix of output files [%s]\n", out_pref);
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
        fprintf(fp_help, "    --max-bubble INT     maximum bubble size for assembly graph clean [%d]\n", bubble_size);
        fprintf(fp_help, "    --max-tip    INT     maximum tip size for assembly graph clean [%d]\n", tip_size);
        fprintf(fp_help, "    --weak-cross FLOAT   maximum relative edge coverage for weak crosslink clean [%.2f]\n", weak_cross);
        fprintf(fp_help, "Example: ./pathfinder -m asm_annot.mito.txt -p asm_annot.pltd.txt -o oatk.asm asm.gfa\n");
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
    
    ret = pathfinder(argv[opt.ind], mito_annot, pltd_annot, n_core, min_len, min_ex_g, max_d_len, max_copy, 
            max_eval, min_score, min_cf, seq_cf, no_trn, do_graph_clean, bubble_size, tip_size, weak_cross,
            out_s, out_pref, VERBOSE);
    
    if (ret) {
        fprintf(stderr, "[E::%s] failed to analysis the GFA file\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (fflush(stdout) == EOF) {
        fprintf(stderr, "[E::%s] failed to write the results\n", __func__);
        exit(EXIT_FAILURE);
    }

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
#endif
