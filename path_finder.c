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
#include "version.h"

#undef DEBUG_MINICIRCLE_REPEAT_UNIT

KHASHL_MAP_INIT(KH_LOCAL, kh_s32_t, kh_s32, khint32_t, uint32_t, kh_hash_dummy, kh_eq_generic);

typedef struct { size_t n, m; uint32_t *a; } v_u32_t;
typedef struct { size_t n, m; v_u32_t *a; } lv_u32_t;

static void lv_u32_destroy(lv_u32_t *lv)
{
    size_t i;
    for (i = 0; i < lv->n; ++i)
        kv_destroy(lv->a[i]);
    kv_destroy(*lv);
}

static lv_u32_t parse_subgraph(asg_t *asg)
{
    uint32_t i, j, n_seg, *vlist, nv;
    lv_u32_t lv_list;
    asmg_t *g;
    uint8_t *visited;

    g = asg->asmg;
    n_seg = asg->n_seg;
    kv_init(lv_list);
    MYCALLOC(visited, n_seg);
    for (i = 0; i < n_seg; ++i) {
        if (visited[i] || g->vtx[i].del) continue;
        vlist = asmg_subgraph(g, &i, 1, 0, 0, &nv, 0);
        assert(nv > 0);
        v_u32_t lv = {nv, nv, vlist};
        kv_push(v_u32_t, lv_list, lv);
        for (j = 0; j < nv; ++j)
            visited[vlist[j]] = 1;
    }
    free(visited);

    return lv_list;
}

FILE *fopen1(char *fn)
{
    FILE *fo;
    fo = fopen(fn, "w");
    if (!fo) {
        fprintf(stderr, "[E::%s] failed to open file %s to write\n", __func__, fn);
        exit(EXIT_FAILURE);
    }
    return fo;
}

static void parse_organelle_component(asg_t *asg, hmm_annot_db_t *annot_db, og_component_v *og_components, int min_s_len,
        int max_copy, int min_ext_g, double seq_cf, int do_clean, double min_cf, double min_score, double max_eval,
        int bubble_size, int tip_size, double weak_cross, char *out_pref, int out_opt, OG_TYPE_t og_type, int VERBOSE)
{
    assert(og_type == OG_MITO || og_type == OG_PLTD || og_type == OG_MINI);

    char *fn_str;
    MYMALLOC(fn_str, strlen(out_pref) + 64);
    //sprintf(fn_str, "%s.%s.stats", out_pref, OG_TYPES[og_type]);
    //FILE *out_stats = fopen1(fn_str);
    sprintf(fn_str, "%s.%s.ctg.fasta", out_pref, OG_TYPES[og_type]);
    FILE *out_ctg = fopen1(fn_str);
    sprintf(fn_str, "%s.%s.ctg.bed", out_pref, OG_TYPES[og_type]);
    FILE *out_ctg_bed = fopen1(fn_str);
    sprintf(fn_str, "%s.%s.gfa", out_pref, OG_TYPES[og_type]);
    FILE *out_gfa = fopen1(fn_str);
    sprintf(fn_str, "%s.%s.bed", out_pref, OG_TYPES[og_type]);
    FILE *out_gfa_bed = fopen1(fn_str);
    free(fn_str);

    uint32_t i, j, v, opt_circ, b_length;
    int c, ext_g, all_g, absent;
    uint64_t n_seg;
    double opt_coverage, g_diff, c_diff, h_score, b_score, score, score1;
    og_component_t *component;
    kh_s32_t *h_genes, *b_genes;
    khint32_t k;
    kvec_t(uint32_t) sub_v;
    asmg_t *o_asmg;
    hmm_annot_bed6_db_t *bed_annots;

    o_asmg = asg->asmg; // original copy of asmg
    n_seg = asg->n_seg;
    h_genes = kh_s32_init();
    b_genes = kh_s32_init();
    kv_init(sub_v);
    c = 0;
    opt_circ = 0;
    opt_coverage = 0.;
    g_diff = 0.85; // score difference for new genes
    c_diff = 0.6; // gene density difference
    bed_annots = hmm_annot_bed6_db_init();

    // build gene table
    // best score for each gene
    for (i = 0; i < og_components->n; ++i) {
        component = &og_components->a[i];
        if (component->type != og_type)
            continue;
        for (j = 0; j < component->ng; ++j) {
            if ((component->g[j]>>32 & 0x3) != (uint32_t) og_type)
                continue;
            k = kh_s32_put(h_genes, component->g[j]>>32, &absent);
            if (absent || kh_val(h_genes, k) < (uint32_t) component->g[j])
                kh_val(h_genes, k) = (uint32_t) component->g[j];
        }
    }
    h_score = 0; // total score
    for (k = (khint32_t) 0; k < kh_end(h_genes); ++k)
        if (kh_exist(h_genes, k))
            h_score += kh_val(h_genes, k);
    if (VERBOSE > 0)
        fprintf(stderr, "[M::%s] total gene score for the organelle: type, %s; score, %.1f\n",
                __func__, OG_TYPES[og_type], h_score);

    b_score = 0;
    b_length = 0;
    for (i = 0; i < og_components->n; ++i) {
        component = &og_components->a[i];
        if (component->type != og_type)
            continue;
        if (VERBOSE > 0)
            fprintf(stderr, "[M::%s] processing subgraph seeding from %s: "
                    "type, %s; score, %.1f; sscore, %.1f; len, %u; nv, %u; ng, %u\n", 
                    __func__, asg->seg[component->v[0]].name, OG_TYPES[component->type], 
                    component->score, component->sscore, component->len, component->nv, component->ng);
        
        // find number of genes included in the subgraph
        ext_g = all_g = 0;
        for (j = 0; j < component->ng; ++j) {
            if ((component->g[j]>>32 & 0x3) != (uint32_t) og_type)
                continue;
            k = kh_s32_get(b_genes, component->g[j]>>32);
            score = k == kh_end(b_genes)? 0 : kh_val(b_genes, k);
            score1 = (uint32_t) component->g[j];
            if (score1 >= min_score && score1 >= score)
                ++ext_g; // extra gene count
            if (score1 >= score * g_diff)
                ++all_g; // all gene count
        }
        
        // skip component with no enough extra genes
        if (ext_g < min_ext_g && all_g < kh_size(b_genes) * c_diff) {
            if (VERBOSE > 0)
                fprintf(stderr, "[M::%s] subgraph seeding from %s SKIPPED due to insufficient gene gain (%d)\n",
                        __func__, asg->seg[component->v[0]].name, ext_g);
            continue;
        }

        // do a specical check for OG_PLTD
        // require enough gene density
        // TODO is this a good strategy
        if (og_type == OG_PLTD &&
                b_length + component->len > COMMON_AVG_PLTD_SIZE && // the size exceeds the expected PLTD size
                component->score * b_length < b_score * component->len * c_diff) {
            if (VERBOSE > 0)
                fprintf(stderr, "[M::%s] subgraph seeding from %s SKIPPED due to low PLTD gene density (%.1f/%u < %.1f/%u*%.1f)\n",
                        __func__, asg->seg[component->v[0]].name, component->score, component->len, b_score, b_length, c_diff);
            continue;
        }

        // add subgraph genes to the organelle gene table
        for (j = 0; j < component->ng; ++j) {
            if ((component->g[j]>>32 & 0x3) != (uint32_t) og_type)
                continue;
            k = kh_s32_put(b_genes, component->g[j]>>32, &absent);
            if (absent || kh_val(b_genes, k) < (uint32_t) component->g[j])
                kh_val(b_genes, k) = (uint32_t) component->g[j];
        }
        b_score += component->score;
        b_length += component->len;

        // processing the subgraph
        // make a copy of the orginal asmg
        // asg->asmg = asg_make_asmg_copy(o_asmg, 0);
        asg->asmg = component->asmg;
        // make subgraph
        // asmg_subgraph(asg->asmg, component->v, 1, 0, 0, 0, 1);

        if (VERBOSE > 1) {
            if (VERBOSE > 2) {
                fprintf(stderr, "[M::%s] subgraph created\n", __func__);
                asg_print(asg, stderr, 1);
            }    
            fprintf(stderr, "[M::%s] subgraph stats\n", __func__);
            asg_stat(asg, stderr);
        }

        if (do_clean) {
            // do basic cleanup
            uint64_t cleaned = 1;
            while (cleaned) {
                cleaned = 0;
                cleaned += asmg_pop_bubble(asg->asmg, bubble_size, 0, 0, 1, 0, VERBOSE);
                cleaned += asmg_remove_weak_crosslink(asg->asmg, weak_cross, 10, 0, VERBOSE);
                cleaned += asmg_drop_tip(asg->asmg, INT32_MAX, tip_size, 1, 0, VERBOSE);
            }
            
            if (VERBOSE > 1) {
                if (VERBOSE > 2) {
                    fprintf(stderr, "[M::%s] subgraph after basic clean\n", __func__);
                    asg_print(asg, stderr, 1);
                }
                fprintf(stderr, "[M::%s] subgraph stats after basic clean\n", __func__);
                asg_stat(asg, stderr);
            }
            // TODO filter low coverage arcs
        }
        
        if (asmg_vtx_n1(asg->asmg) == 0) {
            asg->asmg = o_asmg; // roll back to the original copy
            continue;
        }
        
        uint32_t clen = asg_seg_len(asg);

        int *copy_number = 0;
        double avg_coverage;
        avg_coverage = graph_sequence_coverage_precise(asg, min_cf, 1, max_copy, &copy_number);
        if (VERBOSE > 0)
            fprintf(stderr, "[M::%s] estimated per-copy sequence coverage: %.3f\n", __func__, avg_coverage);

        if (og_type == OG_MITO && opt_coverage > 0 &&
                (avg_coverage < opt_coverage * min_cf ||
                 avg_coverage * min_cf > opt_coverage)) {
            // the sequence coverage of this subgraph is too low or too high
            // clean and roll back
            free(copy_number);
            asg->asmg = o_asmg; // roll back to the original copy
            continue;
        }
        if (opt_coverage == 0.) opt_coverage = avg_coverage;

        // make a copy of the graph
        asg_t *asg_copy = asg_make_copy(asg);
        kh_u32_t *seg_dups = sequence_duplication_by_copy_number(asg_copy, copy_number, 0);
        if (VERBOSE > 1) {
            if (VERBOSE > 2) {
                fprintf(stderr, "[M::%s] subgraph after making seg copies\n", __func__);
                asg_print(asg_copy, stderr, 1);
            }
            fprintf(stderr, "[M::%s] seg copies included subgraph stats\n", __func__);
            asg_stat(asg_copy, stderr);
        }
        path_v paths;
        kv_init(paths);
        graph_path_finder(asg_copy, seg_dups, &paths, seq_cf, og_type == OG_PLTD);
        kh_u32_destroy(seg_dups);
        asg_destroy(asg_copy);

        if (paths.n == 0) {
            // not possible to solve the graph
            // likely because the graph is too complicated
            if (VERBOSE > 0)
                fprintf(stderr, "[M::%s] subgraph seeding from %s is unresolvable, output unitigs as unassembled\n",
                        __func__, asg->seg[component->v[0]].name);
            // FIXME sequences removed by graph clean are not included
            // asg_print_fa(asg, out_ctg, 60);            
            for (j = 0; j < component->nv; ++j) {
                v = component->v[j];
                if (asg->asmg->vtx[v].del) continue;
                // write sequence
                v <<= 1;
                path_t path = {0, 1, 0, 1, &v, asg->seg[v>>1].len,
                    (double) asg->seg[v>>1].len * asg->seg[v>>1].cov, 0};
                print_seq(asg, &path, out_ctg, ++c, 0, 60, 100);
                path_add_hmm_annot_bed6(bed_annots, annot_db, asg, &path, c, 0, 100, og_type, max_eval);
            }
            kv_push(uint32_t, sub_v, i);
        } else {
            uint32_t b, is_circ;
            kvec_t(uint32_t) v_pb;
            double f;
            if (og_type == OG_PLTD) // only for pltd
                for (j = 0; j < paths.n; ++j)
                    path_rotate(asg, &paths.a[j], annot_db, 2);
            path_sort(&paths);
            // TODO improve the selection criteria
            kv_init(v_pb);
            b = select_best_seq(asg, &paths, 0, out_opt, seq_cf, 0, og_type == OG_PLTD);
            f = sequence_covered_by_path(asg, &paths.a[b], clen);
            is_circ = paths.a[b].circ;
            kv_push(uint32_t, v_pb, b);
            if (VERBOSE > 0)
                fprintf(stderr, "[M::%s] best path after first pass: type, %s; coverage, %.3f\n", 
                        __func__, is_circ? "circular" : "linear", f);
            if (!is_circ || f < 1.0) {
                // adjust the sequence copy number and redo path finding
                asg_copy = asg_make_copy(asg);
                double adjusted_avg_coverage;
                int updated = adjust_sequence_copy_number_by_graph_layout(asg_copy, avg_coverage, &adjusted_avg_coverage, copy_number, max_copy, 10);
                if (updated) {
                    if (VERBOSE > 0)
                        fprintf(stderr, "[M::%s] adjusted per-copy sequence coverage: %.3f\n", __func__, adjusted_avg_coverage);
                    // make another copy of the graph
                    asg_t *asg_copy1 = asg_make_copy(asg_copy);
                    // do path finding
                    // allow deletion of segs here if the adjusted copy number is zero
                    kh_u32_t *seg_dups1 = sequence_duplication_by_copy_number(asg_copy1, copy_number, 1);

                    // subgraph may be divided into multiple subgraphs
                    path_v paths1;
                    uint32_t b1, is_circ1;
                    kvec_t(uint32_t) v_pb1;
                    double f1;
                    asmg_t *o_g1, *g1;

                    // vertex list of subgraphs
                    lv_u32_t vlist = parse_subgraph(asg_copy1);
                    is_circ1 = 1;
                    f1 = .0;
                    o_g1 = asg_copy1->asmg;
                    kv_init(paths1);
                    kv_init(v_pb1);
                    // do path finding for each subgraph
                    for (j = 0; j < vlist.n; ++j) {
                        // initialise the subgraph
                        g1 = asg_make_asmg_copy(o_g1, 0);
                        for (v = 0; v < asg_copy1->n_seg; ++v)
                            g1->vtx[v].del = 1;
                        for (v = 0; v < vlist.a[j].n; ++v)
                            g1->vtx[vlist.a[j].a[v]].del = 0;
                        // remove dirty arcs
                        for (v = 0; v < g1->n_arc; ++v)
                            if (g1->vtx[g1->arc[v].v>>1].del || g1->vtx[g1->arc[v].w>>1].del)
                                g1->arc[v].del = 1;

                        path_v tmp_paths1 = {0, 0, 0};
                        asg_copy1->asmg = g1;
                        if (VERBOSE > 1) {
                            if (VERBOSE > 2) {
                                fprintf(stderr, "[M::%s] subgraph after making seg copies\n", __func__);
                                asg_print(asg_copy1, stderr, 1);
                            }
                            fprintf(stderr, "[M::%s] seg copies included subgraph stats\n", __func__);
                            asg_stat(asg_copy1, stderr);
                        }
                        graph_path_finder(asg_copy1, seg_dups1, &tmp_paths1, seq_cf, og_type == OG_PLTD);

                        if (og_type == OG_PLTD) // only for pltd
                            for (j = 0; j < tmp_paths1.n; ++j)
                                path_rotate(asg_copy1, &tmp_paths1.a[j], annot_db, 2);
                        path_sort(&tmp_paths1);
                        // TODO improve the selection criteria
                        b1 = select_best_seq(asg_copy1, &tmp_paths1, 0, out_opt, seq_cf, 0, og_type == OG_PLTD);    
                        f1 += sequence_covered_by_path(asg_copy1, &tmp_paths1.a[b1], clen);
                        is_circ1 &= tmp_paths1.a[b1].circ;
                        // add all paths
                        kv_push(uint32_t, v_pb1, b1 + paths1.n);
                        kv_pushn(path_t, paths1, tmp_paths1.a, tmp_paths1.n);
                        // do not destroy path_v
                        // path_v_destroy(&tmp_paths1);
                        // only need to free tmp_paths1.a
                        kv_destroy(tmp_paths1);
                        asmg_destroy(g1);
                    }
                    asg_copy1->asmg = o_g1;

                    kh_u32_destroy(seg_dups1);
                    asg_destroy(asg_copy1);

                    if (VERBOSE > 0) fprintf(stderr, "[M::%s] best path in second pass: type, %s; coverage, %.3f\n",
                            __func__, is_circ1? "circular" : "linear", f1);

                    // compare new path to old to make final choice
                    // TODO improve comparison criteria
                    // if ((is_circ1 && f1 >= seq_cf) || ((is_circ1 || !is_circ) && f1 > f)) {
                    if ((is_circ1 == is_circ && f1 > f) || // second pass covered more sequences 
                            (is_circ1 > is_circ && f1 >= f * seq_cf) || // second pass is circle and covered at least seq_cf of first pass
                            (is_circ1 < is_circ && f1 * seq_cf >= f)) { // second pass is not circle but covered significantly more sequence
                        // choose the new one
                        f = f1;
                        is_circ = is_circ1;
                        kv_copy(uint32_t, v_pb, v_pb1);
                        path_v_destroy(&paths);
                        paths = paths1;
                    } else {
                        path_v_destroy(&paths1);
                    }
                    
                    kv_destroy(v_pb1);
                    lv_u32_destroy(&vlist);

                    if (VERBOSE > 0)
                        fprintf(stderr, "[M::%s] best path after second pass: type, %s; coverage, %.3f\n",
                                __func__, is_circ? "circular" : "linear", f);

                }
                asg_destroy(asg_copy);
            }

            if (is_circ || !opt_circ || clen >= min_s_len) {
                if (!opt_circ) opt_circ = is_circ;
                
                kv_push(uint32_t, sub_v, i);

                // path_stats(asg, &paths, out_stats);
                int *incl;
                MYCALLOC(incl, n_seg);
                for (j = 0; j < component->nv; ++j) {
                    v = component->v[j];
                    if (!asg->asmg->vtx[v].del)
                        incl[v] = 1;
                }
                for (j = 0; j < v_pb.n; ++j) {
                    path_t *path = &paths.a[v_pb.a[j]];
                    print_seq(asg, path, out_ctg, ++c, 0, 60, 100);
                    path_add_hmm_annot_bed6(bed_annots, annot_db, asg, path, c, 0, 100, og_type, max_eval);
                    for (v = 0; v < path->nv; ++v)
                        incl[path->v[v]>>1] = 0;
                }
                
                // revisit to write sequences not included in the paths
                for (j = 0; j < component->nv; ++j) {
                    v = component->v[j];
                    if (!incl[v] || asg->seg[v].len < min_s_len) continue;
                    // write sequence
                    v <<= 1;
                    path_t path = {0, 1, 0, 1, &v, asg->seg[v>>1].len, 
                        (double) asg->seg[v>>1].len * asg->seg[v>>1].cov, 0};
                    print_seq(asg, &path, out_ctg, ++c, 0, 60, 100);
                    path_add_hmm_annot_bed6(bed_annots, annot_db, asg, &path, c, 0, 100, og_type, max_eval);
                }
                free(incl);

                if (VERBOSE > 0)
                    fprintf(stderr, "[M::%s] processing subgraph seeding from %s DONE, %d better genes gained, total score %.1f\n",
                            __func__, asg->seg[component->v[0]].name, ext_g, b_score);
            }

            kv_destroy(v_pb);
        }

        path_v_destroy(&paths);

        free(copy_number);

        asg->asmg = o_asmg; // roll back to the original copy
    }
    
    hmm_annot_print_bed6(bed_annots, out_ctg_bed, 1);

    // extract and output organelle gfa components
    if (sub_v.n > 0) {
        // in subgraph all vtx and arcs deleted during expansion are back
        uint32_t cov, del;
        asmg_t *g, *g1;
        g = asg_make_asmg_copy(og_components->a[sub_v.a[0]].asmg, 0);
        for (i = 1; i < sub_v.n; ++i) {
            g1 = og_components->a[sub_v.a[i]].asmg;
            // merge g1 to g
            for (j = 0; j < g->n_vtx; ++j) {
                cov = 0;
                del = 1;
                if (!g->vtx[j].del)  del = 0, cov += g->vtx[j].cov;
                if (!g1->vtx[j].del) del = 0, cov += g1->vtx[j].cov;
                if (del) continue;
                if (cov > o_asmg->vtx[j].cov) cov = o_asmg->vtx[j].cov;
                g->vtx[j].del = del;
                g->vtx[j].cov = cov;
            }
            for (j = 0; j < g->n_arc; ++j) {
                cov = 0;
                del = 1;
                if (!g->arc[j].del)  del = 0, cov += g->arc[j].cov;
                if (!g1->arc[j].del) del = 0, cov += g1->arc[j].cov;
                if (del) continue;
                if (cov > o_asmg->arc[j].cov) cov = o_asmg->arc[j].cov;
                g->arc[j].del = del;
                g->arc[j].cov = cov;
            }
        }
        asg->asmg = g;
        uint64_t nv = 0;
        char **vlist = asg_vtx_name_list(asg, &nv);
        hmm_annot_formatted_print_sname_list(annot_db, vlist, nv, out_gfa_bed, og_type, max_eval, 1);
        free(vlist);
        asg_print(asg, out_gfa, 0);
        asmg_destroy(asg->asmg);
        asg->asmg = o_asmg;
    }
    kv_destroy(sub_v);
    kh_s32_destroy(h_genes);
    kh_s32_destroy(b_genes);
    hmm_annot_bed6_db_destroy(bed_annots);

    //fclose(out_stats);
    fclose(out_ctg);
    fclose(out_ctg_bed);
    fclose(out_gfa);
    fclose(out_gfa_bed);
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

static inline void rev_array(uint32_t *arr, uint32_t n)
{
    uint32_t tmp, i, j;
    if (n == 0) return;
    for (i = 0, j = n - 1; i < j; ++i, --j) {
        tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
}

static int path_cmpfunc(const void *a, const void *b)
{
    path_t *x, *y;
    x = (path_t *) a;
    y = (path_t *) b;
    if (x->nv != y->nv)
        return (x->nv > y->nv) - (x->nv < y->nv);

    uint32_t i, nv, *vx, *vy;
    nv = x->nv;
    vx = x->v;
    vy = y->v;

    for (i = 0; i < nv; ++i)
        if (vx[i] != vy[i])
            return (vx[i] > vy[i]) - (vx[i] < vy[i]);

    return 0;
}

static int extract_minicircles_with_anchor(scg_ra_v *ra, scg_t *scg, uint64_t anchor_sid, int n_threads, path_v *paths)
{
    uint64_t i, j, k, nr, *mcs;
    
    nr = ra->n;
    MYMALLOC(mcs, nr);
    mc_shared_t mc_shared = {ra, anchor_sid, mcs};
    kt_for(n_threads, minicircle_analysis_thread, &mc_shared, nr);

    for (i = 0; i < nr; ++i) {
        if (mcs[i] == UINT64_MAX) continue;
        scg_ra_t *r = &ra->a[i];
#ifdef DEBUG_MINICIRCLE_REPEAT_UNIT
        scg_ra_print(r, stderr);
        fprintf(stderr, "[DEBUG_MINICIRCLE_REPEAT_UNIT::%s] repeat unit: beg %lu; end %u; %c\n", __func__, 
                mcs[i]>>33, (uint32_t) (mcs[i]>>1), "+-"[mcs[i]&1]);
#endif
        // collect paths
        kvec_t(uint32_t) vt;
        kv_init(vt);
        uint32_t n;
        for (j = (mcs[i]>>33), n = (uint32_t) (mcs[i]>>1); j <= n; ++j)
            kv_push(uint32_t, vt, (uint32_t) (r->a[j].uid));
        if (mcs[i]&1) {
            // do reverse complement
            rev_array(vt.a+1, vt.n-1);
            for (j = 0; j < vt.n; ++j)
                vt.a[j] ^= 1;
        }
    
        path_t path = {0, vt.n, 1, 0, vt.a, 0, 0., 0.};
        kv_push(path_t, *paths, path);
    }
    free(mcs);

    if (paths->n == 0) return 0;

    // sort paths and remove duplicates
    qsort(paths->a, paths->n, sizeof(path_t), path_cmpfunc);
    j = 0;
    k = 1;
    for (i = 1; i < paths->n; ++i) {
        int c = path_cmpfunc(&paths->a[i], &paths->a[j]);
        if (c == 0) {
            free(paths->a[i].v);
            --k;
        } else {
            if (k != i)
                paths->a[k] = paths->a[i];
            j = k;
        }
        ++k;
    }
    paths->n = k;

    // now update path info
    uint32_t nv, len, len1, cov, *vt;
    double wlen;
    asmg_t *g;
    asmg_arc_t *a;
    
    g = scg->utg_asmg;
    for (i = 0; i < paths->n; ++i) {
        path_t *path = &paths->a[i];
        nv = path->nv;
        vt = path->v;
        
        a = asmg_arc1(g, vt[nv-1], vt[0]);
        assert(!!a); // should always be a circle
        len = g->vtx[vt[0]>>1].len;
        cov = g->vtx[vt[0]>>1].cov;
        wlen = (double) cov * len;
        len -= a->ls, wlen -= cov * a->ls; // circular path
        
        for (j = 1; j < nv; ++j) {
            len1 = g->vtx[vt[j]>>1].len;
            cov = g->vtx[vt[j]>>1].cov;
            len += len1;
            wlen += (double) cov * len1;
            a = asmg_arc1(g, vt[j-1], vt[j]);
            assert(!!a);
            len -= a->ls;
            wlen -= (double) cov * a->ls;
        }

        path->len = len;
        path->wlen = wlen;
    }

    return paths->n;
}

static int parse_organelle_minicircle(asg_t *asg, hmm_annot_db_t *annot_db, og_component_v *og_components,
        double *seg_annot_score, scg_meta_t *scg_meta, int n_threads, char *out_pref, int out_opt, 
        double max_eval, double seq_cf, int VERBOSE)
{
    if (og_components->n == 0) {
        if (VERBOSE > 0)
            fprintf(stderr, "[M::%s] no OG component found\n", __func__);
        return 1;
    }

    char *fn_str;
    MYMALLOC(fn_str, strlen(out_pref) + 64);
    //sprintf(fn_str, "%s.%s.stats", out_pref, OG_TYPES[OG_MINI]);
    //FILE *out_stats = fopen1(fn_str);
    sprintf(fn_str, "%s.%s.ctg.fasta", out_pref, OG_TYPES[OG_MINI]);
    FILE *out_ctg = fopen1(fn_str);
    sprintf(fn_str, "%s.%s.ctg.bed", out_pref, OG_TYPES[OG_MINI]);
    FILE *out_ctg_bed = fopen1(fn_str);
    sprintf(fn_str, "%s.%s.gfa", out_pref, OG_TYPES[OG_MINI]);
    FILE *out_gfa = fopen1(fn_str);
    sprintf(fn_str, "%s.%s.bed", out_pref, OG_TYPES[OG_MINI]);
    FILE *out_gfa_bed = fopen1(fn_str);
    free(fn_str);

    uint32_t i, anchor_sid;
    uint64_t n_seg;
    double max_s;
    og_component_t *component;
    asmg_t *asmg;
    hmm_annot_bed6_db_t *bed_annots;

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
        // WARN assembly graph and anchor sequence will be updated here
        // remember to update asg
        // path_exists = extend_and_assemble_from_anchor_sequence(scg_meta, asg, &anchor_sid);
        if (VERBOSE > 0)
            fprintf(stderr, "[M::%s] circlar path %s found between anchor sequence in the reassembled graph\n",
                    __func__, path_exists? "WAS" : "NOT");
    }

    path_v paths = {0, 0, 0};

    if (path_exists) {
        // align reads and find minicircles
        asmg_clean_consensus(scg_meta->scg->utg_asmg);
        scg_read_alignment(scg_meta->sr_db, scg_meta->ra_db, scg_meta->scg, n_threads, 0);

#ifdef DEBUG_MINICIRCLE_REPEAT_UNIT
        scg_rv_print(scg_meta->ra_db, stderr);
#endif

        scg_consensus(scg_meta->sr_db, scg_meta->scg, 0, 0, 0);
        extract_minicircles_with_anchor(scg_meta->ra_db, scg_meta->scg, anchor_sid, n_threads, &paths);
    }

    // prepare assembly subgraph for output
    asmg_t *o_asmg = asg->asmg;
    asg->asmg = asg_make_asmg_copy(o_asmg, 0);
    asmg_subgraph(asg->asmg, &anchor_sid, 1, 0, 0, 0, 1);
    bed_annots = hmm_annot_bed6_db_init();

    if (paths.n == 0) {
        // not possible to solve the graph
        // likely because the graph is not circularisable or too complicated
        if (VERBOSE > 0)
            fprintf(stderr, "[M::%s] subgraph seeding from %s is unresolvable, output unitigs as unassembled\n",
                    __func__, asg->seg[anchor_sid].name);
        // FIXME sequences removed by graph clean are not included
        asg_print_fa(asg, stdout, 60);
        uint32_t c, v;
        c = 0;
        for (i = 0; i < component->nv; ++i) {
            v = component->v[i];
            if (asg->asmg->vtx[v].del) continue;
            // write sequence
            v <<= 1;
            path_t path = {0, 1, 0, 1, &v, asg->seg[v>>1].len,
                (double) asg->seg[v>>1].len * asg->seg[v>>1].cov, 0};
            print_seq(asg, &path, out_ctg, ++c, 0, 60, 100);
            path_add_hmm_annot_bed6(bed_annots, annot_db, asg, &path, c, 0, 100, OG_MINI, max_eval);
        }
    } else {
        uint32_t b;
        path_sort(&paths);
        // TODO improve the selection criteria
        b = select_best_seq(asg, &paths, 0, out_opt, seq_cf, 0, 0);
        // path_stats(asg, &paths, out_stats);
        print_seq(asg, &paths.a[b], out_ctg, 1, 0, 60, 100);
        path_add_hmm_annot_bed6(bed_annots, annot_db, asg, &paths.a[b], 1, 0, 100, OG_MINI, max_eval);
    }
    
    hmm_annot_print_bed6(bed_annots, out_ctg_bed, 1);

    uint64_t nv = 0;
    char **vlist = asg_vtx_name_list(asg, &nv);
    hmm_annot_formatted_print_sname_list(annot_db, vlist, nv, out_gfa_bed, OG_MINI, max_eval, 1);
    free(vlist);
    asg_print(asg, out_gfa, 0);
    asmg_destroy(asg->asmg);
    asg->asmg = o_asmg; // roll back to the original copy

    path_v_destroy(&paths);
    hmm_annot_bed6_db_destroy(bed_annots);

    //fclose(out_stats);
    fclose(out_ctg);
    fclose(out_ctg_bed);
    fclose(out_gfa);
    fclose(out_gfa_bed);

    return 0;
}

int pathfinder_minicircle(char *asg_file, char *mini_annot, scg_meta_t *scg_meta, int min_len,
        int min_ext_g, int max_copy, double max_eval, double min_score, double min_cf, double seq_cf,
        int no_trn, int no_rrn, int do_graph_clean, int bubble_size, int tip_size, double weak_cross,
        int out_opt, char *out_pref, int n_threads, int VERBOSE)
{
    asg_t *asg;
    hmm_annot_db_t *annot_db;
    og_component_v *og_components;
    double *seg_annot_score;
    int ret = 0;

    asg = 0;
    annot_db = 0;
    og_components = 0;
    seg_annot_score = 0;

    asg = asg_read(asg_file);
    if (asg == 0) {
        fprintf(stderr, "[E::%s] failed to read the graph: %s\n", __func__, asg_file);
        ret = 1;
        goto do_clean;
    }

    annot_db = hmm_annot_read(mini_annot, annot_db, OG_MINI);
    if (annot_db == 0) {
        fprintf(stderr, "[E::%s] failed to read the annotation file\n", __func__);
        ret = 1;
        goto do_clean;
    }

    if (VERBOSE > 2) hmm_annot_db_print(annot_db, stderr);

    seg_annot_score = get_sequence_annot_score(annot_db, asg, no_trn, no_rrn, max_eval, 0, VERBOSE);
    og_components = annot_subgraph_og_type(annot_db, asg, no_trn, no_rrn, max_eval, 0, min_len, min_score, 1, VERBOSE);
    if (!og_components) {
        fprintf(stderr, "[E::%s] no organelle component found\n", __func__);
        ret = 1;
        goto do_clean;
    }
    if (VERBOSE > 1) print_og_classification_summary(asg, annot_db, og_components, stderr);

    parse_organelle_minicircle(asg, annot_db, og_components, seg_annot_score, scg_meta, n_threads, out_pref, out_opt, max_eval, seq_cf, VERBOSE);

do_clean:
    asg_destroy(asg);
    hmm_annot_db_destroy(annot_db);
    og_component_v_destroy(og_components);
    free(seg_annot_score);

    return ret;
}

int pathfinder(char *asg_file, char *mito_annot, char *pltd_annot, int min_len, int ext_p, int ext_m,
        int max_copy, double max_eval, double min_score, double min_cf, double seq_cf,
        int no_trn, int no_rrn, int do_graph_clean, int bubble_size, int tip_size, double weak_cross,
        int out_opt, char *out_pref, int VERBOSE)
{
    asg_t *asg;
    hmm_annot_db_t *annot_db;
    og_component_v *og_components;
    int ret = 0;

    asg = 0;
    annot_db = 0;
    og_components = 0;

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

    if (mito_annot) annot_db = hmm_annot_read(mito_annot, annot_db, OG_MITO);
    if (pltd_annot) annot_db = hmm_annot_read(pltd_annot, annot_db, OG_PLTD);
    
    if (VERBOSE > 2) hmm_annot_db_print(annot_db, stderr);

    // this is not necessary
    // if (do_graph_clean && min_cf > .0) clean_graph_by_sequence_coverage(asg, min_cf, max_copy, VERBOSE);
    og_components = asg_annotation(annot_db, asg, no_trn, no_rrn, max_eval, 0, min_len, min_score, 1, VERBOSE);
    // og_components = annot_subgraph_og_type(annot_db, asg, no_trn, no_rrn, max_eval, 0, min_len, min_score, 1, VERBOSE);
    if (!og_components) {
        fprintf(stderr, "[E::%s] no organelle component found\n", __func__);
        ret = 1;
        goto do_clean;
    }

    if (VERBOSE > 1) print_og_classification_summary(asg, annot_db, og_components, stderr);

    // graph will be changed with extra copies of sequences added
    if (mito_annot)
        parse_organelle_component(asg, annot_db, og_components, min_len, max_copy, ext_m, seq_cf, do_graph_clean, 
                min_cf, min_score, max_eval, bubble_size, tip_size, weak_cross, out_pref, out_opt, OG_MITO, VERBOSE);
    if (pltd_annot)
        parse_organelle_component(asg, annot_db, og_components, min_len, max_copy, ext_p, seq_cf, do_graph_clean, 
                min_cf, min_score, max_eval, bubble_size, tip_size, weak_cross, out_pref, out_opt, OG_PLTD, VERBOSE);

do_clean:
    asg_destroy(asg);
    hmm_annot_db_destroy(annot_db);
    og_component_v_destroy(og_components);

    return ret;
}

#ifdef PATHFINDER_MAIN
#include "ketopt.h"

int VERBOSE = 0;

static ko_longopt_t long_options[] = {
    { "longest",        ko_no_argument,       301 },
    { "circular",       ko_no_argument,       302 },
    { "all",            ko_no_argument,       303 },
    { "edge-c-tag",     ko_required_argument, 304 },
    { "kmer-c-tag",     ko_required_argument, 305 },
    { "seq-c-tag",      ko_required_argument, 306 },
    { "include-trn",    ko_no_argument,       307 },
    { "include-rrn",    ko_no_argument,       308 },
    { "max-bubble",     ko_required_argument, 307 },
    { "max-tip",        ko_required_argument, 310 },
    { "weak-cross",     ko_required_argument, 311 },
    { "no-graph-clean", ko_no_argument,       312 },
    { "mito-annot",     ko_required_argument, 'm' },
    { "pltd-annot",     ko_required_argument, 'p' },
    { "min-score",      ko_required_argument, 's' },
    { "min-gain",       ko_required_argument, 'g' },
    { "min-s-cov",      ko_required_argument, 'q' },
    { "max-copy",       ko_required_argument, 'c' },
    { "max-eval",       ko_required_argument, 'e' },
    { "min-s-len",      ko_required_argument, 'l' },
    { "verbose",        ko_required_argument, 'v' },
    { "version",        ko_no_argument,       'V' },
    { "help",           ko_no_argument,       'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    const char *opt_str = "c:e:f:g:hl:m:o:p:q:s:v:V";
    ketopt_t opt = KETOPT_INIT;
    int c, out_s, out_c, max_copy, ret = 0;
    FILE *fp_help;
    char *out_pref, *mito_annot, *pltd_annot, *ec_tag, *kc_tag, *sc_tag;
    int no_trn, no_rrn, min_len, ext_p, ext_m, bubble_size, tip_size, do_graph_clean;
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
    no_rrn = 1;
    do_graph_clean = 1;
    bubble_size = 100000;
    tip_size = 10000;
    weak_cross = 0.3;
    max_eval = 1e-6;
    min_len = 10000;
    min_score = 300;
    ext_p = 3;
    ext_m = 1;
    seq_cf = .90;
    min_cf = .20;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >=0 ) {
        if (c == 'm') mito_annot = opt.arg;
        else if (c == 'p') pltd_annot = opt.arg;
        else if (c == 'f') seq_cf = atof(opt.arg);
        else if (c == 's') min_score = atof(opt.arg);
        else if (c == 'c') max_copy = atoi(opt.arg);
        else if (c == 'g') {
            char *p;
            ext_p = strtol(opt.arg, &p, 10);
            if (*p == ',') ext_m = strtol(p + 1, &p, 10);
        }
        else if (c == 'q') min_cf = atof(opt.arg);
        else if (c == 'e') max_eval = atof(opt.arg);
        else if (c == 'l') min_len = atoi(opt.arg);
        else if (c == 'o') out_pref = opt.arg;
        else if (c == 301) out_s = 0, ++out_c;
        else if (c == 302) out_s = 1, ++out_c;
        else if (c == 303) out_s = 2, ++out_c;
        else if (c == 304) ec_tag = opt.arg;
        else if (c == 305) kc_tag = opt.arg;
        else if (c == 306) sc_tag = opt.arg;
        else if (c == 307) no_trn = 0;
        else if (c == 308) no_rrn = 1;
        else if (c == 309) bubble_size = atoi(opt.arg);
        else if (c == 310) tip_size = atoi(opt.arg);
        else if (c == 311) weak_cross = atof(opt.arg);
        else if (c == 312) do_graph_clean = 0;
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
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Usage: pathfinder [options] <file>[.gfa[.gz]]\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "  Input/Output:\n");
        fprintf(fp_help, "    -m FILE              mitochondria core gene annotation file [NULL]\n");
        fprintf(fp_help, "    -p FILE              plastid core gene annotation file [NULL]\n");
        fprintf(fp_help, "    -o FILE              prefix of output files [%s]\n", out_pref);
        fprintf(fp_help, "    -f FLOAT             prefer circular path to longest if >= FLOAT sequence covered [%.2f]\n", seq_cf);
        //fprintf(fp_help, "    --longest            output only the longest path [default]\n");
        //fprintf(fp_help, "    --circular           output only the longest circular path\n");
        //fprintf(fp_help, "    --all                output all best paths\n");
        fprintf(fp_help, "    --edge-c-tag STR     edge coverage tag in the GFA file [EC:i] \n");
        fprintf(fp_help, "    --kmer-c-tag STR     kmer coverage tag in the GFA file [KC:i] \n");
        fprintf(fp_help, "    --seq-c-tag  STR     sequence coverage tag in the GFA file [SC:f]\n");
        fprintf(fp_help, "    -v INT               verbose level [%d]\n", VERBOSE);
        fprintf(fp_help, "    --version            show version number\n");
        fprintf(fp_help, "  Classification:\n");
        fprintf(fp_help, "    -s FLOAT             minimum annotation score of a core gene [%.1f]\n", min_score);
        fprintf(fp_help, "    -e FLOAT             maximum E-value of a core gene [%.3e]\n", max_eval);
        fprintf(fp_help, "    -g INT[,INT]         minimum number of core gene gain; the second INT for mitochondria [%d,%d]\n", ext_p, ext_m);
        fprintf(fp_help, "    -l INT               minimum length of a singleton sequence to keep [%d]\n", min_len);
        fprintf(fp_help, "    --include-trn        include tRNA genes for sequence classification\n");
        fprintf(fp_help, "    --include-rrn        include rRNA genes for sequence classification\n");
        fprintf(fp_help, "  Path-finding:\n");
        fprintf(fp_help, "    -q FLOAT             minimum coverage of a sequence compared to the subgraph average [%.2f]\n", min_cf);
        fprintf(fp_help, "    -c INT               maximum copy number to consider [%d]\n", max_copy);
        fprintf(fp_help, "    --max-bubble INT     maximum bubble size for assembly graph clean [%d]\n", bubble_size);
        fprintf(fp_help, "    --max-tip    INT     maximum tip size for assembly graph clean [%d]\n", tip_size);
        fprintf(fp_help, "    --weak-cross FLOAT   maximum relative edge coverage for weak crosslink clean [%.2f]\n", weak_cross);
        fprintf(fp_help, "    --no-graph-clean     do not do assembly graph clean\n");
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Example: ./pathfinder -m asm_annot.mito.txt -p asm_annot.pltd.txt -o oatk.asm asm.gfa\n\n");
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
    
    if (mito_annot == 0 && pltd_annot == 0) {
        fprintf(stderr, "[E::%s] provide at least one HMM profile annotation file (-m and/or -p)\n", __func__);
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
    
    ret = pathfinder(argv[opt.ind], mito_annot, pltd_annot, min_len, ext_p, ext_m, max_copy, 
            max_eval, min_score, min_cf, seq_cf, no_trn, no_rrn, do_graph_clean, bubble_size, tip_size, weak_cross,
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
