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
 * 03/02/23 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "khashl.h"
#include "kstring.h"
#include "kseq.h"
#include "kvec.h"
#include "kdq.h"

#include "path.h"
#include "graph.h"
#include "hmmannot.h"
#include "misc.h"

#undef DEBUG_SRCC
#undef DEBUG_SEG_COV_EST
#undef DEBUG_SEG_COPY
#undef DEBUG_SEG_COPY_EST
#undef DEBUG_SEG_COV_BOUND
#undef DEBUG_SEG_COV_ADJUST
#undef DEBUG_BRUTE_FORCE_OPTIM
#undef DEBUG_SIM_ANNEAL_OPTIM
#undef DEBUG_PATH_FINDER

typedef struct list_node {
    uint32_t v; // sid << 31 | oriv
    struct list_node *prev;
    size_t n_m, n_n;
    struct list_node **next;
} llnode;

typedef llnode* llnodep;

static llnode *new_node(int v) {
    llnode *node;
    MYMALLOC(node, 1);
    node->v = v;
    node->n_m = 4;
    node->n_n = 0;
    node->prev = NULL;
    MYMALLOC(node->next, node->n_m);
    return node;
}

static void add_next(llnode *node, struct list_node *next)
{
    if (node->n_m == node->n_n)
        MYEXPAND(node->next, node->n_m);
    node->next[node->n_n] = next;
    node->n_n++;
}

static void llnode_destroy(llnode *node)
{
    size_t i;
    for (i = 0; i < node->n_n; ++i)
        llnode_destroy(node->next[i]);
    free(node->next);
    free(node);
}

void path_destroy(path_t *path)
{
    if (path->sid)
        free(path->sid);
    if (path->v)
        free(path->v);
}

void path_v_destroy(path_v *path)
{
    size_t i;
    for (i = 0; i < path->n; ++i)
        path_destroy(&path->a[i]);
    free(path->a);
}

static int u64_cmpfunc(const void *a, const void *b)
{
    return (*(uint64_t *) a > *(uint64_t *) b) - (*(uint64_t *) a < *(uint64_t *) b);
}

uint64_t asg_seg_len(asg_t *asg)
{
    asmg_t *g = asg->asmg;
    uint64_t i, seg_len = 0;
    for (i = 0; i < g->n_vtx; ++i) {
        if (g->vtx[i].del) continue;
        seg_len += g->vtx[i].len;
    }
    return seg_len;
}

static double graph_sequence_coverage_lower_bound(asg_t *asg, double cov_nq)
{
    uint32_t i, n_seg, len, cov;
    uint64_t *seqcovs, tot_seg_len, tot_cov, len_thresh;
    double cov_bound;
    asmg_t *g;

    g = asg->asmg;
    n_seg = 0;
    tot_seg_len = 0;
    MYMALLOC(seqcovs, g->n_vtx);
    MYBONE(seqcovs, g->n_vtx);
    for (i = 0; i < g->n_vtx; ++i) {
        if (g->vtx[i].del) continue;
        len = g->vtx[i].len;
        cov = g->vtx[i].cov;
        
        ++n_seg;
        tot_seg_len += len;
        seqcovs[i] = (uint64_t) cov << 32 | len;
    }

    qsort(seqcovs, g->n_vtx, sizeof(uint64_t), u64_cmpfunc);
    len_thresh = tot_seg_len * cov_nq;
    i = 0;
    len = (uint32_t) seqcovs[i];
    tot_seg_len = tot_cov = 0;
    while (tot_seg_len + len <= len_thresh) {
        tot_cov += (seqcovs[i] >> 32) * len;
        tot_seg_len += len;
        len = (uint32_t) seqcovs[++i];
    }
    if (tot_seg_len < len_thresh)
        tot_cov += (seqcovs[i] >> 32) * (len_thresh - tot_seg_len);
    cov_bound = (double) tot_cov / len_thresh;
    free(seqcovs);

    cov_bound *= 1 - cov_nq;

#ifdef DEBUG_SEG_COV_BOUND
    fprintf(stderr, "[DEBUG_SEG_COV_BOUND::%s] estimated sequence coverage lower boundary: %.3f\n", __func__, cov_bound);
#endif

    return cov_bound;
}

/*** old implementation
 *
static double graph_sequence_coverage_rough(asg_t *asg)
{
    uint32_t i, len, cov, deg_in, deg_out;
    double tot_seg_len, tot_cov, avg_cov;
    asmg_t *g;

    g = asg->asmg;
    tot_seg_len = tot_cov = 0;
    for (i = 0; i < g->n_vtx; ++i) {
        if (g->vtx[i].del) continue;
        len = g->vtx[i].len;
        cov = g->vtx[i].cov;
        avg_cov = cov;
        deg_out = MAX(1, asmg_arc_n1(g, i<<1));
        deg_in  = MAX(1, asmg_arc_n1(g, i<<1|1));
        avg_cov = avg_cov * 2 / (deg_out + deg_in);
#ifdef DEBUG_SEG_COV_EST
        // seg and vtx indices are always interchangeable as we never do graph clean
        fprintf(stderr, "[DEBUG_SEG_COV_EST::%s] %s %u %u [adj out deg: %u] [adj in deg: %u] [normalized: %.3f]\n",
                __func__, asg->seg[i].name, len, cov, deg_out, deg_in, avg_cov);
#endif
        tot_seg_len += len;
        tot_cov += avg_cov * len;
    }
    
    if (tot_seg_len == 0) return 0;

    avg_cov = tot_cov / tot_seg_len;

#ifdef DEBUG_SEG_COV_EST
    fprintf(stderr, "[DEBUG_SEG_COV_EST::%s] estimated sequence coverage: %.3f\n", __func__, avg_cov);
#endif

    return avg_cov;
}
**/

static double graph_sequence_coverage_rough(asg_t *asg, double min_cf)
{
    uint32_t i, j, len, cov, best1;
    double tot_len, tot_len_c, tot_rm, avg_cov, near1, diff1;
    asmg_t *g;
    kvec_t(uint64_t) lc_p;

    g = asg->asmg;
    kv_init(lc_p);
    for (i = 0; i < g->n_vtx; ++i) {
        if (g->vtx[i].del) continue;
        kv_push(uint64_t, lc_p, (uint64_t) (g->vtx[i].cov) << 32 | (g->vtx[i].len));
    }
    if (lc_p.n == 0) return .0;

    qsort(lc_p.a, lc_p.n, sizeof(uint64_t), u64_cmpfunc);

    best1 = 0;
    near1 = DBL_MAX;
    for (i = 0; i < lc_p.n; ++i) {
        avg_cov = (double) (lc_p.a[i] >> 32);
        if (avg_cov == 0)
            continue;
        tot_len = tot_len_c = tot_rm = 0;
        for (j = 0; j < lc_p.n; ++j) {
            len = (uint32_t) lc_p.a[j];
            cov = lc_p.a[j] >> 32;
            if (cov / avg_cov >= min_cf) {
                tot_len += len;
                tot_len_c += (double) len * cov / avg_cov;
            } else {
                tot_rm += len;
            }
        }

        // TODO is this good?
        if (tot_rm / (tot_rm + tot_len) > 0.7) break;

        if (tot_len > 0) {
            diff1 = fabs(tot_len_c / tot_len - 1.0);
            if (diff1 < near1) {
                near1 = diff1;
                best1 = i;
            }
#ifdef DEBUG_SEG_COV_EST
            fprintf(stderr, "[DEBUG_SEG_COV_EST::%s] seg %u [%u %lu] - len: %.0f; len_rm: %.0f; len_c: %.3f; diff1: %.3f; best: [%u %.3f]\n",
                    __func__, i, (uint32_t) lc_p.a[i], lc_p.a[i]>>32, tot_len, tot_rm, tot_len_c, diff1, best1, near1);
#endif
        }
    }
    if (near1 == DBL_MAX) {
        kv_destroy(lc_p);
        return .0;
    }

    avg_cov = (double) (lc_p.a[best1] >> 32);

#ifdef DEBUG_SEG_COV_EST
    fprintf(stderr, "[DEBUG_SEG_COV_EST::%s] estimated sequence coverage: %.3f\n", __func__, avg_cov);
#endif

    kv_destroy(lc_p);
    return avg_cov;
}

static void make_seg_dups(asg_t *asg, kh_u32_t *seg_dups, uint32_t s, uint32_t copy)
{
    // copy number include sequence itself
    // i.e., make 'copy' copies and remove the original
    // FIXME arcs involving self cycles are not copied since it will potentially cause problems
    // for example
    // in a graph containing path (X+)<->(a+,b+,b+,c+)<->(Y+)
    // the self cycle b+,b+ will increase the ambiguities of the path selection
    // valid path including
    // (X+)->(a+,b+,b+,c+)->(Y+)->(a-,b-,b-,c-)
    // (X+)->(a+,b+,b+,b+,c+)->(Y+)->(a-,b-,c-)
    // TODO solve the copy number of self cycles in postprocessing
    // now the self arcs is processed in this way
    // (X+)->(a+)[4]->(Y+) => (X+)->(a0+) ->(a1+) ->(a2+) ->(a3+) ->(Y+)
    //                            ->(a1+) ->(a2+) ->(a0+) ->(a1+) ->(a0+)
    //                            ->(a2+) ->(a3+) ->(a3+) ->(a2+) ->(a1+)
    //                            ->(a3+) ->(Y+)  ->(Y+)  ->(Y+)  ->(a2+)
    // (a+)->(a-) arcs are not copied
    uint32_t sid;
    uint64_t i, j, v, nv, pv, self_arc;
    kvec_t(uint64_t) arcs_diff;
    asmg_t *g;
    asmg_arc_t *av;
    asmg_vtx_t *vt;
    char *name;
    asg_seg_t *seg, *seg_c;
    khint32_t k;
    int absent;

    g = asg->asmg;
    // mark arcs from vtx
    kv_init(arcs_diff);
    self_arc = UINT64_MAX;

    for (i = 0; i < 2; ++i) {
        v = s<<1 | i;
        pv = g->idx_p[v];
        nv = asmg_arc_n(g, v);
        av = asmg_arc_a(g, v);
        for (j = 0; j < nv; ++j) {
            if (av[j].del) continue;
            if ((av[j].v >> 1) != (av[j].w >> 1))
                kv_push(uint64_t, arcs_diff, pv+j);
            else if (av[j].v == av[j].w && i == 0)
                // (a+)->(a+)
                // avoid adding v+->v+ and v-->v- twice
                self_arc = pv + j;
        }
    }

    for (i = 0; i < copy; ++i) {
        // make a copy of the segment
        seg = &asg->seg[s];
        MYMALLOC(name, strlen(seg->name) + floor(log10(abs(i + 1))) + 7);
        sprintf(name, "%s_copy%lu", seg->name, i);
        sid = asg_add_seg(asg, name, 0);
        free(name);
        // add to seg dups map
        k = kh_u32_put(seg_dups, sid, &absent);
        kh_val(seg_dups, k) = s;
        // only copy essential fields
        seg = &asg->seg[s]; // need to redo this as the seg might be reallocated
        seg_c = &asg->seg[sid];
        seg_c->len = seg->len;
        seg_c->cov = seg->cov;

        // add vtx
        asmg_vtx_addp(g, &vt);
        MYBZERO(vt, 1);
        vt->len = seg->len;
        // using vtx coverage instead of seg coverage
        //vt->cov = seg->cov / copy;
        vt->cov = g->vtx[s].cov / copy;

        // make copies of links
        // arc cov already changed
        for (j = 0; j < arcs_diff.n; ++j) {
            av = &g->arc[arcs_diff.a[j]];
            asmg_arc_add2(g, sid << 1 | (av->v&1), av->w, av->ln, av->ls, UINT64_MAX, av->cov / copy, av->comp);
        }
        
        if (self_arc != UINT64_MAX) {
            // tandem repeat
            // self arcs are added between copies
            av = &g->arc[self_arc];
            for (j = 0; j < i; ++j) {
                asmg_arc_add2(g, (sid - i + j) << 1, sid << 1, av->ln, av->ls, UINT64_MAX, av->cov / copy, 0);
                asmg_arc_add2(g, sid << 1, (sid - i + j) << 1, av->ln, av->ls, UINT64_MAX, av->cov / copy, 0);
            }
        }
    }

    kv_destroy(arcs_diff);

    // need to redo sort and index but not graph clean
    asmg_finalize(g, 0);

    // delete the original vtx
    asmg_vtx_del(g, s, 1);

    return;
}

#define EM_MAX_ITER 1000

double graph_sequence_coverage_precise(asg_t *asg, double min_cf, int min_copy, int max_copy, int **_copy_number)
{
    uint32_t i, n_seg, iter;
    int *copy_number;
    double total_covs, total_lens, avg_cov, new_avg_cov, min_avg_cov;
    asmg_t *g;
    
    n_seg = asg->n_seg;
    g = asg->asmg;
    min_avg_cov = graph_sequence_coverage_lower_bound(asg, 0.3);
    avg_cov = graph_sequence_coverage_rough(asg, min_cf);
#ifdef DEBUG_SEG_COPY_EST
    fprintf(stderr, "[DEBUG_SEG_COPY_EST::%s] min coverage: %.3f; rough avg coverage: %.3f\n",
            __func__, min_avg_cov, avg_cov);
#endif
    avg_cov = MAX(avg_cov, min_avg_cov);
    MYCALLOC(copy_number, n_seg);
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        copy_number[i] = MIN(MAX(min_copy, lround((double) g->vtx[i].cov / avg_cov)), max_copy);
    }

    iter = 0;
    while (iter++ < EM_MAX_ITER) {
#ifdef DEBUG_SEG_COPY_EST
        fprintf(stderr, "[DEBUG_SEG_COPY_EST::%s] iteration %u avg coverage: %.3f\n", __func__, iter, avg_cov);
#endif
        total_covs = total_lens = 0;
        for (i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del) continue;
            total_lens += (double) g->vtx[i].len * copy_number[i];
            total_covs += (double) g->vtx[i].len * g->vtx[i].cov;
        }
        // FIXME total_lens could be zero
        new_avg_cov = total_lens < FLT_EPSILON? DBL_MAX : (total_covs / total_lens);
        new_avg_cov = MAX(new_avg_cov, min_avg_cov);
        if (fabs(new_avg_cov - avg_cov) < FLT_EPSILON) 
            break; // converged
        avg_cov = new_avg_cov;
        for (i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del) continue;
            copy_number[i] = MIN(MAX(min_copy, lround((double) g->vtx[i].cov / avg_cov)), max_copy);
        }
    }

#ifdef DEBUG_SEG_COPY_EST
    fprintf(stderr, "[DEBUG_SEG_COPY_EST::%s] sequence copy number estimation finished in %u iterations with an average sequence coverage %.3f\n",
            __func__, iter, avg_cov);
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        fprintf(stderr, "[DEBUG_SEG_COPY_EST::%s] %s %u %d\n", __func__, asg->seg[i].name, g->vtx[i].cov, copy_number[i]);
    }
#endif
    
    if (_copy_number)
        *_copy_number = copy_number;
    else
        free(copy_number);

    return avg_cov;
}

#define DVAL(var) ((var)->D)
// with BALANCE_IN_OUT defined
// the objective function considers balanced indegree and outdegree
#define BALANCE_IN_OUT
#ifdef BALANCE_IN_OUT
#define FVAL(fun) do { \
    int __i, __n = (fun)->N; \
    double __val[2] = {.0, .0}; \
    for (__i = 0; __i < __n; ++__i) \
        __val[(fun)->V[__i]&1] += DVAL((fun)->VAR[(fun)->V[__i]>>1]); \
    (fun)->VAL = (fun)->weight * (fabs((fun)->v_exp - __val[0]) / 2 + \
            fabs((fun)->v_exp - __val[1]) / 2 + \
            fabs(__val[0] - __val[1])); \
} while (0)
#else
#define FVAL(fun) do { \
    int __i, __n = (fun)->N; \
    double __val = (fun)->v_exp; \
    for (__i = 0; __i < __n; ++__i) \
        __val -= DVAL((fun)->VAR[(fun)->V[__i]]); \
    (fun)->VAL = (fun)->weight * __val * __val; \
} while (0)
#endif

typedef struct var {
    int B, D;
    struct var *prev, *next;
} var_t;

typedef struct {
    double weight; // weight
    double v_exp; // expected value
    var_t **VAR;
    int N, *V;
    double VAL;
} func_t;

#define BRUTE_FORCE_N_LIM 100000000

static inline double fvals(func_t *funcs, int n_func)
{
    int i;
    double fval = 0.;
    for (i = 0; i < n_func; ++i) {
        FVAL(&funcs[i]);
        fval += funcs[i].VAL;
    }
    return fval;
}

static inline void copy_results(var_t **vars, int n_var, int *res)
{
    int i;
    for (i = 0; i < n_var; ++i)
        res[i] = vars[i]->D;
}

// brute force optimization
static void estimate_arc_copy_number_brute_force_impl(func_t *funcs, int n_func, var_t **vars, int n_var, int *res, int64_t sol_space_size)
{
    double fval, m_fval;
    int64_t sol;
    int a, v;
    fval = fvals(funcs, n_func);
    m_fval = fval;
    copy_results(vars, n_var, res);
    sol = 0;
    while (++sol < sol_space_size) {
        a = 1, v = 0;
        while (a) {
            vars[v] = vars[v]->next;
            a = !vars[v]->B;
            ++v;
        }
        fval = fvals(funcs, n_func);
        if (fval < m_fval) {
            m_fval = fval;
            copy_results(vars, n_var, res);
        }
        if (fabs(m_fval) < FLT_EPSILON)
            break;
    }
#ifdef DEBUG_BRUTE_FORCE_OPTIM
    fprintf(stderr, "[DEBUG_BRUTE_FORCE_OPTIM::%s] brute force search finished after %ld/%ld attempts with a minimum fval: %.6f\n",
            __func__, sol, sol_space_size, m_fval);
#endif
}

#define SA_TEMPERATURE  1000
#define SA_COOLING_RATE .999
#define SA_MAX_ATTEMPTS 100
#define SA_RESTART_TEMP .99

static inline void random_walk_to_neighbour(var_t **var)
{
    if (rand() < RAND_MAX>>1) {
        // move to prev
        *var = (*var)->B == 0? (*var)->next : (*var)->prev;
    } else {
        // move to next
        *var = (*var)->next->B == 0? (*var)->prev : (*var)->next;
    }
}

static void set_vars(var_t **vars, int n_var, int *res)
{
    int i;
    for (i = 0; i < n_var; ++i) {
        while (vars[i]->D != res[i])
            vars[i] = vars[i]->next;
    }
}

// simulated annealing optimization
static void estimate_arc_copy_number_siman_impl(func_t *funcs, int n_func, var_t **vars, int n_var, int *res)
{
    int i, iter, n_iter;
    double optim_cost, current_cost, new_cost, temp0, temp, p;
    var_t *var;

    current_cost = fvals(funcs, n_func);
    optim_cost = current_cost;
    copy_results(vars, n_var, res);
    
    srand(1234);
    temp0 = SA_TEMPERATURE;
    n_iter = 0;
    for (iter = 0; iter < SA_MAX_ATTEMPTS; ++iter) {
        temp = temp0;
        while (temp > 1e-6) {
            // random select a var to update 
            i = rand() % n_var;
            // record the old var
            var = vars[i];
            // take a random walk to either prev or next
            random_walk_to_neighbour(&vars[i]);
            // calculate new cost
            new_cost = fvals(funcs, n_func);
            // record optim solution
            if (new_cost < optim_cost) {
                optim_cost = new_cost;
                copy_results(vars, n_var, res);
            }
            // acceptance probability
            p = exp(-(new_cost - current_cost) / temp);
            if (new_cost < current_cost || (double) rand() / RAND_MAX < p) {
                // accept update and update cost
                current_cost = new_cost;
            } else {
                // rollback to reject the update
                vars[i] = var;
            }
            // cooling down
            temp *= SA_COOLING_RATE;
            ++n_iter;
        }
        if (optim_cost == 0) break;
        // continue searching from the best solution so far
        temp0 *= SA_RESTART_TEMP;
        set_vars(vars, n_var, res);
#ifdef DEBUG_SIM_ANNEAL_OPTIM
        fprintf(stderr, "[DEBUG_SIM_ANNEAL_OPTIM::%s] optimization after %d iteration [%d attempts] with a minimum fval: %.6f\n",
                __func__, iter, n_iter, optim_cost);
#endif

    }
#ifdef DEBUG_SIM_ANNEAL_OPTIM
    fprintf(stderr, "[DEBUG_SIM_ANNEAL_OPTIM::%s] simulated annealing search finished after %d attempts with a minimum fval: %.6f\n",
            __func__, n_iter, optim_cost);
#endif
}

int adjust_sequence_copy_number_by_graph_layout(asg_t *asg, double seq_coverage, double *_adjusted_cov, int *copy_number, int max_copy, int max_round)
{
    uint32_t i, j, n_seg, n_group, a_g, *arc_group;
    uint64_t link_id;
    int updated, round;
    asmg_t *g;
    asmg_arc_t *a;

    if (_adjusted_cov)
        *_adjusted_cov = seq_coverage;
    if (max_round == 0)
        max_round = 1;

    updated = 0;

    g = asg->asmg;
    n_seg = asg->n_seg;
    arc_group = asmg_uext_arc_group(g, &n_group);

#ifdef DEBUG_SEG_COV_ADJUST
    fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] identified %u arc groups\n", __func__, n_group);
    for (i = 0; i < g->n_arc; ++i) {
        a = &g->arc[i];
        if (a->del) continue;
        link_id = a->link_id;
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] %s%c -> %s%c arc group %u\n", __func__,
                asg->seg[a->v>>1].name, "+-"[a->v&1],
                asg->seg[a->w>>1].name, "+-"[a->w&1],
                arc_group[link_id]);
    }
#endif

    if (n_group == 0) {
        free(arc_group);
        return 0;
    }

#ifdef DEBUG_SEG_COV_ADJUST
    fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] sequence copy number BEFORE adjusted by graph layout\n", __func__);
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] %s %lu %u %d\n", __func__, 
                asg->seg[i].name, g->vtx[i].len, g->vtx[i].cov, copy_number[i]);
    }
#endif

    uint32_t *arc_copy_lb, *arc_copy_ub, vlb, wlb, lb, ub;
    MYCALLOC(arc_copy_lb, n_group); // arc copy number lower bound
    MYCALLOC(arc_copy_ub, n_group); // arc copy number upper bound

    // calculate lower and upper boundary of arc copy number for each arc group
    for (i = 0; i < g->n_arc; ++i) {
        a = &g->arc[i];
        if (a->del) continue;
        link_id = a->link_id;
        a_g = arc_group[link_id];
        // if a is the only outgoing/incoming edge from/to v/w
        // then lower bound is the v/w copy number
        // otherwise 0
        // upper bound is always the v/w copy number
        vlb = asmg_arc_n1(g, a->v) == 1? copy_number[a->v>>1] : 0;
        wlb = asmg_arc_n1(g, a->w^1) == 1? copy_number[a->w>>1] : 0;
        lb = MIN(vlb, wlb);
        ub = MAX(copy_number[a->v>>1], copy_number[a->w>>1]);
        // relax boundary by a factor of 1/3 - at least one copy
        lb = (uint32_t) ((double) lb * 2 / 3);
        ub = (uint32_t) ((double) ub * 4 / 3) + 1;
        ub = MIN(ub, (uint32_t) max_copy);
        arc_copy_lb[a_g] = MIN(lb, arc_copy_lb[a_g]);
        arc_copy_ub[a_g] = MAX(ub, arc_copy_ub[a_g]);
    }

    // make variable list
    var_t **VAR;
    MYCALLOC(VAR, n_group);
    for (i = 0; i < n_group; ++i) {
        var_t *var0;
        MYCALLOC(var0, 1);
        lb = arc_copy_lb[i];
        ub = arc_copy_ub[i];
        var0->B = 0;
        var0->D = lb;
        VAR[i] = var0;
        for (j = 1; j <= ub - lb; ++j) {
            var_t *var1;
            MYCALLOC(var1, 1);
            var1->B = j;
            var1->D = lb + j;
            var0->next = var1;
            var1->prev = var0;
            var0 = var1;
        }
        var0->next = VAR[i];
        VAR[i]->prev = var0;
    }

    // make objective function list
    kvec_t(func_t) funcs;
    kv_init(funcs);
    func_t *FUN;
    uint32_t k, v, na;
    int *funcmap; // the objective function index for each seg
    asmg_arc_t *av;
    MYMALLOC(funcmap, n_seg);
    for (i = 0; i < n_seg; ++i) {
        funcmap[i] = -1;
        if (g->vtx[i].del)
            continue;
#ifdef BALANCE_IN_OUT
        // with BALNCE_IN_OUT defined
        // the incoming and outgoing arcs of a vertex
        // are put in the same objective function unit
        // and marked by the last bit of V
        kvec_t(int) V = {0, 0 ,0};
        for (k = 0; k < 2; ++k) {
            v = i << 1 | k;
            na = asmg_arc_n(g, v);
            av = asmg_arc_a(g, v);
            for (j = 0; j < na; ++j) {
                if (av[j].del)
                    continue;
                link_id = av[j].link_id;
                a_g = arc_group[link_id];
                assert(a_g != (uint32_t) -1);
                kv_push(int, V, (int) (a_g << 1 | k));
            }
        }
        if (V.n) {
            MYREALLOC(V.a, V.n);
            funcmap[i] = funcs.n;
            kv_pushp(func_t, funcs, &FUN);
            FUN->weight = log10(g->vtx[i].len);
            FUN->v_exp = g->vtx[i].cov / seq_coverage; // copy_number[i];
            FUN->VAR = VAR;
            FUN->N = V.n;
            FUN->V = V.a;
            FUN->VAL = 0.;
        }
#else
        for (k = 0; k < 2; ++k) {
            v = i << 1 | k;
            na = asmg_arc_n(g, v);
            av = asmg_arc_a(g, v);
            kvec_t(int) V = {0, 0 ,0};
            for (j = 0; j < na; ++j) {
                if (av[j].del)
                    continue;
                link_id = av[j].link_id;
                a_g = arc_group[link_id];
                assert(a_g != (uint32_t) -1);
                kv_push(int, V, (int) a_g);
            }
            if (V.n) {
                MYREALLOC(V.a, V.n);
                funcmap[i] = funcs.n;
                kv_pushp(func_t, funcs, &FUN);
                FUN->weight = log10(g->vtx[i].len);
                FUN->v_exp = g->vtx[i].cov / seq_coverage; // copy_number[i];
                FUN->VAR = VAR;
                FUN->N = V.n;
                FUN->V = V.a;
                FUN->VAL = 0.;
            }
        }
#endif
    }

    // do optimization
    int *arc_copy;
    double adjusted_cov, min_avg_cov;
    int64_t sol_space_size;
    
    min_avg_cov = graph_sequence_coverage_lower_bound(asg, 0.3);
    adjusted_cov = seq_coverage;

    MYCALLOC(arc_copy, n_group);
    sol_space_size = 1;
    for (i = 0; i < n_group && sol_space_size <= BRUTE_FORCE_N_LIM; ++i)
        sol_space_size *= (arc_copy_ub[i] - arc_copy_lb[i] + 1);
    
    round = 0;
    while (round++ < max_round) {
#ifdef DEBUG_SEG_COV_ADJUST
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] adjusting copy number round %d of %d\n", __func__, round, max_round);
#endif
        if (sol_space_size <= BRUTE_FORCE_N_LIM) {
            // do brute force optimization
#ifdef DEBUG_SEG_COV_ADJUST
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] run brute force searching: %ld\n",
                    __func__, sol_space_size);
#endif
            estimate_arc_copy_number_brute_force_impl(funcs.a, funcs.n, VAR, n_group, arc_copy, sol_space_size);
        } else {
            // do simulated annealing optimization
#ifdef DEBUG_SEG_COV_ADJUST
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] run simulated annealing optimization: %ld\n",
                    __func__, sol_space_size);
#endif
            estimate_arc_copy_number_siman_impl(funcs.a, funcs.n, VAR, n_group, arc_copy);
        }

#ifdef DEBUG_SEG_COV_ADJUST
        for (i = 0; i < n_group; ++i)
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] arc group %u optimum copy number %d [%u %u]\n",
                    __func__, i, arc_copy[i], arc_copy_lb[i], arc_copy_ub[i]);
#endif
        // update average sequence coverage
        double total_covs, total_lens, new_adjusted_cov, copies;
        total_covs = total_lens = 0;
        for (i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del)
                continue;
            copies = 0;
            for (k = 0; k < 2; ++k) {
                v = i << 1 | k;
                na = asmg_arc_n(g, v);
                av = asmg_arc_a(g, v);
                for (j = 0; j < na; ++j) {
                    if (av[j].del)
                        continue;
                    link_id = av[j].link_id;
                    a_g = arc_group[link_id];
                    assert(a_g != (uint32_t) -1);
                    copies += arc_copy[a_g];
                }
            }
            total_lens += (double) g->vtx[i].len * copies / 2;
            total_covs += (double) g->vtx[i].len * g->vtx[i].cov;
        }
    
        if (total_lens < FLT_EPSILON) {
#ifdef DEBUG_SEG_COV_ADJUST
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] zero copies for all sequences\n", __func__);
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] adjusting copy number terminated in round %d of %d\n",
                    __func__, round, max_round);
#endif
            goto do_clean;
        }

        new_adjusted_cov = total_covs / total_lens;
        new_adjusted_cov = MAX(new_adjusted_cov, min_avg_cov);

        if (fabs(new_adjusted_cov - adjusted_cov) < FLT_EPSILON) {
#ifdef DEBUG_SEG_COV_ADJUST
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] adjusting copy number converaged in round %d of %d\n",
                    __func__, round, max_round);
#endif
            break; // converaged
        }

        // update adjusted sequence coverage
        adjusted_cov = new_adjusted_cov;
        // update objective functions
        // more precisely the expected copy number
        for (i = 0; i < n_seg; ++i) {
            if (funcmap[i] == -1)
                continue;
            funcs.a[funcmap[i]].v_exp = g->vtx[i].cov / adjusted_cov;
            funcs.a[funcmap[i]].VAL = 0.;
        }
        // reset VARs
        for (i = 0; i < n_group; ++i) {
            var_t *var = VAR[i];
            while (var->B) var = var->next;
            VAR[i] = var;
        }

#ifdef DEBUG_SEG_COV_ADJUST
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] adjusting copy number round %d of %d DONE with adjusted coverage: %.3f\n",
                __func__, round, max_round, adjusted_cov);
#endif
    }

    // update sequence copy number using the arc copy number information
    int new_copy[2];
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del)
            continue;
#ifdef DEBUG_SEG_COV_ADJUST
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] #### sequence %s ####\n", __func__, asg->seg[i].name);
#endif
        for (k = 0; k < 2; ++k) {
            v = i << 1 | k;
            na = asmg_arc_n(g, v);
            av = asmg_arc_a(g, v);
            new_copy[k] = 0;
#ifdef DEBUG_SEG_COV_ADJUST
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] #### %sing arcs\n", __func__, k? "incom" : "outgo");
#endif
            for (j = 0; j < na; ++j) {
                if (av[j].del)
                    continue;
                link_id = av[j].link_id;
                a_g = arc_group[link_id];
                assert(a_g != (uint32_t) -1);
                new_copy[k] += arc_copy[a_g];
#ifdef DEBUG_SEG_COV_ADJUST
                fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] %s%c -> %s%c: %d\n", __func__, 
                        asg->seg[av[j].v>>1].name, "+-"[av[j].v&1],
                        asg->seg[av[j].w>>1].name, "+-"[av[j].w&1],
                        arc_copy[a_g]);
#endif
            }
        }
#ifdef DEBUG_SEG_COV_ADJUST
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] indegree: %d;  outdegree: %d; diff: %d\n", __func__, 
                new_copy[1], new_copy[0], abs(new_copy[1] - new_copy[0]));
#endif
        // only update copy number if indegree matches outdegree
        // TODO better strategy to update sequence copy number
        if (new_copy[0] == new_copy[1] && copy_number[i] != new_copy[0]) {
#ifdef DEBUG_SEG_COV_ADJUST
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] sequence %s copy number updated %d -> %d\n",
                    __func__, asg->seg[i].name, copy_number[i], new_copy[0]);
#endif
            copy_number[i] = new_copy[0];
            updated = 1;
        }
    }
#ifdef DEBUG_SEG_COV_ADJUST
    fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] sequence copy number AFTER adjusted by graph layout\n", __func__);
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] %s %lu %u %d\n", __func__, 
                asg->seg[i].name, g->vtx[i].len, g->vtx[i].cov, copy_number[i]);
    }
#endif
    
    if (_adjusted_cov)
        *_adjusted_cov = adjusted_cov;

do_clean:
    for (i = 0; i < n_group; ++i) {
        var_t *var0, *var1;
        var0 = VAR[i];
        while (var0->B != 0)
            var0 = var0->next;
        var0 = var0->next;
        v = 0;
        while (!v) {
            var1 = var0->next;
            v = var0->B == 0;
            free(var0);
            var0 = var1;
        }
    }
    free(VAR);
    for (i = 0; i < funcs.n; ++i)
        free(funcs.a[i].V);
    free(funcs.a);
    free(funcmap);
    free(arc_group);
    free(arc_copy_lb);
    free(arc_copy_ub);
    free(arc_copy);

    return updated;
}

kh_u32_t *sequence_duplication_by_copy_number(asg_t *asg, int *copy_number, int allow_del)
{
    uint32_t i, n_seg;
    int copy;
    kh_u32_t *seg_dups;
    asmg_t *g;

    n_seg = asg->n_seg;
    g = asg->asmg;
    seg_dups = kh_u32_init();
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        copy = copy_number[i];
        if (copy > 1) {
            // do copy
            make_seg_dups(asg, seg_dups, i, copy);
#ifdef DEBUG_SEG_COPY
            fprintf(stderr, "[DEBUG_SEG_COPY::%s] make %d extra cop%s of %s [%lu %u]\n", __func__, 
                    copy, copy > 1? "ies" : "y", asg->seg[i].name, g->vtx[i].len, g->vtx[i].cov);
#endif
        } else if (copy == 0 && allow_del) {
            /*** TODO is this safe? **/
            asmg_vtx_del(g, i, 1);
#ifdef DEBUG_SEG_COPY
            fprintf(stderr, "[DEBUG_SEG_COPY::%s] delete seg %s [%lu %u]\n", __func__,
                    asg->seg[i].name, g->vtx[i].len, g->vtx[i].cov);
#endif
        }
    }

#ifdef DEBUG_SEG_COPY
    fprintf(stderr, "[DEBUG_SEG_COPY::%s] graph after expansion\n", __func__);
    asg_print(asg, stderr, 1);
    fprintf(stderr, "[DEBUG_SEG_COPY::%s] expanded graph stats\n", __func__);
    asg_stat(asg, stderr);
#endif
    return seg_dups;
}

KDQ_INIT(llnodep)

static int path_contain_vertex(llnode *node, uint32_t v)
{
    int contained = 0;
    while (node != NULL) {
        if ((node->v >> 1) == (v >> 1)) {
            contained = 1;
            break;
        }
        node = node->prev;
    }

    return contained;
}

#define PATH_N_LIM 1000000

static llnode **graph_path_extension(asmg_t *g, llnode *root, kh_u32_t *seg_dups, uint64_t *n_leaf, int *exceed_limit_)
{
    uint64_t i, j, v, w, nv;
    int skip, exceed_limit;
    llnode *next_node, *node;
    kdq_llnodep_t *node_q;
    void *km;
    asmg_arc_t *av;
    khint32_t k32;

    kvec_t(llnodep) leaf_node;
    kvec_t(uint32_t) dups;

    kv_init(leaf_node);
    kv_init(dups);
    km = km_init();
    node_q = kdq_init_llnodep(km);

    exceed_limit = 0;
    *kdq_pushp_llnodep(node_q) = root;

    while (kdq_size(node_q) > 0) {
        node = *kdq_shift_llnodep(node_q);
        v = node->v;
        nv = asmg_arc_n(g, v);
        av = asmg_arc_a(g, v);
        dups.n = 0;
        for (i = 0; i < nv; ++i) {
            if (av[i].del) continue;

            w = av[i].w;
            skip = 0;
            k32 = kh_u32_get(seg_dups, w >> 1);
            if (k32 < kh_end(seg_dups)) {
                // dup seg
                // check if already processed
                for (j = 0; j < dups.n; ++j) {
                    if (dups.a[j] == kh_val(seg_dups, k32)) {
                        skip = 1;
                        break;
                    }
                }
            }

            if (!skip && !path_contain_vertex(node, w)) {
                next_node = new_node(w);
                next_node->prev = node;
                add_next(node, next_node);
                *kdq_pushp_llnodep(node_q) = next_node;
                
                if (k32 < kh_end(seg_dups))
                    kv_push(uint32_t, dups, kh_val(seg_dups, k32));
            }
        }

        if (node->n_n == 0) {
            // leaf node
            kv_push(llnodep, leaf_node, node);
        }
        
        if (kdq_size(node_q) + leaf_node.n > PATH_N_LIM) {
            exceed_limit = 1;
            break;
        }
    }

    kdq_destroy_llnodep(node_q);
    km_destroy(km);
    kv_destroy(dups);

    if (exceed_limit) {
        *n_leaf = 0;
        *exceed_limit_ = 1;
        kv_destroy(leaf_node);
        return 0;
    } else {
        *n_leaf = leaf_node.n;
        *exceed_limit_ = 0;
        return leaf_node.a;
    }
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

static inline void rev_path(path_t *path)
{
    size_t i;
    rev_array(path->v, path->nv);
    for (i = 0; i < path->nv; ++i)
        path->v[i] ^= 1;
}

static uint64_t find_source_vtx(asmg_t *g, int use_max_scc)
{
    uint64_t i, s, n_seg;
    s = UINT64_MAX;
    
    if (!use_max_scc) {
        // naively select the longest sequence
        uint64_t m_len, len;
        n_seg = g->n_vtx;
        m_len = 0;
        for (i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del) continue;
            len = g->vtx[i].len * g->vtx[i].cov;
            if (m_len < len) {
                m_len = len;
                s = i;
            }
        }
    } else {
        // select the largest sequence from the largest SCC
        int c, m_c, n_scc, *scc;
        uint64_t m_len, len, *lens;
        n_seg = g->n_vtx * 2;

        MYMALLOC(scc, n_seg);
        // find strongly connected components
        n_scc = asmg_tarjans_scc(g, scc);
        
        MYCALLOC(lens, n_scc);
        // find larges SCC
        for (i = 0; i < n_seg; ++i) {
            if (scc[i] < 0) continue; // deleted vtx
            if (scc[i] != scc[i^1] || (i & 1))
                lens[scc[i]] += g->vtx[i>>1].len * g->vtx[i>>1].cov;
        }
        m_len = 0;
        m_c = -1;
        for (c = 0; c < n_scc; ++c) {
            if (m_len < lens[c]) {
                m_len = lens[c];
                m_c = c;
            }
        }
        if (m_c >= 0) {
            m_len = 0;
            for (i = 0; i < n_seg; ++i) {
                if (scc[i] != m_c) continue;
                len = g->vtx[i>>1].len * g->vtx[i>>1].cov;
                if (m_len < len) {
                    m_len = len;
                    s = i;
                }
            }
            s >>= 1;
        }
        free(scc);
        free(lens);
    }

    return s;
}

void graph_path_finder(asg_t *asg, kh_u32_t *seg_dups, path_v *paths, double sub_circ_minf, int is_pltd)
{
    uint64_t i, j, s, n_leaf;
    int circ, exceed_limit;
    asmg_t *g;
    llnode *root, *next_node, *node, **leaves;
    kvec_t(llnodep) leaf_node, root_node;

    g = asg->asmg;
    
    // find source vertex
    s = find_source_vtx(g, 1);

    if (s == UINT64_MAX) return;

    kv_init(root_node);
    kv_init(leaf_node);

    root = new_node(s << 1);
    kv_push(llnodep, root_node, root);

    n_leaf = 0;
    exceed_limit = 0;
    leaves = graph_path_extension(g, root, seg_dups, &n_leaf, &exceed_limit);
    // for linear paths do extension from the other direction of root node
    // no - should do this even if the path is circular
    // n_leaf will be zero if path number exceeds limit
    for (i = 0; i < n_leaf; ++i) {
        node = leaves[i];
        // circ = asmg_arc1(g, node->v, s << 1) != 0;
        // if (circ) {
        //     kv_push(llnodep, leaf_node, node);
        //     continue;
        // }
        
        // make a new path tracing back from leaf to root
        root = new_node(node->v^1);
        kv_push(llnodep, root_node, root);
        while (node->prev != NULL) {
            next_node = new_node(node->prev->v^1);
            next_node->prev = root;
            add_next(root, next_node);
            root = next_node;
            node = node->prev;
        }

        assert(root->v == (s<<1 | 1)); // should be always ended with seq s
        uint64_t n = 0;
        exceed_limit = 0;
        llnode **tmp_nodes = graph_path_extension(g, root, seg_dups, &n, &exceed_limit);
        
        for (j = 0; j < n; ++j)
            kv_push(llnodep, leaf_node, tmp_nodes[j]);
        free(tmp_nodes);

        if (exceed_limit || leaf_node.n > PATH_N_LIM) {
            exceed_limit = 1;
            break;
        }
    }
    free(leaves);

    if (exceed_limit)
        goto final_clean;

#ifdef DEBUG_PATH_FINDER
    fprintf(stderr, "[DEBUG_PATH_FINDER::%s] number leaf nodes: %lu\n", __func__, leaf_node.n);
#endif

    for (i = 0; i < leaf_node.n; ++i) {
        kvec_t(uint32_t) path;
        kv_init(path);
        node = leaf_node.a[i];
        path.n = 0;
        while (node != NULL) {
            kv_push(uint32_t, path, node->v);
            node = node->prev;
        }
        rev_array(path.a, path.n);

        circ = asmg_arc1(g, path.a[path.n-1], path.a[0]) != 0;

#ifdef DEBUG_PATH_FINDER
        uint64_t n = path.n;
        fprintf(stderr, "[DEBUG_PATH_FINDER::%s] Path %lu [%s] (%lu): %s%c", __func__, i, 
                circ? "circle" : "linear", n, asg->seg[path.a[0]>>1].name, "+-"[path.a[0]&1]);
        for (j = 1; j < n; ++j)
            fprintf(stderr, ",%s%c", asg->seg[path.a[j]>>1].name, "+-"[path.a[j]&1]);
        fprintf(stderr, "\n");
#endif

        // calculate sequence length
        asmg_arc_t *a;
        uint64_t l, l1, *l_seg;
        uint32_t *l_beg, *l_end;
        double cov, wl;
        
        MYMALLOC(l_seg, path.n);
        l_seg[0] = g->vtx[path.a[0]>>1].len;
        l = l_seg[0];
        cov = g->vtx[path.a[0]>>1].cov;
        wl = cov * l;
        for (j = 1; j < path.n; ++j) {
            a = asmg_arc1(g, path.a[j-1], path.a[j]);
            assert(!!a);
            l_seg[j-1] = l_seg[j-1] << 32 | a->ls;
            l_seg[j] = g->vtx[path.a[j]>>1].len;
            l1 = l_seg[j] - a->ls;
            cov = g->vtx[path.a[j]>>1].cov;
            l += l1;
            wl += cov * l1;
        }
        l_seg[path.n-1] <<= 32;
        
        l_beg = l_end = 0;
        if (circ) {
            a = asmg_arc1(g, path.a[path.n-1], path.a[0]);
            assert(!!a);
            l1 = a->ls;
            cov = g->vtx[path.a[0]>>1].cov;
            l -= l1;
            wl -= cov * l1;
        } else {
            MYMALLOC(l_beg, path.n);
            MYMALLOC(l_end, path.n);
            l_beg[0] = 0;
            for (j = 1; j < path.n; ++j)
                l_beg[j] = l_beg[j-1] + (l_seg[j-1] >> 32) - (uint32_t) l_seg[j-1];
            for (j = 0; j < path.n; ++j) 
                l_end[j] = l - l_beg[j] - (l_seg[j] >> 32);
            assert(l_end[path.n-1] == 0);
        }
        
        // replace copied vertices with originals
        khint32_t k32;
        for (j = 0; j < path.n; ++j) {
            k32 = kh_u32_get(seg_dups, path.a[j]>>1);
            if (k32 < kh_end(seg_dups))
                path.a[j] = kh_val(seg_dups, k32)<<1 | (path.a[j]&1);
        }

        path_t p = {0, path.n, circ, 0, path.a, l, wl, .0};
        kv_push(path_t, *paths, p);
        
        if (!circ) {
            // for linear paths add longest circular subpaths if at least sub_circ_minf sequences covered
            double drop, max_drop, min_drop;
            int beg, end, beg1, end1, n;
            
            uint64_t L = is_pltd? MIN(l, COMMON_AVG_PLTD_SIZE) : l;
            max_drop = l - (double) L * sub_circ_minf;
            n = path.n;
            beg1 = end1 = -1;
            min_drop = FLT_MAX;
            for (beg = 0; beg < n; ++beg) {
                if (l_beg[beg] > max_drop || l_beg[beg] >= min_drop)
                    break;
                for (end = n - 1; end >= beg; --end) {
                    drop = l_beg[beg] + l_end[end];
                    if (drop > max_drop || drop >= min_drop)
                        break;
                    // check if the path [beg, end] is a circular
                    a = asmg_arc1(g, path.a[end], path.a[beg]);
                    if (a) {
                        beg1 = beg;
                        end1 = end;
                        min_drop = drop;
                        break; // only need the longest subpath
                    }
                }
            }

            if (beg1 >= 0) {
                // circular subpath found
                // add subpath
                uint32_t *v;
                n = end1 - beg1 + 1;
                MYMALLOC(v, n);
                memcpy(v, &path.a[beg1], sizeof(uint32_t) * n);
                wl = (l_seg[beg1] >> 32) * g->vtx[path.a[beg1]>>1].cov;
                for (beg = beg1 + 1; beg <= end1; ++beg)
                    wl += ((l_seg[beg] >> 32) - (uint32_t) l_seg[beg-1]) * g->vtx[path.a[beg]>>1].cov;
                l -= l_beg[beg1] + l_end[end1];
                a = asmg_arc1(g, v[n-1], v[0]);
                l -= a->ls;
                wl -= a->ls * g->vtx[v[0]>>1].cov;
                path_t p1 = {0, n, 1, 0, v, l, wl, .0};
                kv_push(path_t, *paths, p1);
            }
        }
        
        free(l_seg);
        free(l_beg);
        free(l_end);
    }

final_clean:
    for (i = 0; i < root_node.n; ++i)
        llnode_destroy(root_node.a[i]);

    kv_destroy(root_node);
    kv_destroy(leaf_node);

    return;
}

static int pcmpfunc(const void *a, const void *b)
{
    path_t x, y;
    x = *(path_t *) a;
    y = *(path_t *) b;
    if (x.wlen != y.wlen)
        return (x.wlen < y.wlen) - (x.wlen > y.wlen);
    if (x.len != y.len)
        return (x.len < y.len) - (x.len > y.len);
    if (x.circ != y.circ)
        return x.circ? -1 : 1;
    if (x.srcc != y.srcc)
        return (x.srcc < y.srcc) - (x.srcc > y.srcc);
    if (x.nv != y.nv)
        return (x.nv < y.nv) - (x.nv > y.nv);

    return 0;
}

static char *strdup1(char *src, uint32_t n)
{
    char *dst;
    MYMALLOC(dst, n+1);
    memcpy(dst, src, n);
    dst[n] = '\0';
    return dst;
}

path_t make_path_from_str(asg_t *asg, char *path_str, char *sid)
{
    int circ;
    uint32_t v, cov;
    uint64_t len, len1;
    double wlen;
    kvec_t(uint32_t) vt;
    char *ptr, *s;

    kv_init(vt);

    while (isspace(*path_str) && *path_str != '\0')
        ++path_str;

    while (*path_str != '\0') {
        ptr = path_str;
        while (!isspace(*ptr) && *ptr != ',' && *ptr !='\0')
            ++ptr;

        if (*(ptr-1) == '+' || *(ptr-1) == '-') {
            s = strdup1(path_str, ptr - path_str - 1);
            v = asg_name2id(asg, s);
            if (v == UINT32_MAX) {
                fprintf(stderr, "[E::%s] sequence does not exist: %s\n", __func__, s);
                exit(EXIT_FAILURE);
            }
            v = (v<<1) | (*(ptr-1)=='-');
            kv_push(uint32_t, vt , v);
            free(s);
        } else {
            fprintf(stderr, "[E::%s] invalid path string: %s\n", __func__, path_str);
            exit(EXIT_FAILURE);
        }

        if (isspace(*ptr) || *ptr == '\0') break;
        path_str = ptr + 1;
    }
    
    if (vt.n == 0) {
        fprintf(stderr, "[E::%s] invalid path string: %s\n", __func__, path_str);
        exit(EXIT_FAILURE);
    }

    size_t i;
    asmg_t *g;
    asmg_arc_t *a;
    g = asg->asmg;
    a = asmg_arc1(g, vt.a[vt.n-1], vt.a[0]);
    circ = !!a;
    len = g->vtx[vt.a[0]>>1].len;
    cov = g->vtx[vt.a[0]>>1].cov;
    wlen = (double) cov * len;
    if (circ) len -= a->ls, wlen -= cov * a->ls;
    for (i = 1; i < vt.n; ++i) {
        len1 = g->vtx[vt.a[i]>>1].len;
        cov = g->vtx[vt.a[i]>>1].cov;
        len += len1;
        wlen += (double) cov * len1;
        a = asmg_arc1(g, vt.a[i-1], vt.a[i]);
        if (!a) {
            fprintf(stderr, "[W::%s] gap introduced as link does not exist: %s%c -> %s%c\n", __func__,
                    asg->seg[vt.a[i-1]>>1].name, "+-"[vt.a[i-1]&1],
                    asg->seg[vt.a[i]>>1].name, "+-"[vt.a[i]&1]);
        } else {
            len -= a->ls;
            wlen -= (double) cov * a->ls;
        }
    }

    path_t path = {sid? strdup1(sid, strlen(sid)) : 0, vt.n, circ, 0, vt.a, len, wlen, .0};
    
    return path;
}

void path_sort(path_v *paths)
{
    // sort by wlen -> len -> circ -> nv
    qsort(paths->a, paths->n, sizeof(path_t), pcmpfunc);

    size_t i;
    double b_ll, b_cl;
    b_ll = b_cl = .0;
    for (i = 0 ; i < paths->n; ++i) {
        if (!paths->a[i].circ && paths->a[i].wlen > b_ll)
            b_ll = paths->a[i].wlen;
        if (paths->a[i].circ && paths->a[i].wlen > b_cl)
            b_cl = paths->a[i].wlen;
    }
    if (b_cl >= b_ll)
        b_ll = DBL_MAX;
    // find best
    for (i = 0 ; i < paths->n; ++i) {
        if (!paths->a[i].circ && paths->a[i].wlen >= b_ll)
            paths->a[i].best = 1;
        if (paths->a[i].circ && paths->a[i].wlen >= b_cl)
            paths->a[i].best = 1;
    }
}

static void swap(uint32_t arr[], uint32_t f, uint32_t s, uint32_t d)
{
    // swaps d elements starting at index f with d elements starting at index s
    uint32_t i, temp;
    for(i = 0; i < d; i++) {
        temp = arr[f + i];
        arr[f + i] = arr[s + i];
        arr[s + i] = temp;
    }
}

static void array_left_rotate(uint32_t arr[], uint32_t d, uint32_t n)
{
    uint32_t i, j;
    // if number of elements to be rotated is more than array size
    if(d > n)
        d = d % n;
    if(d == 0 || d == n)
        return;
    i = d;
    j = n - d;
    while (i != j) {
        if(i < j) { // left block arr[0..d-1] is shorter
            swap(arr, d-i, d+j-i, i);
            j -= i;
        } else { // rigth block arr[d..n-1] is shorter
            swap(arr, d-i, d, j);
            i -= j;
        }
    }
    // finally, block swap A and B
    swap(arr, d-i, d, i);
}

KHASHL_MAP_INIT(KH_LOCAL, kh_str_t, kh_str, kh_cstr_t, uint64_t, kh_hash_str, kh_eq_str)
KHASHL_MAP_INIT(KH_LOCAL, kh_s32_t, kh_s32, khint32_t, uint32_t, kh_hash_dummy, kh_eq_generic);

static double path_rotate_core(asg_t *g, path_t *path, hmm_annot_db_t *annots, OG_TYPE_t og_type)
{
    uint32_t i, j;
    const uint32_t g_n = athaliana_pltd_g71n;
    const char *const *genes = athaliana_pltd_g71;
    kh_str_t *gene_db;
    kh_s32_t *segs;
    khint32_t k, k1;
    uint64_t v;
    hmm_annot_t *annot, *annot1;
    double coeff = .0;
    int absent;

    assert(g_n < 256);

    // add core genes to table
    gene_db = kh_str_init();
    for (i = 0; i < g_n; ++i) {
        k = kh_str_put(gene_db, genes[i], &absent);
        kh_val(gene_db, k) = (uint64_t) i << 32 | UINT32_MAX;
    }

    // only need segs in the path
    segs = kh_s32_init();
    for (i = 0; i < path->nv; ++i) {
        k = kh_s32_put(segs, path->v[i]>>1, &absent);
        if (absent)
            kh_val(segs, k) = 1;
        else
            kh_val(segs, k) += 1;
    }

    // find (best) positions for core genes
    // on segs included in the path
    for (i = 0; i < annots->n; ++i) {
        annot = &annots->a[i];
        if (annot->og_type == (uint32_t) og_type) {
            k = kh_str_get(gene_db, annot->gname);
            k1 = kh_s32_get(segs, asg_name2id(g, annot->sname)); // sid might not set
            if (k < kh_end(gene_db) && k1 < kh_end(segs) && kh_val(segs, k1) == 1) {
                v = kh_val(gene_db, k);
                annot1 = (uint32_t) v == UINT32_MAX? 0 : &annots->a[(uint32_t) v];
                if (annot1 == 0 || annot1->score < annot->score)
                    kh_val(gene_db, k) = v >> 32 << 32 | i;
            }
        }
    }

    // do rotation if path is circular
    if (path->circ) {
        uint32_t s;
        for (i = 0, s = UINT32_MAX; i < g_n; ++i) {
            k = kh_str_get(gene_db, genes[i]);
            v = kh_val(gene_db, k);
            if ((uint32_t) v != UINT32_MAX) {
                s = (uint32_t) v;
                break;
            }
        }

        if (s != UINT32_MAX) {
            s = asg_name2id(g, annots->a[s].sname);
            // find start seg
            uint32_t t;
            for (i = 0, t = UINT32_MAX; i < path->nv; ++i) {
                if (path->v[i] >> 1 == s) {
                    t = i;
                    break;
                }
            }
            assert(t != UINT32_MAX);
            // do rotation to make t the first seg
            array_left_rotate(path->v, t, path->nv);
        }
    }

    // make gene orders in path
    kvec_t(uint64_t) g_ord;
    uint64_t w;
    kv_init(g_ord);
    for (k = (khint32_t) 0; k < kh_end(gene_db); ++k) {
        if (kh_exist(gene_db, k) && (uint32_t) (v = kh_val(gene_db, k)) != UINT32_MAX) {
            annot = &annots->a[(uint32_t) v];
            w  = (uint64_t) asg_name2id(g, annot->sname) << 40;
            w |= (uint64_t) (annot->alifrom + annot->alito) >> 1 << 8;
            w |= v >> 32;
            kv_push(uint64_t, g_ord, w);
        }
    }

    if (g_ord.n == 0) goto rotate_done;

    qsort(g_ord.a, g_ord.n, sizeof(uint64_t), u64_cmpfunc);

    uint64_t *idx;
    uint32_t last;
    MYCALLOC(idx, g->n_seg);
    for (i = 1, last = 0; i <= g_ord.n; ++i) {
        if (i == g_ord.n || (g_ord.a[i-1] >> 40) != (g_ord.a[i] >> 40))
            idx[g_ord.a[i-1] >> 40] = (uint64_t) last << 32 | (i - last), last = i;
    }

    uint32_t *p_ord, s, p, n, m;
    MYMALLOC(p_ord, g_ord.n);
    m = 0;
    for (i = 0; i < path->nv; ++i) {
        s = path->v[i] >> 1;
        p = idx[s] >> 32;
        n = (uint32_t) idx[s];
        if (n == 0) continue;
        if (path->v[i] & 1) {
            for (j = 0; j < n; ++j)
                p_ord[m++] = g_ord.a[p+n-1-j] & 0xFF;
        } else {
            for (j = 0; j < n; ++j)
                p_ord[m++] = g_ord.a[p+j] & 0xFF;
        }
    }

    assert(m == g_ord.n);

    // wrap gaps
    uint32_t *p_gap;
    MYCALLOC(p_gap, g_n);
    for (i = 0; i < m; ++i)
        ++p_gap[p_ord[i]];
    for (i = 1; i < g_n; ++i)
        p_gap[i] += p_gap[i-1];
    for (i = 0; i < m; ++i)
        p_ord[i] -= p_ord[i] - p_gap[p_ord[i]] + 1;
    free(p_gap);

    // calculate Spearman's rank correlation coefficient
    double ds, srcc;
    ds = .0;
    for (i = 0; i < g_ord.n; ++i)
        ds += ((double) p_ord[i] - i) * ((double) p_ord[i] - i);
    srcc = 1. - 6 * ds / g_ord.n / ((double) g_ord.n * g_ord.n - 1);

#ifdef DEBUG_SRCC
    // calculate Kendall's tau rank correlation coefficient
    double tau;
    ds = .0;
    for (i = 0; i < g_ord.n; ++i)
        for (j = i + 1; j < g_ord.n; ++j)
            ds += p_ord[i] > p_ord[j];
    tau = 1. - 2 * ds / g_ord.n / (g_ord.n - 1);
    fprintf(stderr, "[DEBUG_SRCC::%s] SRCC, %.6f; TAU, %.6f; N, %d\n", __func__, srcc, tau, m);
    for (i = 0; i < m; ++i)
        fprintf(stderr, "[DEBUG_SRCC::%s] %d %d\n", __func__, i, p_ord[i]);
#endif

    coeff = srcc;

    free(idx);
    free(p_ord);

rotate_done:
    kv_destroy(g_ord);
    kh_str_destroy(gene_db);
    kh_s32_destroy(segs);

    return coeff;
}

void path_rotate(asg_t *g, path_t *path, hmm_annot_db_t *annots, OG_TYPE_t og_type)
{
    double coeff, coeff_rev;
    coeff = path_rotate_core(g, path, annots, og_type);

    rev_path(path);
    coeff_rev = path_rotate_core(g, path, annots, og_type);

    if (coeff > coeff_rev) {
        // redo rev complement
        rev_path(path);
        if (path->circ)
            // make the last element the first
            array_left_rotate(path->v, path->nv-1, path->nv);
    } else {
        coeff = coeff_rev;
    }

    path->srcc = coeff;
}

void path_stats(asg_t *asg, path_v *paths, FILE *fo)
{
    uint32_t i, j, n;
    int nd_i, nd_n, nd_l, nd_w;
    uint64_t max_l;
    double max_wl;
    path_t path;

    nd_n = 0;
    max_l = 0;
    max_wl = .0;
    for (i = 0; i < paths->n; ++i) {
        if (paths->a[i].nv > nd_n)
            nd_n = paths->a[i].nv;
        if (paths->a[i].len > max_l)
            max_l = paths->a[i].len;
        if (paths->a[i].wlen > max_wl)
            max_wl = paths->a[i].wlen;
    }
    nd_i = floor(log10(fabs(paths->n))) + 1;
    nd_n = floor(log10(fabs(nd_n))) + 1;
    nd_l = floor(log10(fabs(max_l))) + 1;
    nd_w = floor(log10(fabs(max_wl))) + 3;

    for (i = 0; i < paths->n; ++i) {
        path = paths->a[i];
        n = path.nv;
        fprintf(fo, "%c %-*u %s %-*u %-*u %-*.1f %.3f %s%c", path.best? '*' : '#', nd_i, i, path.circ? "circle" : "linear",
                nd_n, n, nd_l, path.len, nd_w, path.wlen, path.srcc, asg->seg[path.v[0]>>1].name, "+-"[path.v[0]&1]);
        for (j = 1; j < n; ++j)
            fprintf(fo, ",%s%c", asg->seg[path.v[j]>>1].name, "+-"[path.v[j]&1]);
        fprintf(fo, "\n");
    }
}

static char comp_table[] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
    32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
    48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', 123, 124, 125, 126, 127
};

static void put_chars(char *seq, int len, int rv, int ow, FILE *fo, uint64_t *l, int line_wd)
{
    int i;
    if (!rv) {
        // forward
        for (i = ow; i < len; ++i) {
            fputc(seq[i], fo);
            if (++(*l) % line_wd == 0)
                fputc('\n', fo);
        }
    } else {
         // reverse
        for (i = len - ow - 1; i >= 0; --i) {
            fputc(comp_table[(int) seq[i]], fo);
            if (++(*l) % line_wd == 0)
                fputc('\n', fo);
        }
    }
}

static void make_gap(FILE *fo, uint64_t *l, int line_wd, int gap_size)
{
    int i;
    for (i = 0; i < gap_size; ++i) {
        fputc('N', fo);
        if (++(*l) % line_wd == 0)
            fputc('\n', fo);
    }
}

void print_seq(asg_t *asg, path_t *path, FILE *fo, int id, int force_linear, int line_wd, int gap_size)
{
    uint32_t i, n, v, lo, cov, n_gap;
    uint64_t l;
    asmg_t *g;
    asmg_arc_t *a;
    
    n = path->nv;
    if (n == 0) return;

    for (i = 0; i < n; ++i) {
        v = path->v[i];
        if (!asg->seg[v>>1].seq) {
            fprintf(stderr, "[E::%s] cannot make FASTA output: sequence not included in the GFA file\n", __func__);
            return;
        }
    }

    g = asg->asmg;
    lo = 0;
    cov = 0;
    if (path->circ && force_linear) {
        a = asmg_arc1(g, path->v[n-1], path->v[0]);
        assert(!!a);
        lo = a->ls;
        cov = g->vtx[path->v[0]>>1].cov;
    }
    
    if (path->sid)
        fprintf(fo, ">%s\tlength=%u wlength=%.1f nv=%u circular=%s path=%s%c", path->sid, path->len + lo, path->wlen + (double) cov * lo,
                path->nv, (force_linear || !path->circ)? "false" : "true", asg->seg[path->v[0]>>1].name, "+-"[path->v[0]&1]);
    else
        fprintf(fo, ">ctg%06d%c\tlength=%u wlength=%.1f nv=%u circular=%s path=%s%c", id, (force_linear || !path->circ)? 'l' : 'c', 
                path->len + lo, path->wlen + (double) cov * lo, path->nv, (force_linear || !path->circ)? "false" : "true", 
                asg->seg[path->v[0]>>1].name, "+-"[path->v[0]&1]);

    for (i = 1; i < n; ++i)
        fprintf(fo, ",%s%c", asg->seg[path->v[i]>>1].name, "+-"[path->v[i]&1]);
    fprintf(fo, "\n");

    l = 0;
    v = path->v[0];
    lo = (force_linear || !path->circ)? 0 : asmg_arc1(g, path->v[n-1], v)->ls;
    put_chars(asg->seg[v>>1].seq, asg->seg[v>>1].len, v&1, lo, fo, &l, line_wd);

    n_gap = 0;
    for (i = 1; i < n; ++i) {
        v = path->v[i];
        a = asmg_arc1(g, path->v[i-1], v);
        if (!!a) {
            put_chars(asg->seg[v>>1].seq, asg->seg[v>>1].len, v&1, a->ls, fo, &l, line_wd);
        } else {
            make_gap(fo, &l, line_wd, gap_size);
            put_chars(asg->seg[v>>1].seq, asg->seg[v>>1].len, v&1, 0, fo, &l, line_wd);
            ++n_gap;
        }
    }

    if (!path->circ || !force_linear)
        assert(l - (uint64_t) n_gap * gap_size == path->len);
    if (l % line_wd != 0)
        fputc('\n', fo);
}

void path_add_hmm_annot_bed6(hmm_annot_bed6_db_t *bed_annots, hmm_annot_db_t *annot_db, asg_t *asg, path_t *path, 
        int id, int force_linear, int gap_size, OG_TYPE_t og_type, double max_evalue)
{
    uint32_t i, n, v, lo, n_gap;
    uint64_t l;
    asmg_t *g;
    asmg_arc_t *a;
    char *cname;
    
    n = path->nv;
    if (n == 0) return;

    cname = 0;
    if (path->sid) {
        cname = strdup(path->sid);
    } else {
        MYMALLOC(cname, 12);
        sprintf(cname, "ctg%06d%c", id, (force_linear || !path->circ)? 'l' : 'c');
    }
    
    if (bed_annots->n_seg == bed_annots->m_seg) {
        ++bed_annots->m_seg;
        kroundup32(bed_annots->m_seg);
        MYREALLOC(bed_annots->snames, bed_annots->m_seg);
    }
    bed_annots->snames[bed_annots->n_seg++] = cname;

    g = asg->asmg;
    l = 0;
    v = path->v[0];
    lo = (force_linear || !path->circ)? 0 : asmg_arc1(g, path->v[n-1], v)->ls;
    hmm_annot_bed6_sname_add(bed_annots, annot_db, cname, asg->seg[v>>1].name, asg->seg[v>>1].len, lo, v&1, l, og_type, max_evalue);
    l += asg->seg[v>>1].len - lo;

    n_gap = 0;
    for (i = 1; i < n; ++i) {
        v = path->v[i];
        a = asmg_arc1(g, path->v[i-1], v);
        if (!a) l += gap_size, ++n_gap;
        hmm_annot_bed6_sname_add(bed_annots, annot_db, cname, asg->seg[v>>1].name, asg->seg[v>>1].len, a->ls, v&1, l, og_type, max_evalue);
        l += asg->seg[v>>1].len - a->ls;
    }

    if (!path->circ || !force_linear)
        assert(l - (uint64_t) n_gap * gap_size == path->len);
}

uint32_t select_best_seq(asg_t *g, path_v *paths, FILE *fo, int type, double seq_cf, int seq_id, int is_pltd)
{
    if (paths->n == 0)
        return UINT32_MAX;

    uint64_t l;
    size_t i, j, k;
    l = 0, j = 0;
    for (i = 0; i < paths->n; ++i) {
        if ((paths->a[i].circ || !type) && paths->a[i].len > l) {
            l = paths->a[i].len;
            j = i;
        }
    }
    
    if (!paths->a[j].circ) {
        k = UINT32_MAX, l = 0;
        for (i = 0; i < paths->n; ++i) {
            if (paths->a[i].circ && paths->a[i].len > l) {
                l = paths->a[i].len;
                k = i;
            }
        }
        if (k != UINT32_MAX) {
            uint64_t L = paths->a[j].len;
            if (is_pltd) L = MIN(L, COMMON_AVG_PLTD_SIZE);
            if ((double) l / L >= seq_cf) j = k;
        }
    }

    // TODO need more sophisticated rules for cmp_coeff
    if (is_pltd) {
        double coeff;
        uint32_t circ;
        circ = (paths->a[j].circ || type);
        k = UINT32_MAX, coeff = .0;
        for (i = 0; i < paths->n; ++i) {
            if ((paths->a[i].circ || !circ) && paths->a[i].srcc > coeff) {
                coeff = paths->a[i].srcc;
                k = i;
            }
        }
        if (k != UINT32_MAX && paths->a[k].len + 1000 >= paths->a[j].len)
            j = k;
    }
    
    if (fo) print_seq(g, &paths->a[j], fo, seq_id > 0? seq_id : 1, 0, 60, 100);

    return j;
}

void print_all_best_seqs(asg_t *g, path_v *paths, FILE *fo)
{
    if (paths->n == 0)
        return;

    size_t i;
    for (i = 0; i < paths->n; ++i)
        if (!!paths->a[i].best)
            print_seq(g, &paths->a[i], fo, i, 0, 60, 100);
}



/******************
 * Assembly Graph *
 * ****************/

KSTREAM_INIT(gzFile, gzread, 65536)
KHASHL_MAP_INIT(KH_LOCAL, kh_sdict_t, kh_sdict, kh_cstr_t, uint32_t, kh_hash_str, kh_eq_str)
typedef kh_sdict_t sdhash_t;

asg_t *asg_init()
{
    asg_t *g;
    MYCALLOC(g, 1);
    g->h_seg = kh_sdict_init();
    MYCALLOC(g->asmg, 1);
    return g;
}

static void asg_seg_destroy(asg_seg_t *seg)
{
    if (seg->name) free(seg->name);
    if (seg->seq) free(seg->seq);
}

void asg_destroy(asg_t *g)
{
    if (!g) return;
    uint64_t i;
    for (i = 0; i < g->n_seg; ++i)
        asg_seg_destroy(&g->seg[i]);
    free(g->seg);
    // key has been freed by asg_seg_destroy
    if (g->h_seg) kh_sdict_destroy(g->h_seg);
    if (g->asmg) asmg_destroy(g->asmg);
    free(g);
}

uint32_t asg_add_seg(asg_t *g, char *name, int allow_dups)
{
    if (!name)
        return UINT32_MAX;
    sdhash_t *h = g->h_seg;
    khint_t k;
    int absent;
    k = kh_sdict_put(h, name, &absent);
    if (absent) {
        asg_seg_t *s;
        if (g->n_seg == g->m_seg)
            MYEXPAND(g->seg, g->m_seg);
        s = &g->seg[g->n_seg];
        s->len = 0;
        s->seq = 0;
        s->cov = 0;
        kh_key(h, k) = s->name = strdup(name);
        kh_val(h, k) = g->n_seg++;
    } else if (!allow_dups) {
        fprintf(stderr, "[E::%s] duplicate segment '%s'\n", __func__, name);
        exit(EXIT_FAILURE);
    }
    return kh_val(h, k);
}

uint32_t asg_add_seg1(asg_t *g, char *name, char *seq, uint32_t len, uint64_t cov, int allow_dups)
{
    uint32_t k = asg_add_seg(g, name, allow_dups);
    asg_seg_t *s = &g->seg[k];
    s->seq = strdup(seq);
    s->len = len;
    s->cov = cov;
    return k;
}

uint32_t asg_name2id(asg_t *g, char *name)
{
    sdhash_t *h = g->h_seg;
    khint_t k;
    k = kh_sdict_get(h, name);
    return k == kh_end(h)? UINT32_MAX : kh_val(h, k);
}

static void asg_update_seg_seq(asg_seg_t *seg, uint32_t l, char *s)
{
    if (seg->seq) free(seg->seq);
    MYMALLOC(seg->seq, l+1);
    memcpy(seg->seq, s, l);
    seg->seq[l] = 0;
    seg->len = l;
    seg->cov = 0;
}

asmg_t *asg_make_asmg_copy(asmg_t *g, asmg_t *_g)
{
    asmg_t *g1;
    
    if (_g)
        g1 = _g;
    else
        MYCALLOC(g1, 1);
    
    g1->n_vtx = g1->m_vtx = g->n_vtx;
    MYCALLOC(g1->vtx, g1->n_vtx);
    // this works as a and seq are always 0
    memcpy(g1->vtx, g->vtx, sizeof(asmg_vtx_t) * g1->n_vtx);
    g1->n_arc = g1->m_arc = g->n_arc;
    MYCALLOC(g1->arc, g1->n_arc);
    memcpy(g1->arc, g->arc, sizeof(asmg_arc_t) * g1->n_arc);

    uint64_t n_vtx = asmg_vtx_n(g1);
    MYCALLOC(g1->idx_p, n_vtx);
    memcpy(g1->idx_p, g->idx_p, sizeof(uint64_t) * n_vtx);
    MYCALLOC(g1->idx_n, n_vtx);
    memcpy(g1->idx_n, g->idx_n, sizeof(uint64_t) * n_vtx);
    
    return g1;
}

asg_t *asg_make_copy(asg_t *asg)
{
    uint64_t i;
    khint_t k;
    int absent;

    // seg sequence not copied
    asg_t *asg1 = asg_init();
    asg_seg_t *seg, *seg1;
    asg1->m_seg = asg1->n_seg = asg->n_seg;
    MYCALLOC(asg1->seg, asg1->n_seg);
    kh_sdict_t *h_seg = (kh_sdict_t *) asg1->h_seg;
    for (i = 0; i < asg1->n_seg; ++i) {
        seg = &asg->seg[i];
        seg1 = &asg1->seg[i];
        seg1->name = strdup(seg->name);
        seg1->len = seg->len;
        seg1->cov = seg->cov;
        k = kh_sdict_put(h_seg, seg1->name, &absent);
        kh_val(h_seg, k) = i;
    }
    asg_make_asmg_copy(asg->asmg, asg1->asmg);

    return asg1;
}

int clean_graph_by_sequence_coverage(asg_t *asg, double min_cf, int max_copy, int verbose)
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
        if (visited[i] || g->vtx[i].del) continue;
        asg->asmg = asg_make_asmg_copy(g, 0);
        asmg_subgraph(asg->asmg, &i, 1, 0, 0, 0, 1);
        avg_cov = graph_sequence_coverage_precise(asg, min_cf, 0, max_copy, 0);
        if (verbose > 1)
            fprintf(stderr, "[M::%s] subgraph seeding from %s per copy average coverage: %.3f\n", 
                    __func__, asg->seg[i].name, avg_cov);
        for (j = 0; j < n_seg; ++j) {
            if (asg->asmg->vtx[j].del) continue;
            if (avg_cov && asg->asmg->vtx[j].cov / avg_cov < min_cf)
                kv_push(uint32_t, rm_v, j);
            visited[j] = 1;
        }
        asmg_destroy(asg->asmg);
    }
    asg->asmg = g; // roll back

    nv = rm_v.n;
    for (i = 0; i < nv; ++i)
        asmg_vtx_del(g, rm_v.a[i], 1);
    asmg_finalize(g, 0);

    if (verbose > 1) {
        if (verbose > 2) {
            fprintf(stderr, "[M::%s] graph after cleaning\n", __func__);
            asg_print(asg, stderr, 1);
        }
        fprintf(stderr, "[M::%s] number sequence cleaned: %u\n", __func__, nv);
        for (i = 0; i < nv; ++i)
            fprintf(stderr, "[M::%s] %s removed\n", __func__, asg->seg[rm_v.a[i]].name);
        fprintf(stderr, "[M::%s] graph stats after cleaning\n", __func__);
        asg_stat(asg, stderr);
    }
    
    kv_destroy(rm_v);
    free(visited);

    return nv;
}

double sequence_covered_by_path(asg_t *asg, path_t *path, uint32_t len)
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

/****************
 * GFA File I/O *
 * **************/
// adapted from https://github.com/lh3/gfatools

#define aux_decimal_val(type, p) *((type *) (p))
// always as double
static double gfa_aux_decimal_value(uint8_t *aux)
{
    double val = .0;
    if (*aux == 'c') val = aux_decimal_val(int8_t, aux + 1);
    else if (*aux == 'C') val = aux_decimal_val(uint8_t, aux + 1);
    else if (*aux == 's') val = aux_decimal_val(int16_t, aux + 1);
    else if (*aux == 'S') val = aux_decimal_val(uint16_t, aux + 1);
    else if (*aux == 'i') val = aux_decimal_val(int64_t, aux + 1);
    else if (*aux == 'I') val = aux_decimal_val(uint64_t, aux + 1);
    else if (*aux == 'f') val = aux_decimal_val(double, aux + 1);
    return val;
}

static inline int gfa_aux_type2size(int x)
{
    if (x == 'C' || x == 'c' || x == 'A') return 1;
    else if (x == 'S' || x == 's') return 2;
    else if (x == 'I' || x == 'i' || x == 'f') return 8;
    else return 0;
}

#define __skip_tag(s) do { \
    int type = *(s); \
    ++(s); \
    if (type == 'Z') { while (*(s)) ++(s); ++(s); } \
    else if (type == 'B') (s) += 5 + gfa_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
    else (s) += gfa_aux_type2size(type); \
} while(0)

static uint8_t *gfa_aux_get(int l_data, const uint8_t *data, const char tag[2])
{
    const uint8_t *s = data;
    int y = tag[0]<<8 | tag[1];
    while (s < data + l_data) {
        int x = (int)s[0]<<8 | s[1];
        s += 2;
        if (x == y) return (uint8_t*)s;
        __skip_tag(s);
    }
    return 0;
}

// GFA tag
char TAG_ARC_COV[4]; // arc coverage EC:i
char TAG_SEQ_COV[4]; // seq coverage SC:f
char TAG_SBP_COV[4]; // seq total base coverage KC:i FC:i

int is_valid_gfa_tag(const char *tag)
{
    int is_valid = 0;
    if (strlen(tag) != 4)
        return is_valid;
    is_valid = isalpha(tag[0]) && 
        (isalpha(tag[1]) || isdigit(tag[1])) && 
        (tag[2] == ':') &&
        (tag[3] == 'A' || tag[3] == 'i' || tag[3] == 'f' || tag[3] == 'Z' || tag[3] == 'B');
    return is_valid;
}

static int gfa_aux_parse(char *s, uint8_t **data, int *max)
{
    char *q, *p;
    kstring_t str;
    if (s == 0) return 0;
    str.l = 0, str.m = *max, str.s = (char*)*data;
    if (*s == '\t') ++s;
    for (p = q = s;; ++p) {
        if (*p == 0 || *p == '\t') {
            int c = *p;
            *p = 0;
            if (p - q >= 5 && q[2] == ':' && q[4] == ':' && (q[3] == 'A' || q[3] == 'i' || q[3] == 'f' || q[3] == 'Z' || q[3] == 'B')) {
                int type = q[3];
                kputsn_(q, 2, &str);
                q += 5;
                if (type == 'A') {
                    kputc_('A', &str);
                    kputc_(*q, &str);
                } else if (type == 'i') {
                    int64_t x;
                    x = strtoll(q, &q, 10);
                    kputc_(type, &str); kputsn_((char*)&x, 8, &str);
                } else if (type == 'f') {
                    double x;
                    x = strtod(q, &q);
                    kputc_('f', &str); kputsn_(&x, 8, &str);
                } else if (type == 'Z') {
                    kputc_('Z', &str); kputsn_(q, p - q + 1, &str); // note that this include the trailing NULL
                } else if (type == 'B') {
                    type = *q++; // q points to the first ',' following the typing byte
                    if (p - q >= 2 && (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I' || type != 'f')) {
                        int32_t n;
                        char *r;
                        for (r = q, n = 0; *r; ++r)
                            if (*r == ',') ++n;
                        kputc_('B', &str); kputc_(type, &str); kputsn_(&n, 4, &str);
                        // TODO: to evaluate which is faster: a) aligned array and then memmove(); b) unaligned array; c) kputsn_()
                        if (type == 'c')      while (q + 1 < p) { int8_t   x = strtol(q + 1, &q, 0);  kputc_(x, &str); }
                        else if (type == 'C') while (q + 1 < p) { uint8_t  x = strtol(q + 1, &q, 0);  kputc_(x, &str); }
                        else if (type == 's') while (q + 1 < p) { int16_t  x = strtol(q + 1, &q, 0);  kputsn_(&x, 2, &str); }
                        else if (type == 'S') while (q + 1 < p) { uint16_t x = strtol(q + 1, &q, 0);  kputsn_(&x, 2, &str); }
                        else if (type == 'i') while (q + 1 < p) { int64_t  x = strtoll(q + 1, &q, 0); kputsn_(&x, 8, &str); }
                        else if (type == 'I') while (q + 1 < p) { uint64_t x = strtoll(q + 1, &q, 0); kputsn_(&x, 8, &str); }
                        else if (type == 'f') while (q + 1 < p) { double   x = strtod(q + 1, &q);     kputsn_(&x, 8, &str); }
                    }
                } // should not be here, as we have tested all types
            }
            q = p + 1;
            if (c == 0) break;
        }
    }
    if (str.l > 0 && str.l == str.m) ks_resize(&str, str.l + 1);
    if (str.s) str.s[str.l] = 0;
    *max = str.m, *data = (uint8_t*)str.s;
    return str.l;
}

#define PARSE_A_ERR -1
#define PARSE_Q_ERR -2
#define PARSE_S_ERR -3
#define PARSE_L_ERR -4

static int asg_parse_fa_hdr(asg_t *g, char *s, asg_seg_t **seg)
{
    uint64_t i;
    char *c = s;
    while (*c != 0 && !isspace(*c))
        ++c;
    *c = 0;
    if (c - s == 1) // empty
        return PARSE_A_ERR;
    i = asg_add_seg(g, s + 1, 0);
    *seg = &g->seg[i];
    return 0;
}

static inline int gfa_parse_S(asg_t *g, char *s)
{
    if (*s != 'S') return PARSE_S_ERR;

    int i, c, is_ok;
    char *p, *q, *seg, *seq, *rest;
    uint32_t sid, len;
    
    seg = seq = rest = 0;
    len = 0;
    for (i = 0, p = q = s + 2;; ++p) {
        if (*p == 0 || *p == '\t') {
            c = *p;
            *p = 0;
            if (i == 0) seg = q;
            else if (i == 1) {
                seq = q[0] == '*'? 0 : strdup(q);
                is_ok = 1, rest = c? p + 1 : 0;
                break;
            }
            ++i, q = p + 1;
            if (c == 0) break;
        }
    }

    if (is_ok) { // all mandatory fields read
        // parse sequence coverage if presented
        int l_aux, m_aux = 0;
        uint32_t LN = 0;
        uint8_t *aux = 0, *s_LN = 0;
        asg_seg_t *s;
        l_aux = gfa_aux_parse(rest, &aux, &m_aux); // parse optional tags
        s_LN = l_aux? gfa_aux_get(l_aux, aux, "LN") : 0;
        if (s_LN && s_LN[0] == 'i')
            LN = *(int64_t*)(s_LN + 1);
        if (seq == 0) {
            if (LN > 0) len = LN;
        } else {
            len = strlen(seq);
        }
        if (LN > 0 && len != LN)
            fprintf(stderr, "[W::%s] for segment '%s', LN:i:%u tag is different from sequence length %d\n", __func__, seg, LN, len);
        sid = asg_add_seg(g, seg, 0);
        s = &g->seg[sid];
        s->len = len, s->seq = seq;
        if (l_aux > 0) {
            uint8_t *s_SBP_COV = 0, *s_SEQ_COV = 0;
            double dv = 0;
            char tag[2];
            if (TAG_SBP_COV[0] != 0) {
                memcpy(tag, TAG_SBP_COV, 2);
                s_SBP_COV = gfa_aux_get(l_aux, aux, tag);
                if (s_SBP_COV && *s_SBP_COV == TAG_SBP_COV[3]) {
                    dv = gfa_aux_decimal_value(s_SBP_COV);
                    s->cov = len > 0? dv/len : dv;
                } else {
                    fprintf(stderr, "[W::%s] for segment '%s', %s tag is absent\n", __func__, seg, tag);
                }
            } else if (TAG_SEQ_COV[0] != 0) {
                memcpy(tag, TAG_SEQ_COV, 2);
                s_SEQ_COV = gfa_aux_get(l_aux, aux, tag);
                if (s_SEQ_COV && *s_SEQ_COV == TAG_SEQ_COV[3]) {
                    s->cov = gfa_aux_decimal_value(s_SEQ_COV);
                } else {
                    fprintf(stderr, "[W::%s] for segment '%s', %s tag is absent\n", __func__, seg, tag);
                }
            } else {
                // check KC and FC
                s_SBP_COV = gfa_aux_get(l_aux, aux, "KC");
                if (s_SBP_COV && *s_SBP_COV == 'i') {
                    dv = *(int64_t*)(s_SBP_COV + 1);
                } else {
                    s_SBP_COV = gfa_aux_get(l_aux, aux, "FC");
                    if (s_SBP_COV && *s_SBP_COV == 'i')
                        dv = *(int64_t*)(s_SBP_COV + 1);
                }
                s->cov = len > 0? dv/len : dv;
            }
        }
        if (s->cov == 0) {
            fprintf(stderr, "[W::%s] the coverage of segment '%s' is zero\n", __func__, seg);
            s->cov = 1;
        }
        free(aux);
    } else return PARSE_S_ERR;

    return 0;
}

static int gfa_parse_L(asg_t *g, char *s)
{
    if (*s != 'L') return PARSE_L_ERR;

    int i, is_ok = 0;
    char *p, *q, *segv, *segw, *rest;
    uint64_t ov, ow, oriv, oriw;
    
    segv = segw = rest = 0;
    ov = ow = oriv = oriw = UINT64_MAX;
    for (i = 0, p = q = s + 2;; ++p) {
        if (*p == 0 || *p == '\t') {
            int c = *p;
            *p = 0;
            if (i == 0) {
                segv = q;
            } else if (i == 1) {
                if (*q != '+' && *q != '-') return PARSE_L_ERR;
                oriv = (*q != '+');
            } else if (i == 2) {
                segw = q;
            } else if (i == 3) {
                if (*q != '+' && *q != '-') return PARSE_L_ERR;
                oriw = (*q != '+');
            } else if (i == 4) {
                if (*q == '*') {
                    ov = ow = 0;
                } else if (*q == ':') {
                    ov = UINT64_MAX;
                    ow = isdigit(*(q+1))? strtoul(q+1, &q, 10) : UINT64_MAX;
                } else if (isdigit(*q)) {
                    char *r;
                    ov = strtol(q, &r, 10);
                    if (isupper(*r)) { // CIGAR
                        ov = ow = 0;
                        do {
                            long l;
                            l = strtol(q, &q, 10);
                            if (*q == 'M' || *q == 'D' || *q == 'N') ov += l;
                            if (*q == 'M' || *q == 'I' || *q == 'S') ow += l;
                            ++q;
                        } while (isdigit(*q));
                    } else if (*r == ':') { // overlap lengths
                        ow = isdigit(*(r+1))? strtoul(r+1, &r, 10) : UINT64_MAX;
                    } else break;
                } else break;
                is_ok = 1, rest = c? p + 1 : 0;
                break;
            }
            ++i, q = p + 1;
            if (c == 0) break;
        }
    }
    if (i == 4 && is_ok == 0) ov = ow = 0, is_ok = 1; // no overlap field
    if (is_ok) {
        uint64_t v, w;
        int l_aux, m_aux = 0;
        uint8_t *aux = 0;
        asmg_arc_t *arc;
        v = asg_add_seg(g, segv, 1) << 1 | oriv;
        w = asg_add_seg(g, segw, 1) << 1 | oriw;
        arc = asmg_arc_add(g->asmg, v, w, 0, ov, UINT64_MAX, 0, 0);
        l_aux = gfa_aux_parse(rest, &aux, &m_aux); // parse optional tags
        if (l_aux) {
            uint8_t *s_ARC_COV = 0;
            char tag[2];
            if (TAG_ARC_COV[0] != 0) {
                memcpy(tag, TAG_ARC_COV, 2);
                s_ARC_COV = gfa_aux_get(l_aux, aux, tag);
                if (s_ARC_COV && *s_ARC_COV == TAG_ARC_COV[3])
                    arc->cov = gfa_aux_decimal_value(s_ARC_COV);
                else
                    fprintf(stderr, "[W::%s] for arc '%s%c' -> '%s%c', %s tag is absent\n", __func__, segv, "+-"[oriv], segw, "+-"[oriw], tag);
            } else {
                // check EC
                s_ARC_COV = gfa_aux_get(l_aux, aux, "EC");
                if (s_ARC_COV && *s_ARC_COV == 'i')
                    arc->cov = *(int64_t*)(s_ARC_COV + 1);
            }
        }
        if (arc->cov == 0) {
            fprintf(stderr, "[W::%s] the coverage of arc '%s%c' -> '%s%c' is zero\n", __func__, segv, "+-"[oriv], segw, "+-"[oriw]);
            arc->cov = 1;
        }
        free(aux);
    } else return PARSE_L_ERR;
    return 0;
}

static void asg_finalize_asmg(asg_t *g)
{
    if (g->asmg == 0)
        MYCALLOC(g->asmg, 1);
    asmg_t *asmg = g->asmg;
    uint64_t i;
    asmg_vtx_t *vtx;
    asmg->n_vtx = asmg->m_vtx = g->n_seg;
    MYCALLOC(asmg->vtx, g->n_seg);
    for (i = 0; i < g->n_seg; ++i) {
        vtx = &asmg->vtx[i];
        // will never need vtx->a
        // seg and vtx indices are consistent
        //vtx->n = 1;
        //MYMALLOC(vtx->a, 1);
        //vtx->a[0] = i;
        vtx->len = g->seg[i].len;
        vtx->cov = g->seg[i].cov;
    }
    asmg_finalize(asmg, 0);
}

asg_t *asg_read(const char *fn)
{
    int dret, ret, is_fa, is_fq, is_gfa;
    gzFile fp;
    kstream_t *ks;
    asg_seg_t *fa_seg;
    uint64_t lineno;
    asg_t *g;
    kstring_t s = {0, 0, 0}, fa_seq = {0, 0, 0};

    fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
    if (fp == 0) return 0;
    ks = ks_init(fp);
    g = asg_init();
    lineno = 0;
    fa_seg = 0;
    is_fa = is_fq = is_gfa = 0;
    
    while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
        ++lineno;
        if (s.l == 0) // empty line
            continue;

        ret = 0;
        if (!is_gfa && s.s[0] == '>') { // FASTA header
            is_fa = 1;
            if (fa_seg) asg_update_seg_seq(fa_seg, fa_seq.l, fa_seq.s);
            // parse header
            asg_parse_fa_hdr(g, s.s, &fa_seg);
            fa_seq.l = 0;
        } else if (!is_gfa && s.s[0] == '@') { // FASTQ header
            is_fq = 1;
            // parse and write header
            asg_parse_fa_hdr(g, s.s, &fa_seg);
            // parse sequence
            if(ks_getuntil(ks, KS_SEP_LINE, &s, &dret) < 0)
                ret = PARSE_Q_ERR;
            else
                asg_update_seg_seq(fa_seg, s.l, s.s);
            ++lineno;
            // skip quality score lines
            if (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) < 0 ||
                    ks_getuntil(ks, KS_SEP_LINE, &s, &dret) < 0)
                ret = PARSE_Q_ERR;
            lineno += 2;
        } else if (is_fa) { // FASTA sequence line
            // append sequence to seg
            kputsn(s.s, s.l, &fa_seq);
        } else {
            is_gfa = 1;
            if (s.s[0] == 'S')
                ret = gfa_parse_S(g, s.s);
            else if (s.s[0] == 'L')
                ret = gfa_parse_L(g, s.s);
        }

        if (ret < 0) {
            fprintf(stderr, "[E::%s] failed to parse %s file: %c-line at line %lu (error code %d)\n",
                    __func__, is_fa? "FASTA" : (is_fq? "FASTQ" : "GFA"), s.s[0], lineno, ret);
            exit(EXIT_FAILURE);
        }
    }
    if (is_fa && fa_seg) asg_update_seg_seq(fa_seg, fa_seq.l, fa_seq.s);

    free(fa_seq.s);
    free(s.s);
    ks_destroy(ks);
    gzclose(fp);

    // add vtx, fix symmetric arcs, sort and index etc
    asg_finalize_asmg(g);

    return g;
}

void asg_stat(asg_t *asg, FILE *fo)
{
    uint64_t i, nv, n_vtx, n_seg, max_deg, tot_seg_len, n_link, n_arc, tot_deg;
    asmg_t *g;

    n_vtx = n_seg = max_deg = tot_seg_len = n_link = n_arc = tot_deg = 0;
    g = asg->asmg;
    for (i = 0; i < asg->n_seg; ++i) {
        if (g->vtx[i].del) continue;
        tot_seg_len += asg->seg[i].len;
        ++n_seg;
    }
    fprintf(fo, "Number of segments: %lu\n", n_seg);
    fprintf(fo, "Total segment length: %lu\n", tot_seg_len);
    if (n_seg) fprintf(fo, "Average segment length: %.3f\n", (double) tot_seg_len / n_seg);

    for (i = 0; i < g->n_arc; ++i) {
        if (!g->arc[i].del) {
            ++n_arc;
            if (!g->arc[i].comp)
                ++n_link;
        }
    }
    fprintf(fo, "Number of links: %lu\n", n_link);
    fprintf(fo, "Number of arcs: %lu\n", n_arc);

    n_vtx = asmg_vtx_n(g);
    for (i = 0; i < n_vtx; ++i) {
        nv = asmg_arc_n1(g, i);
        if (nv > max_deg)
            max_deg = nv;
        tot_deg += nv;
    }
    fprintf(fo, "Max degree: %lu\n", max_deg);
    if (n_seg > 0) fprintf(fo, "Average degree: %.3f\n", (double) tot_deg / n_seg / 2);
}

void asg_print(asg_t *g, FILE *fo, int no_seq)
{
    uint32_t i, cov;
    uint64_t k;
    asmg_t *asmg = g->asmg;

    fprintf(fo, "H\tVN:Z:1.0\n");
    for (i = 0; i < g->n_seg; ++i) {
        const asg_seg_t *s = &g->seg[i];
        // seg and vtx indices are interchangeable
        if (asmg && asmg->vtx[i].del) continue;
        cov = asmg? asmg->vtx[i].cov : s->cov;
        fprintf(fo, "S\t%s\t", s->name);
        if (s->seq && !no_seq) fputs(s->seq, fo);
        else fputc('*', fo);
        fprintf(fo, "\tLN:i:%u\tKC:i:%lu\tSC:f:%.3f\n", s->len, (uint64_t)s->len*cov, (double)cov);
    }

    if (asmg == 0) return;
    for (k = 0; k < asmg->n_arc; ++k) {
        const asmg_arc_t *a = &asmg->arc[k];
        if (a->del || a->comp) continue;
        fprintf(fo, "L\t%s\t%c\t%s\t%c\t%luM\tEC:i:%u\n", g->seg[a->v>>1].name, "+-"[a->v&1], 
                g->seg[a->w>>1].name, "+-"[a->w&1], a->ls, a->cov);
    }
}

void asg_print_fa(asg_t *g, FILE *fo, int line_wd)
{
    uint64_t i, l;
    for (i = 0; i < g->n_seg; ++i) {
        if (g->asmg && g->asmg->vtx[i].del)
            continue;
        if (g->seg[i].seq == 0)
            fprintf(stderr, "[W::%s] skip empty sequence: %s\n", __func__, g->seg[i].name);
        fprintf(fo, ">%s\n", g->seg[i].name);
        l = 0;
        put_chars(g->seg[i].seq, g->seg[i].len, 0, 0, fo, &l, line_wd);
        if (l % line_wd != 0) fputc('\n', fo);
    }
}

/*****************************
 * Assembly Graph Annotation *
 * ***************************/

static int og_cmpfunc_r(const void *a, const void *b)
{
    return (((og_component_t *) a)->score < ((og_component_t *) b)->score) -
        (((og_component_t *) a)->score > ((og_component_t *) b)->score);
}

static int u64_cmpfunc_r(const void *a, const void *b)
{
    return (*(uint64_t *) a < *(uint64_t *) b) - (*(uint64_t *) a > *(uint64_t *) b);
}

static int dbl_cmpfunc_r(const void *a, const void *b)
{
    return (*(double *) a < *(double *) b) - (*(double *) a > *(double *) b);
}


void og_component_destroy(og_component_t *og_component)
{
    if (og_component->v) free(og_component->v);
    if (og_component->g) free(og_component->g);
    if (og_component->asmg) asmg_destroy(og_component->asmg);
}

void og_component_v_destroy(og_component_v *component_v)
{
    if (!component_v) return;
    size_t i;
    for (i = 0; i < component_v->n; ++i)
        og_component_destroy(&component_v->a[i]);
    if (component_v->a)
        free(component_v->a);
    free(component_v);
}

static uint32_t max2(double *a, uint32_t n, uint32_t *_smax)
{
    uint32_t i, imax, smax;
    double max_a, smax_a;
    imax = smax = 0;
    max_a = smax_a = -DBL_MAX;
    for (i = 0; i < n; ++i) {
        if (a[i] > max_a) {
            smax = imax;
            smax_a = max_a;
            imax = i;
            max_a = a[i];
        } else if (a[i] > smax_a) {
            smax = i;
            smax_a = a[i];
        }
    }
    if (_smax) *_smax = smax;
    return imax;
}

static void fix_og_misclassification(og_component_v *component_v, int verbose)
{
    // TODO fix this by maybe improving the specificity of the HMM database
    // revisit to deal with mito misclassification
    // find the best pltd component
    uint32_t i, j, m, n, p_b, p_b1, genid, *mito_gen, *pltd_gen;
    uint64_t x;
    double p_s, p_s1;
    og_component_t *component;
    kvec_t(uint64_t) gen_list;
    OG_TYPE_t og_type;
    
    n = component_v->n;
    kv_init(gen_list);
    // collect best gene hits in each subgraph
    // pack in the component id
    for (i = 0; i < n; ++i) {
        component = &component_v->a[i];
        m = gen_list.n;
        kv_pushn(uint64_t, gen_list, component->g, component->ng);
        for (j = m; j < gen_list.n; ++j) {
            x = gen_list.a[j];
            gen_list.a[j] = (x & 0xFFFFFFFF00000000ULL) | ((uint32_t) x << 16 | i);
        }
    }
    
    if (gen_list.n == 0)
        return;

    // sort by gid and then score - descending order
    qsort(gen_list.a, gen_list.n, sizeof(uint64_t), u64_cmpfunc_r);
    
    MYCALLOC(mito_gen, n);
    MYCALLOC(pltd_gen, n);
    genid = gen_list.a[0] >> 32;
    for (i = 0, j = 0, m = gen_list.n; i < m; ++i) {
        if (gen_list.a[i] >> 32 != genid ||
                i == m - 1) {
            og_type = (gen_list.a[i] >> 32) & 0x3;
            if (og_type != OG_MITO &&
                    og_type != OG_PLTD)
                continue;
            if (i == j || 
                   (double) ((gen_list.a[j+1] >> 16) & 0xFFFFULL) < 
                   (double) ((gen_list.a[j] >> 16) & 0xFFFFULL) * 0.8) {
                x = gen_list.a[j] & 0xFFFFULL;
                if (og_type == OG_MITO)
                    ++mito_gen[x];
                else
                    ++pltd_gen[x];
            }
            genid = gen_list.a[i] >> 32;
            j = i;
        }
    }

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
            if (pltd_gen[i] > mito_gen[i] * PLTD_TO_MITO_FST[1])
                continue;
            if (component->score > component->sscore * PLTD_TO_MITO_FST[1])
                continue;
            if (component->score < component->sscore * PLTD_TO_MITO_FST[0] ||
                    (component->len < COMMON_MIN_PLTD_SIZE ||
                     component->len > COMMON_MAX_PLTD_SIZE)) {
                double tmp = component->score;
                component->score = component->sscore;
                component->sscore = tmp;
                component->type = component->score > .0? OG_MITO : OG_UNCLASSIFIED;
                if (verbose > 0)
                    fprintf(stderr, "[M::%s] change subgraph organelle type annotation: PLTD [S=%.3f G=%u] -> MITO [S=%.3f G=%u]\n",
                            __func__, component->sscore, pltd_gen[i], component->score, mito_gen[i]);
            }
        }
    }

    free(mito_gen);
    free(pltd_gen);
    kv_destroy(gen_list);
}

double *get_sequence_annot_score(hmm_annot_db_t *annot_db, asg_t *asg, int no_trn, int no_rrn, double max_eval,
        int n_core, int verbose)
{
    uint32_t i, j, k, n, n_seg, n_gene, m_gene, gid;
    uint32_t *asg_snameid_map;
    hmm_annot_t *annots, *annot;
    double *annot_score, *gene_score, *a_s;

    n_gene = annot_db->n;
    n_seg = asg->n_seg;

    if (n_gene == 0) return 0;
    if (n_core == 0) n_core = INT32_MAX;
    
    m_gene = annot_db->n_gene;

    // sort annotations by sid - og_type - gid - score
    hmm_annot_db_sort(annot_db, ORDER_SID_OG);
    // build asg sname to annot_db sname id map
    MYMALLOC(asg_snameid_map, n_seg);
    for (i = 0; i < n_seg; ++i)
        asg_snameid_map[i] = hmm_annot_sname2id(annot_db, asg->seg[i].name);

    // calculate annotation score for each seg
    MYCALLOC(annot_score, n_seg * 4);
    MYMALLOC(gene_score, m_gene * 4);
    for (i = 0; i < n_seg; ++i) {
        annots = hmm_annot_db_index_query(annot_db, asg_snameid_map[i], &n);
        if (!annots) continue;
        
        MYBZERO(gene_score, m_gene * 4);
        for (j = 0; j < n; ++j) {
            annot = &annots[j];
            if (annot->evalue > max_eval ||
                    (no_trn && is_trn(annot)) ||
                    (no_rrn && is_rrn(annot)))
                continue;
            gid = annot->og_type * m_gene + annot->gid;
            if (gene_score[gid] < annot->score)
                gene_score[gid] = annot->score;
        }

        // sort annotation scores for each og_type
        // and get the sum annotation score
        a_s = &annot_score[i * 4];
        for (j = 0; j < 4; ++j) {
            qsort(&gene_score[j * m_gene], m_gene, sizeof(double), dbl_cmpfunc_r);
            n = MIN((uint32_t) n_core, m_gene);
            for (k = 0 ; k < n; ++k)
                a_s[j] += gene_score[j * m_gene + k];
        }
        
        if (verbose > 0)
            fprintf(stderr, "[M::%s] sequence %s: size, %d; mito score, %.3f; pltd score, %.3f; mini score, %.3f\n",
                    __func__, asg->seg[i].name, asg->seg[i].len, a_s[OG_MITO], a_s[OG_PLTD], a_s[OG_MINI]);
    }
    
    free(gene_score);
    free(asg_snameid_map);

    return annot_score;
}

og_component_v *annot_sequence_og_type(hmm_annot_db_t *annot_db, asg_t *asg, int no_trn, int no_rrn, double max_eval,
        int n_core, int min_len, int min_score, int fix_og, int verbose)
{
    uint32_t i, j, k, n, n_seg, n_gene, m_gene, gid, imax, smax;
    uint32_t *asg_snameid_map;
    hmm_annot_t *annots, *annot;
    double a_s[4], *annot_score;
    OG_TYPE_t og_t;
    og_component_t *component;
    og_component_v *component_v;

    n_gene = annot_db->n;
    n_seg = asg->n_seg;

    if (n_gene == 0) return 0;
    if (n_core == 0) n_core = INT32_MAX;

    m_gene = annot_db->n_gene;

    // sort annotations by sid - og_type - gid - score
    hmm_annot_db_sort(annot_db, ORDER_SID_OG);
    // build asg sname to annot_db sname id map
    MYMALLOC(asg_snameid_map, n_seg);
    for (i = 0; i < n_seg; ++i)
        asg_snameid_map[i] = hmm_annot_sname2id(annot_db, asg->seg[i].name);

    MYMALLOC(annot_score, m_gene * 4);
    MYCALLOC(component_v, 1);
    MYCALLOC(component_v->a, n_seg);
    component_v->n = n_seg;
    component_v->m = n_seg;
    for (i = 0; i < n_seg; ++i) {
        component_v->a[i].type = (uint32_t) OG_UNCLASSIFIED;
        if (asg->asmg->vtx[i].del)
            continue;

        MYBZERO(annot_score, m_gene * 4);
        annots = hmm_annot_db_index_query(annot_db, asg_snameid_map[i], &n);
        if (annots) {
            for (k = 0; k < n; ++k) {
                annot = &annots[k];
                if (annot->evalue > max_eval ||
                        (no_trn && is_trn(annot)) ||
                        (no_rrn && is_rrn(annot)))
                    continue;
                gid = annot->og_type * m_gene + annot->gid;
                if (annot_score[gid] < annot->score)
                    annot_score[gid] = annot->score;
            }
        }

        // sort annotation scores for each og_type
        // and get the sum annotation score
        MYBZERO(a_s, 4);
        for (j = 0; j < 4; ++j) {
            qsort(&annot_score[j * m_gene], m_gene, sizeof(double), dbl_cmpfunc_r);
            n = MIN((uint32_t) n_core, m_gene);
            for (k = 0 ; k < n; ++k)
                a_s[j] += annot_score[j * m_gene + k];
        }
        imax = max2(a_s, 4, &smax);
        og_t = OG_UNCLASSIFIED;
        if (a_s[imax] >= min_score)
            // 1, mito; 2, pltd; 3, mini
            og_t = a_s[imax] == a_s[smax]? OG_UNCLASSIFIED : imax;

        if (og_t != OG_UNCLASSIFIED) {
            // make seg list
            uint32_t *comp_v;
            MYMALLOC(comp_v, 1);
            comp_v[0] = i;
            // make gene list
            // get all genes in the component
            kvec_t(uint64_t) comp_g;
            kv_init(comp_g);
            annots = hmm_annot_db_index_query(annot_db, asg_snameid_map[i], &n);
            if (annots) {
                for (k = 0; k < n; ++k) {
                    annot = &annots[k];
                    if (annot->evalue > max_eval ||
                            (no_trn && is_trn(annot)) ||
                            (no_rrn && is_rrn(annot)))
                        continue;
                    // gid:30 | OG_TYPE:2 | score:32
                    kv_push(uint64_t, comp_g,
                            ((uint64_t) annot->gid << 2 | annot->og_type) << 32 |
                            (uint32_t) annot->score);
                }
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
            component = &component_v->a[i];
            component->type = (uint32_t) og_t;
            component->score = a_s[imax];
            component->sscore = a_s[smax];
            component->len = asg->seg[i].len;
            component->nv = 1;
            component->v = comp_v;
            component->ng = comp_g.n;
            component->g = comp_g.a;
            component->asmg = 0;
        }

        if (verbose > 0)
            fprintf(stderr, "[M::%s] sequence %s: size, %d; mito score, %.3f; pltd score, %.3f; mini score, %.3f; classification, %d\n",
                    __func__, asg->seg[i].name, asg->seg[i].len, a_s[OG_MITO], a_s[OG_PLTD], a_s[OG_MINI], og_t);
    }

    free(annot_score);
    free(asg_snameid_map);
    
    if (fix_og) fix_og_misclassification(component_v, verbose);

    return component_v;
}

// this collects top n_core unique genes from each subgraph and calculates total scores
// TODO consider the bias introduced by assembly graph size such as nuclear regions with numts
og_component_v *annot_subgraph_og_type(hmm_annot_db_t *annot_db, asg_t *asg, int no_trn, int no_rrn, double max_eval, 
        int n_core, int min_len, int min_score, int fix_og, int verbose)
{
    uint32_t i, j, k, n, n_seg, n_gene, m_gene, gid, imax, smax;
    uint32_t *asg_snameid_map;
    hmm_annot_t *annots, *annot;
    double a_s[4], *annot_score;
    int len, nv;
    uint8_t *visited;
    OG_TYPE_t og_t;
    og_component_t *component;
    og_component_v *component_v;
    asmg_t *g;

    n_gene = annot_db->n;
    n_seg = asg->n_seg;

    if (n_gene == 0) return 0;
    if (n_core == 0) n_core = INT32_MAX;

    m_gene = annot_db->n_gene;

    // sort annotations by sid - og_type - gid - score
    hmm_annot_db_sort(annot_db, ORDER_SID_OG);
    // build asg sname to annot_db sname id map
    MYMALLOC(asg_snameid_map, n_seg);
    for (i = 0; i < n_seg; ++i)
        asg_snameid_map[i] = hmm_annot_sname2id(annot_db, asg->seg[i].name);

    // extract each subgraph and do classification
    // TODO deal with graphs when mito and pltd are in the same subgraph
    MYMALLOC(annot_score, m_gene * 4);
    MYCALLOC(component_v, 1);
    MYCALLOC(visited, n_seg);
    for (i = 0; i < n_seg; ++i) {
        if (visited[i] || asg->asmg->vtx[i].del)
            continue;
        g = asg_make_asmg_copy(asg->asmg, 0);
        asmg_subgraph(g, &i, 1, 0, 0, 0, 1);
        MYBZERO(annot_score, m_gene * 4);
        len = nv = 0;
        for (j = 0; j < n_seg; ++j) {
            if (!g->vtx[j].del) {
                annots = hmm_annot_db_index_query(annot_db, asg_snameid_map[j], &n);
                if (annots) {
                    for (k = 0; k < n; ++k) {
                        annot = &annots[k];
                        if (annot->evalue > max_eval ||
                                (no_trn && is_trn(annot)) ||
                                (no_rrn && is_rrn(annot)))
                            continue;
                        gid = annot->og_type * m_gene + annot->gid;
                        if (annot_score[gid] < annot->score)
                            annot_score[gid] = annot->score;
                    }
                }
                ++nv;
                len += g->vtx[j].len;
                visited[j] = 1;
            }
        }

        // sort annotation scores for each og_type
        // and get the sum annotation score
        MYBZERO(a_s, 4);
        for (j = 0; j < 4; ++j) {
            qsort(&annot_score[j * m_gene], m_gene, sizeof(double), dbl_cmpfunc_r);
            n = MIN((uint32_t) n_core, m_gene);
            for (k = 0 ; k < n; ++k)
                a_s[j] += annot_score[j * m_gene + k];
        }
        
        imax = max2(a_s, 4, &smax);
        og_t = OG_UNCLASSIFIED;
        if (len >= min_len || a_s[imax] >= min_score)
            // 1, mito; 2, pltd; 3, mini
            og_t = a_s[imax] == a_s[smax]? OG_UNCLASSIFIED : imax;

        if (og_t == OG_UNCLASSIFIED) {
            asmg_destroy(g);
            continue;
        }
        
        // make seg list
        uint32_t *comp_v;
        MYMALLOC(comp_v, nv);
        k = 0;
        for (j = 0; j < n_seg; ++j)
            if (!g->vtx[j].del)
                comp_v[k++] = j;
        // make gene list
        // get all genes in the component
        kvec_t(uint64_t) comp_g;
        kv_init(comp_g);
        for (j = 0; j < nv; ++j) {
            annots = hmm_annot_db_index_query(annot_db, asg_snameid_map[comp_v[j]], &n);
            if (annots) {
                for (k = 0; k < n; ++k) {
                    annot = &annots[k];
                    if (annot->evalue > max_eval ||
                            (no_trn && is_trn(annot)) ||
                            (no_rrn && is_rrn(annot)))
                        continue;
                    // gid:30 | OG_TYPE:2 | score:32
                    kv_push(uint64_t, comp_g,
                            ((uint64_t) annot->gid << 2 | annot->og_type) << 32 |
                            (uint32_t) annot->score);
                }
            }
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
        component->type = (uint32_t) og_t;
        component->score = a_s[imax];
        component->sscore = a_s[smax];
        component->len = len;
        component->nv = nv;
        component->v = comp_v;
        component->ng = comp_g.n;
        component->g = comp_g.a;
        component->asmg = g;
        
        if (verbose > 0)
            fprintf(stderr, "[M::%s] subgraph seeding from %s: segs, %d; size, %d; mito score, %.3f; pltd score, %.3f; mini score, %.3f; classification, %d\n",
                    __func__, asg->seg[i].name, nv, len, a_s[OG_MITO], a_s[OG_PLTD], a_s[OG_MINI], og_t);
    }
    
    free(annot_score);
    free(visited);
    free(asg_snameid_map);

    if (fix_og) fix_og_misclassification(component_v, verbose);

    qsort(component_v->a, component_v->n, sizeof(og_component_t), og_cmpfunc_r);

    return component_v;
}

// this collects top n_core genes from this sequence and calculates total scores for each subgraph
// it might be biased when there are multiple copies of some genes
// especially for subgraphs from repetitive regions
// this function has been replaced by annot_subgraph_og_type
og_component_v *annot_subgraph_og_type1(hmm_annot_db_t *annot_db, asg_t *asg, int no_trn, int no_rrn, double max_eval,
        int n_core, int min_len, int min_score, double **_annot_score, int verbose)
{
    uint32_t i, j, k, n_seg, m_seg, n_gene, sid, gid, imax, smax;
    uint32_t *ann_snameid_map;
    hmm_annot_t *annots, *annot;
    double a_s[4], *annot_score;
    int a_n[4], len, nv;
    uint8_t *visited;
    OG_TYPE_t og_t;
    og_component_t *component;
    og_component_v *component_v;
    asmg_t *g;

    n_gene = annot_db->n;
    annots = annot_db->a;
    n_seg = asg->n_seg;
    m_seg = annot_db->n_seg;

    if (n_gene == 0) return 0;

    // sort annotations by sid - score
    hmm_annot_db_sort(annot_db, ORDER_SID_OG);

    // build annot_db sname to asg sname id map
    MYMALLOC(ann_snameid_map, m_seg);
    for (i = 0; i < m_seg; ++i)
        ann_snameid_map[i] = asg_name2id(asg, annot_db->snames[i]);

    // calculate annotation score for each seg
    sid = annots[0].sid;
    MYBZERO(a_s, 4);
    MYBZERO(a_n, 4);
    MYCALLOC(annot_score, n_seg * 4);
    for (i = 0; i < n_gene; ++i) {
        annot = &annots[i];
        if (annot->sid != sid) {
            // here remember to convert annotation sid to asg sid
            if (ann_snameid_map[sid] != UINT32_MAX)
                memcpy(annot_score+ann_snameid_map[sid]*4, a_s, sizeof(double)*4);
            sid = annot->sid;
            MYBZERO(a_s, 4);
            MYBZERO(a_n, 4);
        }
        if (annot->evalue > max_eval || (no_trn && is_trn(annot)) || (no_rrn && is_rrn(annot)))
            continue;
        if (a_n[annot->og_type]++ < n_core)
            a_s[annot->og_type] += annot->score;
    }
    // here remember to convert annotation sid to asg sid
    if (ann_snameid_map[sid] != UINT32_MAX)
        memcpy(annot_score+ann_snameid_map[sid]*4, a_s, sizeof(double)*4);

    // extract each subgraph and do classification
    // TODO deal with graphs when mito and pltd are in the same subgraph
    MYCALLOC(component_v, 1);
    MYCALLOC(visited, n_seg);
    for (i = 0; i < n_seg; ++i) {
        if (visited[i] || asg->asmg->vtx[i].del)
            continue;
        g = asg_make_asmg_copy(asg->asmg, 0);
        asmg_subgraph(g, &i, 1, 0, 0, 0, 1);
        MYBZERO(a_s, 4);
        len = nv = 0;
        for (j = 0; j < n_seg; ++j) {
            if (!g->vtx[j].del) {
                for (k = 0; k < 4; ++k)
                    a_s[k] += annot_score[j*4+k];
                ++nv;
                len += g->vtx[j].len;
                visited[j] = 1;
            }
        }
        imax = max2(a_s, 4, &smax);
        og_t = OG_UNCLASSIFIED;
        if (len >= min_len || a_s[imax] >= min_score)
            // 1, mito; 2, pltd; 3, mini
            og_t = a_s[imax] == a_s[smax]? OG_UNCLASSIFIED : imax;

        if (og_t == OG_UNCLASSIFIED) {
            asmg_destroy(g);
            continue;
        }
        
        // make seg list
        uint32_t *comp_v;
        MYMALLOC(comp_v, nv);
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
            sid = ann_snameid_map[annot->sid];
            if (sid == UINT32_MAX ||
                    g->vtx[sid].del ||
                    annot->evalue > max_eval ||
                    (no_trn && is_trn(annot)) || 
                    (no_rrn && is_rrn(annot)))
                continue;
            // gid:30 | OG_TYPE:2 | score:32
            kv_push(uint64_t, comp_g,
                    ((uint64_t) annot->gid << 2 | annot->og_type) << 32 | 
                    (uint32_t) annot->score);
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
        component->type = (uint32_t) og_t;
        component->score = a_s[imax];
        component->sscore = a_s[smax];
        component->len = len;
        component->nv = nv;
        component->v = comp_v;
        component->ng = comp_g.n;
        component->g = comp_g.a;
        component->asmg = g;

        if (verbose > 0)
            fprintf(stderr, "[M::%s] subgraph seeding from %s: segs, %d; size, %d; mito score, %.3f; pltd score, %.3f; mini score, %.3f; classification, %d\n",
                    __func__, asg->seg[i].name, nv, len, a_s[OG_MITO], a_s[OG_PLTD], a_s[OG_MINI], og_t);
    }
    
    free(visited);
    if (_annot_score)
        *_annot_score = annot_score;
    else
        free(annot_score);
    free(ann_snameid_map);

    fix_og_misclassification(component_v, verbose);

    qsort(component_v->a, component_v->n, sizeof(og_component_t), og_cmpfunc_r);

    return component_v;
}

typedef struct {
    int index;
    double val;
    int size;
    int clust;
    int gene_num[4];
    OG_TYPE_t og_type;
} clust_dp_t;

typedef struct {
    int clust;
    int n_dps;
    int *dps;
    int size;
    double mean;
    double og_score[4];
    OG_TYPE_t og_type;
} clust_t;

static int dpVal_cmpfunc(const void *a, const void *b)
{
    double x, y;
    x = ((clust_dp_t *) a)->val;
    y = ((clust_dp_t *) b)->val;
    return (x > y) - (x < y);
}

static int dpIdx_cmpfunc(const void *a, const void *b)
{
    int x, y;
    x = ((clust_dp_t *) a)->index;
    y = ((clust_dp_t *) b)->index;
    return (x > y) - (x < y);
}

#define DBSCAN_EPS 0.25
#define CLUSTV_EPS 0.50

static int dbscan_cluster(clust_dp_t *dps, int n_dps, double eps, double v_eps)
{
    if (n_dps <= 0) return 0;
    
    qsort(dps, n_dps, sizeof(clust_dp_t), dpVal_cmpfunc);
    
    int i, n;
    double vals;
    dps[0].clust = 0;
    vals = dps[0].val;
    n = 1;
    for (i = 1; i < n_dps; ++i) {
        if (dps[i].val <= dps[i-1].val * (1 + eps) &&
                dps[i].val <= vals / n * (1 + v_eps)) {
            dps[i].clust = dps[i-1].clust;
            vals += dps[i].val;
            n += 1;
        } else {
            dps[i].clust = dps[i-1].clust + 1;
            vals = dps[i].val;
            n = 1;
        }
    }
    
    int n_clust = dps[n_dps-1].clust + 1;
    
    qsort(dps, n_dps, sizeof(clust_dp_t), dpIdx_cmpfunc);

    return n_clust;
}

static clust_t *make_cluster(clust_dp_t *dps, int n_dps, int n_clust, int min_len)
{
    int i, c;
    clust_t *clusts;
    
    MYCALLOC(clusts, n_clust);
    // get dps number for each cluster
    for (i = 0; i < n_dps; ++i)
        ++clusts[dps[i].clust].n_dps;
    for (i = 0; i < n_clust; ++i) {
        clusts[i].clust = i;
        MYMALLOC(clusts[i].dps, clusts[i].n_dps);
        clusts[i].n_dps = 0;
    }
    for (i = 0; i < n_dps; ++i) {
        c = dps[i].clust;
        clusts[c].dps[clusts[c].n_dps] = dps[i].index;
        clusts[c].mean += dps[i].val;
        clusts[c].size += dps[i].size;
        clusts[c].n_dps++;
    }
    for (i = 0; i < n_clust; ++i)
        clusts[i].mean /= clusts[i].n_dps;
    
    /***
    // for each cluster determine if is a major or minor cluster
    // major cluster tends to be low-copy sequences
    // minor cluster is allowed to have higher copy numbers
    // decision made mainly based on the sequence size
    // mark a major cluster if 
    // the cluster size is larger than 10% of the total size
    // and at least one sequence is larger than min_len
    // mark major cluster using the last bit of 'clust'
    // TODO estiamte size proportion (10%) with sequences of similar coverage
    int j;
    double size = 0.;
    for (i = 0; i < n_clust; ++i)
        size += clusts[i].size;
    size *= 0.1;
    for (i = 0; i < n_clust; ++i) {
        clusts[i].clust <<= 1;
        if (clusts[i].size < size)
            continue;
        for (j = 0; j < clusts[i].n_dps; ++j) {
            if (dps[clusts[i].dps[j]].size >= min_len) {
                clusts[i].clust |= 1; // mark major cluster
                break;
            }
        }
    }
    **/

    return clusts;
}

#define LOG1_5 0.405465108108
#define LOG2_5 0.916290731874
#define LOG3_5 1.252762968495
#define LOG4_5 1.504077396776

KDQ_INIT(uint64_t)

typedef kdq_t(uint64_t) kdq_u64_t;

static void slim_graph(asg_t *asg, og_component_v *sequence_og, og_component_t *component_g, clust_dp_t *comp_dps, OG_TYPE_t og_target, 
        OG_TYPE_t *og_seeds, double c_mean, int max_r_len, og_component_v *components, int verbose)
{
    uint32_t i, j, k, l, v, w, r, nv, na, len, gen, cov, gid, n_vtx, imax, smax, *comp_v;
    uint8_t *dels, *visited;
    double score[4];
    asmg_t *asmg;
    asmg_arc_t *av;
    og_component_t *component_s, *component;

    asmg = asg_make_asmg_copy(component_g->asmg, 0);
    n_vtx = asmg->n_vtx;
    comp_v = component_g->v;
    nv = component_g->nv;

    // do graph slimming
    if (verbose > 1) {
        fprintf(stderr, "[M::%s] subgraph slimming for organelle [%s] with average sequence coverage %.3f\n", 
                __func__, OG_TYPES[og_target], c_mean);
        if (verbose > 2) {
            asmg_t *g = asg->asmg;
            asg->asmg = asmg;
            asg_print(asg, stderr, 1);
            asg->asmg = g;
        }
    }

    // find sequences to remove
    MYCALLOC(dels, n_vtx);
    for (i = 0; i < nv; ++i)
        if (og_seeds[i] != og_target)
            dels[comp_v[i]] = 1;

    // bring back the sequence if in a critical path
    // TODO improve the recalling strategy
    // currently bring repeat sequences back if attached to any non-deleted
    /***
    while (1) {
        int recall = 0;
        for (i = 0; i < nv; ++i) {
            v = comp_v[i];
            if (!dels[v] || asmg->vtx[v].len > max_r_len)
                continue; // not a repeat
            // check incoming and outgoing arcs to see if any non-deleted
            int done = 0;
            for (j = 0; j < 2 && !done; ++j) {
                av = asmg_arc_a(asmg, v<<1 | j);
                na = asmg_arc_n(asmg, v<<1 | j);
                for (k = 0; k < na; ++k) {
                    if (!av[k].del && !dels[av[k].w>>1]) {
                        dels[v] = 0;
                        ++recall;
                        done = 1;
                        break;
                    }
                }
            }
        }
        if (!recall) break;
    }
    **/
    
    // here is an improved strategy
    // bring a repeat sequence back if both ends connect to a sequence
    // through a path of all repeats
    // and the distance is no greater than max_r_len
    // repeat this process until no sequence is added
    int recall, max_r, *dist, *flag;
    kdq_u64_t *q;

    MYMALLOC(dist, n_vtx << 1);
    MYMALLOC(flag, n_vtx << 1);
    q = kdq_init(uint64_t, 0);
    while (1) {
        MYBZERO(dist, n_vtx << 1);
        for (i = 0; i < nv; ++i) {
            if (dels[comp_v[i]])
                // skip dels
                continue;
            // the maximum distance to radiate
            max_r = asmg->vtx[comp_v[i]].len;
            if (max_r > max_r_len)
                max_r = max_r_len;
            for (k = 0; k < 2; ++k) {
                uint64_t source = comp_v[i] << 1 | k;
                // find all reachable repeats from 'source'
                MYBZERO(flag, n_vtx << 1);
                kdq_push(uint64_t, q, source << 32 | 0);
                while (kdq_size(q) > 0) {
                    uint64_t x = *kdq_shift(uint64_t, q);
                    v = x >> 32;
                    r = (uint32_t) x;
                    flag[v] = 1;
                    dist[v] = source << 1 | 1;
                    av = asmg_arc_a(asmg, v);
                    na = asmg_arc_n(asmg, v);
                    for (j = 0; j < na; ++j) {
                        if (av[j].del)
                            continue;
                        w = av[j].w;
                        if (!flag[w] &&
                                r <= av[j].ls + max_r &&
                                asmg->vtx[w>>1].len <= max_r)
                            // add non-visited repeats
                            kdq_push(uint64_t, q, (uint64_t) w << 32 | (r + asmg->vtx[w>>1].len - av[j].ls));
                    }
                }
            }
        }

        recall = 0;
        for (i = 0; i < nv; ++i) {
            v = comp_v[i];
            if (dels[v] &&
                    asmg->vtx[v].len <= max_r_len &&
                    dist[v<<1] &&
                    dist[v<<1 | 1]) {
                dels[v] = 0;
                ++recall;
                if (verbose > 2)
                    fprintf(stderr, "[M::%s] recall repeat %s: %s%c <-> %s%c\n", __func__,
                            asg->seg[v].name,
                            asg->seg[dist[v<<1]>>2].name, "+-"[(dist[v<<1]>>1)&1],
                            asg->seg[dist[v<<1|1]>>2].name, "+-"[(dist[v<<1|1]>>1)&1]);
            }
        }
        if (!recall) break;
    }
    free(dist);
    free(flag);
    kdq_destroy(uint64_t, q);

    for (i = 0; i < nv; ++i) {
        v = comp_v[i];
        if (dels[v])
            asmg_vtx_del(asmg, v, 1);
    }

    if (verbose > 2) {
        asmg_t *g = asg->asmg;
        asg->asmg = asmg;
        asg_print(asg, stderr, 1);
        asg->asmg = g;
    }

    uint64_t cleaned = 1;
    while (cleaned) {
        // simple clean
        cleaned = 0;
        cleaned += asmg_pop_bubble(asmg, max_r_len, 0, 0, 1, 0, verbose);
        cleaned += asmg_remove_weak_crosslink(asmg, 0.3, 10, 0, verbose);
        cleaned += asmg_drop_tip(asmg, INT32_MAX, max_r_len, 1, 0, verbose);
    }
    // add cleaned vtx to dels for integrity
    for (i = 0; i < nv; ++i) {
        v = comp_v[i];
        if (asmg->vtx[v].del)
            dels[v] = 1;
    }

    if (verbose > 2) {
        asmg_t *g = asg->asmg;
        asg->asmg = asmg;
        asg_print(asg, stderr, 1);
        asg->asmg = g;
    }

    // minimum subgraph size to keep
    // set as 10% of the total sequence
    double m_size = 0;
    for (i = 0; i < nv; ++i) {
        v = comp_v[i];
        if (!asmg->vtx[v].del)
            m_size += asmg->vtx[v].len;
    }
    m_size *= 0.1;

    // now process each subgraph
    MYCALLOC(visited, n_vtx);
    for (i = 0; i < nv; ++i) {
        v = comp_v[i];
        if (visited[v] || asmg->vtx[v].del)
            continue;
        
        asmg_t *g = asg_make_asmg_copy(asmg, 0);
        asmg_subgraph(g, &v, 1, 0, 0, 0, 1);

        kvec_t(uint32_t) comp_s;
        kv_init(comp_s);
        len = gen = 0;
        for (j = 0; j < nv; ++j) {
            w = comp_v[j];
            if (g->vtx[w].del)
                continue;
            kv_push(uint32_t, comp_s, w);
            len += g->vtx[w].len;
            gen += comp_dps[j].gene_num[og_target];
            visited[w] = 1;
        }
        if (len < m_size || gen == 0) {
            // not a good subgraph
            asmg_destroy(g);
            kv_destroy(comp_s);
            continue;
        }
        
        // add subgraph component
        // adjust repeat sequence coverage if necessary
        // also fix the arc coverage
        for (j = 0; j < nv; ++j) {
            w = comp_v[j];
            if (g->vtx[w].del ||
                    og_seeds[j] == og_target ||
                    g->vtx[w].len >= max_r_len ||
                    g->vtx[w].cov < c_mean * 3.5)
                continue;
            int n_del, n_arc;
            n_del = n_arc = 0;
            for (k = 0; k < 2; ++k) {
                av = asmg_arc_a(asmg, w << 1 | k);
                na = asmg_arc_n(asmg, w << 1 | k);
                for (l = 0; l < na; ++l) {
                    if (dels[av[l].w>>1])
                        ++n_del;
                    if (!av[l].del)
                        ++n_arc;
                }
            }
            if (!n_del) continue;
            cov = g->vtx[w].cov;
            g->vtx[w].cov = c_mean * n_arc / 2.0;
            // fix arc coverage
            for (k = 0; k < 2; ++k) {
                av = asmg_arc_a(asmg, w << 1 | k);
                na = asmg_arc_n(asmg, w << 1 | k);
                for (l = 0; l < na; ++l) {
                    if (av[l].del) continue;
                    if (av[l].cov > cov)
                        av[l].cov = cov;
                }
            }
            if (verbose > 2)
                fprintf(stderr, "[M::%s] sequence %s coverage adjusted: %u -> %u\n", 
                        __func__, asg->seg[w].name, cov, g->vtx[w].cov);
        }

        // make seg list
        MYREALLOC(comp_s.a, comp_s.n);

        // make gene list
        kvec_t(uint64_t) comp_g;
        kv_init(comp_g);
        // collect genes from sequences in the subgraph
        for (j = 0; j < nv; ++j) {
            w = comp_v[j];
            if (g->vtx[w].del)
                continue;
            component_s = &sequence_og->a[w];
            kv_pushn(uint64_t, comp_g, component_s->g, component_s->ng);
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
        comp_g.n = k;
        MYREALLOC(comp_g.a, k);

        // calculate subgraph annotation scores
        uint64_t *a_s = comp_g.a;
        MYBZERO(score, 4);
        for (j = 0, k = comp_g.n; j < k; ++j)
            score[(a_s[j] >> 32) & 0x3] += (uint32_t) a_s[j];
        imax = max2(score, 4, &smax);

        // update component vector
        kv_pushp(og_component_t, *components, &component);
        component->type = (uint32_t) og_target;
        component->score = score[imax];
        component->sscore = score[smax];
        component->len = len;
        component->nv = comp_s.n;
        component->v = comp_s.a;
        component->ng = comp_g.n;
        component->g = comp_g.a;
        component->asmg = g;

        if (verbose > 0)
            fprintf(stderr, "[M::%s] subgraph seeding from %s: segs, %lu; size, %d; mito score, %.3f; pltd score, %.3f; mini score, %.3f; classification, %d\n",
                    __func__, asg->seg[v].name, comp_s.n, len, score[OG_MITO], score[OG_PLTD], score[OG_MINI], og_target);

    }

    free(dels);
    free(visited);
    asmg_destroy(asmg);
}

static void find_coverage_bound_in_mixture_graph(clust_t *clusts, uint32_t n_clust, clust_dp_t *comp_dps, uint32_t nv, uint32_t max_mito_size, 
        uint32_t max_pltd_size, double fold_thresh, double *min_mito, double *max_mito, double *min_pltd, double *max_pltd, int verbose)
{
    return;
}

static uint32_t find_seeds_in_pure_graph(clust_t *clusts, uint32_t n_clust, clust_dp_t *comp_dps, uint32_t nv, OG_TYPE_t og_t, double min_mean, 
        double max_mean, double fold_thresh, uint32_t min_size, uint32_t max_size, double *_c_mean, OG_TYPE_t *og_seeds, int verbose)
{
    uint32_t i, j, c, v, seed, n_seeds, l_seeds, *seed_clust, *gseq_clust;
    uint64_t *gene_clust;
    double c_mean, min_mean1, max_mean1;

    // count number of best genes in each coverage cluster
    MYCALLOC(gseq_clust, n_clust);
    MYCALLOC(gene_clust, n_clust);
    for (i = 0; i < n_clust; ++i) {
        for (j = 0; j < clusts[i].n_dps; ++j) {
            v = clusts[i].dps[j];
            gene_clust[i] += comp_dps[v].gene_num[og_t];
            if (comp_dps[v].gene_num[og_t] > 0)
                gseq_clust[i] += comp_dps[v].size;
        }
        gene_clust[i] <<= 32;
        gene_clust[i] |= i;
        // a simple rule to decide if
        // all sequences in a cluster should be added together
        if (gseq_clust[i] > 0.5 * clusts[i].size)
            gseq_clust[i] = clusts[i].size;
    }
    // sort clusters by number of best genes in descending order
    qsort(gene_clust, n_clust, sizeof(uint64_t), u64_cmpfunc_r);
    
    // collect seeds in each cluster
    MYCALLOC(seed_clust, n_clust);
    min_mean1 = max_mean1 = 0;
    n_seeds = l_seeds = 0;
    for (i = 0; i < n_clust; ++i) {
        if ((gene_clust[i] >> 32) == 0)
            // no gene
            break;
        c = (uint32_t) gene_clust[i]; // cluster index
        if (clusts[c].og_type != og_t)
            // not the right organelle type
            continue;
        c_mean = clusts[c].mean;
        if (c_mean < min_mean && c_mean > max_mean)
            // coverage not in range
            continue;
        if (l_seeds + gseq_clust[c] > max_size)
            // too much sequence
            continue;
        seed = 0;
        if (n_seeds == 0) {
            // major cluster
            // mean coverage not set yet
            min_mean1 = c_mean;
            max_mean1 = c_mean;
            seed = 1;
        } else {
            // secondary clusters
            if (gseq_clust[c] >= min_size) {
                // this is a large cluster
                // make it a seed only when
                // the fold difference is within the range
                if (c_mean >= min_mean1 && c_mean <= max_mean1) {
                    // mean coverage is within the old range
                    seed = 1;
                } else if (fabs(log(min_mean1/c_mean)) <= fold_thresh &&
                        fabs(log(max_mean1/c_mean)) <= fold_thresh) {
                    // mean coverage is within the new range
                    // update coverage
                    if (c_mean < min_mean1)
                        min_mean1 = c_mean;
                    if (c_mean > max_mean1)
                        max_mean1 = c_mean;
                    seed = 1;
                }
            } else {
                // this is a small cluster
                // make it a seed
                // but do not change the coverage
                seed = 1;
            }
        }
        if (seed) {
            seed_clust[c] = 1;
            n_seeds += 1;
            l_seeds += gseq_clust[c];
        }
    }
    
    MYBZERO(og_seeds, nv);
    // process clusters
    int all_seq;
    for (i = 0; i < n_clust; ++i) {
        if (!seed_clust[i])
            continue;
        // not always to add all sequences
        all_seq = clusts[i].size == gseq_clust[i];
        for (j = 0; j < clusts[i].n_dps; ++j)
            if (all_seq || 
                    comp_dps[clusts[i].dps[j]].gene_num[og_t] > 0)
                og_seeds[clusts[i].dps[j]] = og_t;
    }

    // process each unclassified individual sequence
    uint32_t ext_l, ext_n;
    OG_TYPE_t og_t1, *og_seeds1;
    og_t1 = OG_UNCLASSIFIED;
    if (og_t == OG_MITO)
        og_t1 = OG_PLTD;
    else if (og_t == OG_PLTD)
        og_t1 = OG_MITO;
    MYCALLOC(og_seeds1, nv);
    ext_l = ext_n = 0;
    for (i = 0; i < nv; ++i) {
        c_mean = comp_dps[i].val;
        if (!og_seeds[i] && 
                (og_t1 == OG_UNCLASSIFIED || 
                 comp_dps[i].gene_num[og_t1] == 0 ||
                 comp_dps[i].gene_num[og_t] > 0) &&
                c_mean >= min_mean &&
                c_mean <= max_mean &&
                fabs(log(min_mean1/c_mean)) <= fold_thresh) {
            ext_l += comp_dps[i].size;
            ext_n += 1;
            og_seeds1[i] = og_t;
        }
    }
    // only add extra sequences if total size within threshold
    if (l_seeds + ext_l <= max_size) {
        for (i = 0; i < nv; ++i)
            if (og_seeds1[i])
                og_seeds[i] = og_seeds1[i];
        n_seeds += ext_n;
        l_seeds += ext_l;
    }

    free(gseq_clust);
    free(gene_clust);
    free(seed_clust);
    free(og_seeds1);
    
    if (verbose > 2)
        fprintf(stderr, "[M::%s] [%s] - seeds number: %u; seeds length: %u; seeds coverage: %.3f - %.3f\n", 
                __func__, OG_TYPES[og_t], n_seeds, l_seeds, min_mean1, max_mean1);

    if (_c_mean)
        *_c_mean = min_mean1;
    
    return l_seeds;
}

// this is a improved version of annot_subgraph_og_type
// to deal with graphs with contamination sequences
// e.g., pltd-mito mixtures
// e.g., graphs with sequence groups of abnormal coverage
// the core idea is to
// make sequence clusters based on coverage as implemented in dbscan_cluster
// and then clean graph based on sequence clusters as implemented in slim_graph
og_component_v *asg_annotation(hmm_annot_db_t *annot_db, asg_t *asg, int no_trn, int no_rrn, double max_eval,
        int n_core, int min_len, int min_score, int fix_og, int verbose)
{
    uint32_t i, j, k, n, m, v, nv, gid, imax, smax, score, n_seg, n_gene, m_gene, n_clust;
    uint32_t l_seeds[4], n_seeds[4], *comp_v;
    uint64_t *gene_score;
    double *a_s, *annot_score, *seg_score, g_n[4], c_means[4];
    OG_TYPE_t og_t, **og_seeds;
    og_component_t *component_s, *component_g;
    og_component_v *sequence_og, *subgraph_og, *component_v;
    clust_dp_t *comp_dps;
    clust_t *clusts;

    n_gene = annot_db->n;
    n_seg = asg->n_seg;

    if (n_gene == 0) return 0;
    
    m_gene = annot_db->n_gene;

    seg_score = get_sequence_annot_score(annot_db, asg, no_trn, no_rrn, max_eval, 0, verbose);
    sequence_og = annot_sequence_og_type(annot_db, asg, no_trn, no_rrn, max_eval, n_core, min_len, min_score, 0, verbose);
    subgraph_og = annot_subgraph_og_type(annot_db, asg, no_trn, no_rrn, max_eval, n_core, min_len, min_score, 0, verbose);

    MYCALLOC(annot_score, m_gene * 4);
    // collect best annotation score for each gene across the entire graph
    for (i = 0; i < n_seg; ++i) {
        component_s = &sequence_og->a[i];
        gene_score = component_s->g;
        for (j = 0; j < component_s->ng; ++j) {
            gid = gene_score[j] >> 34;
            og_t = (gene_score[j] >> 32) & 0x3;
            score = (uint32_t) gene_score[j];
            
            gid += m_gene * og_t;
            if (annot_score[gid] < (double) score)
                annot_score[gid] = score;
        }
    }

    // check if each subgraph is a mito-pltd mixture
    // separate organelle genomes for mixtures
    double g_diff = 0.85; // relaxation for gene socres
    MYCALLOC(component_v, 1);
    for (i = 0; i < subgraph_og->n; ++i) {
        component_g = &subgraph_og->a[i];
        comp_v = component_g->v;
        nv = component_g->nv;

        MYMALLOC(comp_dps, nv);
        for (j = 0; j < nv; ++j) {
            comp_dps[j].index = j;
            comp_dps[j].val = component_g->asmg->vtx[comp_v[j]].cov;
            comp_dps[j].size = component_g->asmg->vtx[comp_v[j]].len;
            comp_dps[j].clust = -1;
            comp_dps[j].og_type = OG_UNCLASSIFIED;
            MYBZERO(comp_dps[j].gene_num, 4);
        }

        // count genes with best annotation scores
        for (j = 0; j < nv; ++j) {
            component_s = &sequence_og->a[comp_v[j]];
            gene_score = component_s->g;
            // count genes with best annotation scores
            n = 0;
            for (k = 0; k < component_s->ng; ++k) {
                og_t = (gene_score[k] >> 32) & 0x3;
                gid = gene_score[k] >> 34;
                score = (uint32_t) gene_score[k];
                a_s = &annot_score[m_gene * og_t];
                if (score >= min_score &&
                        (double) score >= a_s[gid] * g_diff)
                    // the number of best genes
                    ++comp_dps[j].gene_num[og_t];
            }
        }
        
        // do clustering by sequence overage
        n_clust = dbscan_cluster(comp_dps, nv, DBSCAN_EPS, CLUSTV_EPS);
        clusts = make_cluster(comp_dps, nv, n_clust, min_len);

        // organelle type classification for each coverage cluster
        MYBZERO(l_seeds, 4);
        MYBZERO(n_seeds, 4);
        for (j = 0; j < n_clust; ++j) {
            a_s = clusts[j].og_score;
            MYBZERO(a_s, 4);
            MYBZERO(g_n, 4);
            for (m = 0, n = clusts[j].n_dps; m < n; ++m) {
                v = clusts[j].dps[m];
                for (k = 0; k < 4; ++k) {
                    a_s[k] += seg_score[comp_v[v] * 4 + k];
                    g_n[k] += comp_dps[v].gene_num[k];
                }
            }
            imax = max2(a_s, 4, &smax);
            og_t = a_s[imax] == a_s[smax]? OG_UNCLASSIFIED : imax;
            if (og_t == OG_PLTD && 
                    smax == OG_MITO &&
                    g_n[OG_MITO] > 0 &&
                    (a_s[OG_PLTD] < a_s[OG_MITO] * PLTD_TO_MITO_FST[0] ||
                     (a_s[OG_PLTD] < a_s[OG_MITO] * PLTD_TO_MITO_FST[1] && 
                      clusts[j].size > COMMON_MAX_PLTD_SIZE))) {
                // fix OG_TYPE
                og_t = OG_MITO;
                if (verbose > 2)
                    fprintf(stderr, "[M::%s] OG_TYPE changed from PLTD to MITO: CLUST=%d [M=%.3f N=%d L=%d] [MITO=%.3f PLTD=%.3f MINI=%.3f]\n",
                            __func__, j, clusts[j].mean, clusts[j].n_dps, clusts[j].size, a_s[OG_MITO], a_s[OG_PLTD], a_s[OG_MINI]);
            }
            for (m = 0, n = clusts[j].n_dps; m < n; ++m) {
                v = clusts[j].dps[m];
                if (comp_dps[v].gene_num[og_t] > 0) {
                    l_seeds[og_t] += component_g->asmg->vtx[comp_v[v]].len;
                    n_seeds[og_t] += 1;
                }
            }
            clusts[j].og_type = og_t;
        }

        if (verbose > 2) {
            fprintf(stderr, "[M::%s] subgraph sequence coverage clustering [CLUST=%d]:\n", __func__, n_clust);
            for (j = 0; j < n_clust; ++j) {
                a_s = clusts[j].og_score;
                fprintf(stderr, "[M::%s] CLUST=%d OG=%s [M=%.3f N=%d L=%d] [MITO=%.3f PLTD=%.3f MINI=%.3f]\n", 
                        __func__, j, OG_TYPES[clusts[j].og_type], clusts[j].mean, clusts[j].n_dps, clusts[j].size, 
                        a_s[OG_MITO], a_s[OG_PLTD], a_s[OG_MINI]);
            }
            fprintf(stderr, "[M::%s] subgraph sequence [NV=%d]:\n", __func__, nv);
            qsort(comp_dps, nv, sizeof(clust_dp_t), dpVal_cmpfunc);
            for (j = 0; j < nv; ++j) {
                v = comp_dps[j].index;
                og_t = sequence_og->a[comp_v[v]].type;
                fprintf(stderr, "[M::%s] %s [%s #MG=%d #PG=%d] COV=%.3f LEN=%d CLUST=%d\n", __func__,
                        asg->seg[comp_v[v]].name, OG_TYPES[og_t], comp_dps[j].gene_num[OG_MITO], comp_dps[j].gene_num[OG_PLTD], 
                        comp_dps[j].val, comp_dps[j].size, comp_dps[j].clust);
            }
            qsort(comp_dps, nv, sizeof(clust_dp_t), dpIdx_cmpfunc);
            fprintf(stderr, "[M::%s] subgraph organelle initial classification:\n", __func__);
            for (j = 0; j < 4; ++j)
                fprintf(stderr, "[M::%s] %s: SEED=%u SIZE=%u\n", __func__, OG_TYPES[j], n_seeds[j], l_seeds[j]);
        }

        if (l_seeds[OG_MITO] > 0 && l_seeds[OG_PLTD] > 0) {
            if (l_seeds[OG_MITO] > min_len && l_seeds[OG_PLTD] < min_len) {
                l_seeds[OG_PLTD] = 0;
                n_seeds[OG_PLTD] = 0;
            } else if (l_seeds[OG_MITO] < min_len && l_seeds[OG_PLTD] > min_len) {
                l_seeds[OG_MITO] = 0;
                n_seeds[OG_MITO] = 0;
            }
        }

        MYMALLOC(og_seeds, 4);
        MYCALLOC(og_seeds[0], nv * 4);
        for (j = 1; j < 4; ++j)
            og_seeds[j] = og_seeds[j-1] + nv;
        MYBZERO(c_means, 4);
        if (l_seeds[OG_MITO] > 0 && l_seeds[OG_PLTD] > 0) {
            // is a organelle mixture
            if (verbose > 2)
                fprintf(stderr, "[M::%s] processing subgraph as a organelle mixture\n", __func__);
            double min_mito, max_mito, min_pltd, max_pltd;
            min_mito = min_pltd = 0;
            max_mito = max_pltd = DBL_MAX;
            // TODO implement this subroutine to determine boundary of MITO and PLTD coverage
            find_coverage_bound_in_mixture_graph(clusts, n_clust, comp_dps, nv, COMMON_MAX_MITO_SIZE, COMMON_MAX_PLTD_SIZE, 
                    LOG4_5, &min_mito, &max_mito, &min_pltd, &max_pltd, verbose);
            l_seeds[OG_MITO] = find_seeds_in_pure_graph(clusts, n_clust, comp_dps, nv, OG_MITO, min_mito, max_mito, LOG4_5, 
                    min_len, COMMON_MAX_MITO_SIZE, &c_means[OG_MITO], og_seeds[OG_MITO], verbose);
            l_seeds[OG_PLTD] = find_seeds_in_pure_graph(clusts, n_clust, comp_dps, nv, OG_PLTD, min_pltd, max_pltd, LOG4_5, 
                    min_len, COMMON_MAX_PLTD_SIZE, &c_means[OG_PLTD], og_seeds[OG_PLTD], verbose);
        } else if (l_seeds[OG_MITO] > 0) {
            // pure mito
            if (verbose > 2)
                fprintf(stderr, "[M::%s] processing subgraph as a pure MITO\n", __func__);
            l_seeds[OG_MITO] = find_seeds_in_pure_graph(clusts, n_clust, comp_dps, nv, OG_MITO, 0, DBL_MAX, LOG4_5, 
                    min_len, COMMON_MAX_MITO_SIZE, &c_means[OG_MITO], og_seeds[OG_MITO], verbose);
        } else if (l_seeds[OG_PLTD] > 0) {
            // pure pltd
            // similar to mito
            // but also consider genome size
            if (verbose > 2)
                fprintf(stderr, "[M::%s] processing subgraph as a pure PLTD\n", __func__);
            l_seeds[OG_PLTD] = find_seeds_in_pure_graph(clusts, n_clust, comp_dps, nv, OG_PLTD, 0, DBL_MAX, LOG4_5, 
                    min_len, COMMON_MAX_PLTD_SIZE, &c_means[OG_PLTD], og_seeds[OG_PLTD], verbose);
        } else if (l_seeds[OG_MINI] > 0) {
            // is minicircle
            if (verbose > 2)
                fprintf(stderr, "[M::%s] processing subgraph as a MINICIRCLE\n", __func__);
            l_seeds[OG_MINI] = find_seeds_in_pure_graph(clusts, n_clust, comp_dps, nv, OG_MINI, 0, DBL_MAX, LOG4_5, 
                    min_len, COMMON_MAX_MINICIRCLE_SIZE, &c_means[OG_MINI], og_seeds[OG_MINI], verbose);
        }

        if (l_seeds[OG_MITO] > 0)
            slim_graph(asg, sequence_og, component_g, comp_dps, OG_MITO, og_seeds[OG_MITO],
                    c_means[OG_MITO], min_len, component_v, verbose); 
        if (l_seeds[OG_PLTD] > 0)
            slim_graph(asg, sequence_og, component_g, comp_dps, OG_PLTD, og_seeds[OG_PLTD],
                    c_means[OG_PLTD], min_len, component_v, verbose);
        if (l_seeds[OG_MINI] > 0)
            slim_graph(asg, sequence_og, component_g, comp_dps, OG_MINI, og_seeds[OG_MINI],
                    c_means[OG_MINI], min_len, component_v, verbose);

        free(og_seeds[0]);
        free(og_seeds);
        free(comp_dps);
        for (j = 0; j < n_clust; ++j)
            free(clusts[j].dps);
        free(clusts);
    }
    
    free(seg_score);
    free(annot_score);
    og_component_v_destroy(sequence_og);
    og_component_v_destroy(subgraph_og);

    if (fix_og) fix_og_misclassification(component_v, verbose);
    
    qsort(component_v->a, component_v->n, sizeof(og_component_t), og_cmpfunc_r);

    return component_v;
}

void print_og_classification_summary(asg_t *asg, hmm_annot_db_t *annot_db, og_component_v *og_components, FILE *fo)
{
    if (!og_components) return;
    uint32_t i, j;
    og_component_t *component;
    for (i = 0; i < og_components->n; ++i) {
        component = &og_components->a[i];
        fprintf(fo, "[M::%s] OG component %u \n", __func__, i);
        fprintf(fo, "[M::%s] OG component %u og_type: %s\n", __func__, i, OG_TYPES[component->type]);
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
                    annot_db->gnames[component->g[j]>>34], (uint32_t) component->g[j]);
    }
}

char **asg_vtx_name_list(asg_t *g, uint64_t *_n)
{
    if (_n) *_n = 0;
    uint64_t i, n = 0;
    uint64_t *vlist = asmg_vtx_list(g->asmg, &n);
    char **names = 0;
    if (n) {
        MYMALLOC(names, n);
        for (i = 0; i < n; ++i)
            names[i] = g->seg[vlist[i]].name;
    }
    free(vlist);
    if (_n) *_n = n;
    return names;
}


