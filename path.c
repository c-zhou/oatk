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

const int athaliana_pltd_g71n = 71;
const char *const athaliana_pltd_g71[] = {
    "psbA",  "matK",  "rps16", "psbK",  "psbI", "atpA",  "atpF",  "atpH", "atpI",  "rps2",
    "rpoC2", "rpoC1", "rpoB",  "ycf6",  "psbM", "psbD",  "psbC",  "ycf9", "rps14", "psaB",
    "psaA",  "ycf3",  "rps4",  "ndhJ",  "psbG", "ndhC",  "atpE",  "atpB", "rbcL",  "accD",
    "psaI",  "ycf4",  "cemA",  "petA",  "psbJ", "psbL",  "psbF",  "psbE", "ORF31", "petG",
    "psaJ",  "rpl33", "rps18", "rpl20", "clpP", "psbB",  "psbT",  "psbN", "psbH",  "petB",
    "petD",  "rpoA",  "rps11", "rpl36", "rps8", "rpl14", "rpl16", "rps3", "rpl22", "rps19",
    "ndhF",  "rpl32", "ycf5",  "ndhD",  "psaC", "ndhE",  "ndhG",  "ndhI", "ndhA",  "ndhH",
    "rps15"
};

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

static int u64_cmpfunc(const void *a, const void *b)
{
    return (*(uint64_t *) a > *(uint64_t *) b) - (*(uint64_t *) a < *(uint64_t *) b);
}

KDQ_INIT(uint64_t)

void asg_subgraph(asg_t *asg, uint32_t *seeds, uint32_t n, uint32_t step)
{
    // make a subgraph from given a vertex set and radius
    // all vertices and arcs in other components will be marked as deleted
    if (n == 0) return;
    
    uint32_t i, v, r, nv;
    uint64_t x;
    int8_t *flag;
    kdq_t(uint64_t) *q;
    asmg_t *g;
    asmg_arc_t *av;
   
    if (step == 0)
        step = INT32_MAX;
    g = asg->asmg;
    q = kdq_init(uint64_t, 0);
    MYCALLOC(flag, asmg_vtx_n(g));
    for (i = 0; i < n; ++i) {
        if (seeds[i] < g->n_vtx) {
            kdq_push(uint64_t, q, ((uint64_t)seeds[i]<<1|0)<<32);
            kdq_push(uint64_t, q, ((uint64_t)seeds[i]<<1|1)<<32);
        }
    }
    for (i = 0; i < g->n_vtx; ++i) // mark all segments to be deleted
        g->vtx[i].del = 1;
    for (i = 0; i < g->n_arc; ++i) // mark all arcs to be deleted
        g->arc[i].del = 1;
    while (kdq_size(q) > 0) {
        x = *kdq_shift(uint64_t, q);
        v = x>>32;
        r = (uint32_t) x;
        if (flag[v]) continue; // already visited
        flag[v] = 1;
        g->vtx[v>>1].del = 0;
        if (r < step) {
            nv = asmg_arc_n(g, v);
            av = asmg_arc_a(g, v);
            for (i = 0; i < nv; ++i) {
                av[i].del = 0;
                if (flag[av[i].w] == 0)
                    kdq_push(uint64_t, q, (uint64_t)av[i].w<<32 | (r + 1));
                if (flag[av[i].w^1] == 0)
                    kdq_push(uint64_t, q, (uint64_t)(av[i].w^1)<<32 | (r + 1));
            }
        }
    }
    
    kdq_destroy(uint64_t, q);
    free(flag);
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

#ifdef DEBUG_SEG_COV_BOUND
    fprintf(stderr, "[DEBUG_SEG_COV_BOUND::%s] estimated sequence coverage lower boundary: %.3f\n", __func__, cov_bound);
#endif

    return cov_bound;
}

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
    uint32_t sid;
    uint64_t i, j, v, nv, pv;
    kvec_t(uint64_t) arcs_diff; //, arcs_self;
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
    // kv_init(arcs_self);
    
    for (i = 0; i < 2; ++i) {
        v = s<<1 | i;
        pv = g->idx_p[v];
        nv = asmg_arc_n(g, v);
        av = asmg_arc_a(g, v);
        for (j = 0; j < nv; ++j) {
            if (av[j].del) continue;
            if ((av[j].v >> 1) != (av[j].w >> 1))
                kv_push(uint64_t, arcs_diff, pv+j);
            /***
            else if (av[j].v != av[j].w)
                kv_push(uint64_t, arcs_self, pv+j);
            else if (i == 0)
                // avoid adding v+->v+ and v-->v- twice
                kv_push(uint64_t, arcs_self, pv+j);
            **/
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
        vt->cov = seg->cov / copy;

        // make copies of links
        // arc cov already changed
        for (j = 0; j < arcs_diff.n; ++j) {
            av = &g->arc[arcs_diff.a[j]];
            asmg_arc_add2(g, sid << 1 | (av->v&1), av->w, av->lo, UINT64_MAX, av->cov / copy, av->comp);
        }
        
        /***
        // self arcs
        for (j = 0; j < arcs_self.n; ++j) {
            av = &g->arc[arcs_self.a[j]];
            asmg_arc_add2(g, sid << 1 | (av->v&1), sid << 1 | (av->w&1), av->lo, UINT64_MAX, av->cov / copy, av->comp);
            uint64_t c;
            for (c = 1; c <= i; c++) {
                asmg_arc_add2(g, sid << 1 | (av->v&1), (sid-c) << 1 | (av->w&1), av->lo, UINT64_MAX, av->cov / copy, av->comp);
                asmg_arc_add2(g, (sid-c) << 1 | (av->v&1), sid << 1 | (av->w&1), av->lo, UINT64_MAX, av->cov / copy, av->comp);
            }
        }
        **/
    }

    kv_destroy(arcs_diff);
    // kv_destroy(arcs_self);

    // need to redo sort and index but not graph clean
    asmg_finalize(g, 0);

    // delete the original vtx
    asmg_vtx_del(g, s, 1);

    return;
}

#define EM_MAX_ITER 1000

double estimate_sequence_copy_number_from_coverage(asg_t *asg, int *copy_number, int max_copy)
{
    uint32_t i, n_seg, iter;
    double total_covs, total_lens, avg_cov, new_avg_cov, min_avg_cov;
    asmg_t *g;
    int null_v;
    
    n_seg = asg->n_seg;
    g = asg->asmg;
    min_avg_cov = graph_sequence_coverage_lower_bound(asg, 0.3);
    avg_cov = graph_sequence_coverage_rough(asg);
    avg_cov = MAX(avg_cov, min_avg_cov);
    null_v = copy_number == 0;
    if (null_v) MYCALLOC(copy_number, n_seg);
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        copy_number[i] = MIN(MAX(1, lround((double) asg->seg[i].cov / avg_cov)), max_copy);
    }

    iter = 0;
    while (iter++ < EM_MAX_ITER) {
        total_covs = total_lens = 0;
        for (i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del) continue;
            total_lens += asg->seg[i].len * copy_number[i];
            total_covs += (double) asg->seg[i].len * asg->seg[i].cov / copy_number[i];
        }
        new_avg_cov = MAX(total_covs / total_lens, min_avg_cov);
        if (fabs(new_avg_cov - avg_cov) < FLT_EPSILON) 
            break; // converged
        avg_cov = new_avg_cov;
        for (i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del) continue;
            copy_number[i] = MIN(MAX(1, lround((double) asg->seg[i].cov / avg_cov)), max_copy);
        }
    }

#ifdef DEBUG_SEG_COPY_EST
    fprintf(stderr, "[DEBUG_SEG_COPY_EST::%s] sequence copy number estimation finished in %u iterations with an average sequence coverage %.3f\n",
            __func__, iter, avg_cov);
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        fprintf(stderr, "[DEBUG_SEG_COPY_EST::%s] %s %u %d\n", __func__, asg->seg[i].name, asg->seg[i].cov, copy_number[i]);
    }
#endif
    
    if(null_v) free(copy_number);

    return avg_cov;
}

#define DVAL(var) ((var)->D)
#define FVAL(fun) do { \
    int __i, __n = (fun)->N; \
    double __val = (fun)->v_exp; \
    for (__i = 0; __i < __n; ++__i) \
        __val -= DVAL((fun)->VAR[(fun)->V[__i]]); \
    (fun)->VAL = (fun)->weight * __val * __val; \
} while (0)

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

int adjust_sequence_copy_number_by_graph_layout(asg_t *asg, int *copy_number, int max_copy)
{
    uint32_t i, j, n_seg, n_group, a_g, *arc_group;
    uint64_t link_id;
    int updated;
    asmg_t *g;
    asmg_arc_t *a;

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

#ifdef DEBUG_SEG_COV_ADJUST
    fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] sequence copy number BEFORE adjusted by graph layout\n", __func__);
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] %s %u %u %d\n", __func__, 
                asg->seg[i].name, asg->seg[i].len, asg->seg[i].cov, copy_number[i]);
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
    asmg_arc_t *av;
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del)
            continue;
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
                kv_pushp(func_t, funcs, &FUN);
                FUN->weight = log10(asg->seg[i].len);
                FUN->v_exp = copy_number[i];
                FUN->VAR = VAR;
                FUN->N = V.n;
                FUN->V = V.a;
                FUN->VAL = 0.;
            }
        }
    }
 
    // do optimization
    int *arc_copy;
    int64_t sol_space_size;
    MYCALLOC(arc_copy, n_group);
    sol_space_size = 1;
    for (i = 0; i < n_group && sol_space_size <= BRUTE_FORCE_N_LIM; ++i)
        sol_space_size *= (arc_copy_ub[i] - arc_copy_lb[i] + 1);
    if (sol_space_size <= BRUTE_FORCE_N_LIM) {
        // do brute force optimization
        estimate_arc_copy_number_brute_force_impl(funcs.a, funcs.n, VAR, n_group, arc_copy, sol_space_size);
    } else {
        // do simulated annealing optimization
        estimate_arc_copy_number_siman_impl(funcs.a, funcs.n, VAR, n_group, arc_copy);
    }

#ifdef DEBUG_SEG_COV_ADJUST
    for (i = 0; i < n_group; ++i)
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] arc group %u optimum copy number %d [%u %u]\n",
                __func__, i, arc_copy[i], arc_copy_lb[i], arc_copy_ub[i]);
#endif

    // update sequence copy number using the arc copy number information
    int new_copy[2];
    updated = 0;
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del)
            continue;
        for (k = 0; k < 2; ++k) {
            v = i << 1 | k;
            na = asmg_arc_n(g, v);
            av = asmg_arc_a(g, v);
            new_copy[k] = 0;
            for (j = 0; j < na; ++j) {
                if (av[j].del)
                    continue;
                link_id = av[j].link_id;
                a_g = arc_group[link_id];
                assert(a_g != (uint32_t) -1);
                new_copy[k] += arc_copy[a_g];
            }
        }
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
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] %s %u %u %d\n", __func__, 
                asg->seg[i].name, asg->seg[i].len, asg->seg[i].cov, copy_number[i]);
    }
#endif

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
    free(arc_group);
    free(arc_copy_lb);
    free(arc_copy_ub);
    free(arc_copy);

    return updated;
}

kh_u32_t *sequence_duplication_by_copy_number(asg_t *asg, int *copy_number)
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
            fprintf(stderr, "[DEBUG_SEG_COPY::%s] make %d extra cop%s of %s [%u %u]\n", __func__, 
                    copy, copy > 1? "ies" : "y", asg->seg[i].name, asg->seg[i].len, asg->seg[i].cov);
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

void graph_path_finder(asg_t *asg, kh_u32_t *seg_dups, path_v *paths)
{
    uint64_t i, j, s, n_seg, n_leaf, len;
    int circ, exceed_limit;
    asmg_t *g;
    llnode *root, *next_node, *node, **leaves;
    kvec_t(llnodep) leaf_node, root_node;

    n_seg = asg->n_seg;
    g = asg->asmg;
    len = 0;
    s = UINT64_MAX;
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        if (g->vtx[i].len > len) {
            s = i;
            len = g->vtx[i].len;
        }
    }

    if (s == UINT64_MAX) return;

    kv_init(root_node);
    kv_init(leaf_node);

    root = new_node(s << 1);
    kv_push(llnodep, root_node, root);

    n_leaf = 0;
    exceed_limit = 0;
    leaves = graph_path_extension(g, root, seg_dups, &n_leaf, &exceed_limit);
    // for linear paths do extension from the other direction of root node
    // n_leaf will be zero if path number exceeds limit
    for (i = 0; i < n_leaf; ++i) {
        node = leaves[i];
        circ = asmg_arc1(g, node->v, s << 1) != 0;
        if (circ) {
            kv_push(llnodep, leaf_node, node);
            continue;
        }
        
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
        uint64_t l, l1;
        double cov, wl;

        l = g->vtx[path.a[0]>>1].len;
        cov = g->vtx[path.a[0]>>1].cov;
        wl = cov * l;
        for (j = 1; j < path.n; ++j) {
            a = asmg_arc1(g, path.a[j-1], path.a[j]);
            assert(!!a);
            l1 = g->vtx[path.a[j]>>1].len - a->lo;
            cov = g->vtx[path.a[j]>>1].cov;
            l += l1;
            wl += cov * l1;
        }
        if (circ) {
            a = asmg_arc1(g, path.a[path.n-1], path.a[0]);
            assert(!!a);
            l1 = a->lo;
            cov = g->vtx[path.a[0]>>1].cov;
            l -= l1;
            wl -= cov * l1;
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
    strncpy(dst, src, n);
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
    if (circ) len -= a->lo, wlen -= cov * a->lo;
    for (i = 1; i < vt.n; ++i) {
        len1 = g->vtx[vt.a[i]>>1].len;
        cov = g->vtx[vt.a[i]>>1].cov;
        len += len1;
        wlen += (double) cov * len1;
        a = asmg_arc1(g, vt.a[i-1], vt.a[i]);
        if (!a) {
            fprintf(stderr, "[E::%s] link does not exist: %s%c -> %s%c\n", __func__,
                    asg->seg[vt.a[i-1]>>1].name, "+-"[vt.a[i-1]&1],
                    asg->seg[vt.a[i]>>1].name, "+-"[vt.a[i]&1]);
            exit(EXIT_FAILURE);
        }
        len -= a->lo;
        wlen -= (double) cov * a->lo;
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

static double path_rotate_core(asg_t *g, path_t *path, hmm_annot_v *annots, uint8_t og_type)
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
        if ((annot->mito && og_type == OG_MITO) || (annot->pltd && og_type == OG_PLTD)) {
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

void path_rotate(asg_t *g, path_t *path, hmm_annot_v *annots, uint8_t og_type)
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

char comp_table[] = {
      0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
     16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
     32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
     48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
     64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
     64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
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

void print_seq(asg_t *asg, path_t *path, FILE *fo, int id, int force_linear, int line_wd)
{
    uint32_t i, n, v, lo, cov;
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
        lo = a->lo;
        cov = g->vtx[path->v[0]>>1].cov;
    }

    if (path->sid)
        fprintf(fo, ">%s\tlength=%u wlength=%.1f nv=%u circular=%s path=%s%c", path->sid, path->len + lo, path->wlen + (double) cov * lo,
                path->nv, (force_linear || !path->circ)? "false" : "true", asg->seg[path->v[0]>>1].name, "+-"[path->v[0]&1]);
    else
        fprintf(fo, ">%d\tlength=%u wlength=%.1f nv=%u circular=%s path=%s%c", id, path->len + lo, path->wlen + (double) cov * lo,
                path->nv, (force_linear || !path->circ)? "false" : "true", asg->seg[path->v[0]>>1].name, "+-"[path->v[0]&1]);

    for (i = 1; i < n; ++i)
        fprintf(fo, ",%s%c", asg->seg[path->v[i]>>1].name, "+-"[path->v[i]&1]);
    fprintf(fo, "\n");

    l = 0;
    v = path->v[0];
    put_chars(asg->seg[v>>1].seq, asg->seg[v>>1].len, v&1, (force_linear || !path->circ)? 0 : asmg_arc1(g, path->v[n-1], v)->lo, fo, &l, line_wd);

    for (i = 1; i < n; ++i) {
        v = path->v[i];
        a = asmg_arc1(g, path->v[i-1], v);
        assert(!!a);
        put_chars(asg->seg[v>>1].seq, asg->seg[v>>1].len, v&1, a->lo, fo, &l, line_wd);
    }

    if (!path->circ || !force_linear)
        assert(l == path->len);
    if (l % line_wd != 0)
        fputc('\n', fo);
}

uint32_t select_best_seq(asg_t *g, path_v *paths, FILE *fo, int type, double seq_cf, int seq_id, int cmp_coeff)
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
        if (k != UINT32_MAX && (double) l / paths->a[j].len >= seq_cf)
            j = k;
    }

    // TODO need more sophisticated rules for cmp_coeff
    if (cmp_coeff) {
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
    
    if (fo) print_seq(g, &paths->a[j], fo, seq_id > 0? seq_id : 1, 0, 60);

    return j;
}

void print_all_best_seqs(asg_t *g, path_v *paths, FILE *fo)
{
    if (paths->n == 0)
        return;

    size_t i;
    for (i = 0; i < paths->n; ++i)
        if (!!paths->a[i].best)
            print_seq(g, &paths->a[i], fo, i, 0, 60);
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
    uint64_t i;
    for (i = 0; i < g->n_seg; ++i)
        asg_seg_destroy(&g->seg[i]);
    // key has been freed by asg_seg_destroy
    if (g->h_seg) kh_sdict_destroy(g->h_seg);
    if (g->asmg) asmg_destroy(g->asmg);
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
    asmg_t *g, *g1;
    g = asg->asmg;
    g1 = asg1->asmg;
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

    return asg1;
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
    else if (*aux == 'i') val = aux_decimal_val(int32_t, aux + 1);
    else if (*aux == 'I') val = aux_decimal_val(uint32_t, aux + 1);
    else if (*aux == 'f') val = aux_decimal_val(float, aux + 1);
    return val;
}

static inline int gfa_aux_type2size(int x)
{
    if (x == 'C' || x == 'c' || x == 'A') return 1;
    else if (x == 'S' || x == 's') return 2;
    else if (x == 'I' || x == 'i' || x == 'f') return 4;
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
                    int32_t x;
                    x = strtol(q, &q, 10);
                    kputc_(type, &str); kputsn_((char*)&x, 4, &str);
                } else if (type == 'f') {
                    float x;
                    x = strtod(q, &q);
                    kputc_('f', &str); kputsn_(&x, 4, &str);
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
                        if (type == 'c')      while (q + 1 < p) { int8_t   x = strtol(q + 1, &q, 0); kputc_(x, &str); }
                        else if (type == 'C') while (q + 1 < p) { uint8_t  x = strtol(q + 1, &q, 0); kputc_(x, &str); }
                        else if (type == 's') while (q + 1 < p) { int16_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
                        else if (type == 'S') while (q + 1 < p) { uint16_t x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
                        else if (type == 'i') while (q + 1 < p) { int32_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
                        else if (type == 'I') while (q + 1 < p) { uint32_t x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
                        else if (type == 'f') while (q + 1 < p) { float    x = strtod(q + 1, &q);    kputsn_(&x, 4, &str); }
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
            LN = *(int32_t*)(s_LN + 1);
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
                    s->cov = len>0? dv/len : dv;
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
                    dv = *(int32_t*)(s_SBP_COV + 1);
                } else {
                    s_SBP_COV = gfa_aux_get(l_aux, aux, "FC");
                    if (s_SBP_COV && *s_SBP_COV == 'i')
                        dv = *(int32_t*)(s_SBP_COV + 1);
                }
                s->cov = len > 0? dv/len : dv;
            }
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
        arc = asmg_arc_add(g->asmg, v, w, ov, UINT64_MAX, 0, 0);
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
                    arc->cov = *(int32_t*)(s_ARC_COV + 1);
            }        
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
    uint32_t i;
    uint64_t k;
    asmg_t *asmg = g->asmg;

    fprintf(fo, "H\tVN:Z:1.0\n");
    for (i = 0; i < g->n_seg; ++i) {
        const asg_seg_t *s = &g->seg[i];
        // seg and vtx indices are interchangeable
        if (asmg && asmg->vtx[i].del) continue;
        fprintf(fo, "S\t%s\t", s->name);
        if (s->seq && !no_seq) fputs(s->seq, fo);
        else fputc('*', fo);
        fprintf(fo, "\tLN:i:%u\tKC:i:%lu\tSC:f:%.3f\n", s->len, (uint64_t)s->len*s->cov, (double)s->cov);
    }

    if (asmg == 0) return;
    for (k = 0; k < asmg->n_arc; ++k) {
        const asmg_arc_t *a = &asmg->arc[k];
        if (a->del || a->comp) continue;
        fprintf(fo, "L\t%s\t%c\t%s\t%c\t%luM\tEC:i:%u\n", g->seg[a->v>>1].name, "+-"[a->v&1], 
                g->seg[a->w>>1].name, "+-"[a->w&1], a->lo, a->cov);
    }
}

void asg_print_fa(asg_t *g, FILE *fo, int line_wd)
{
    uint64_t i, l = 0;
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

