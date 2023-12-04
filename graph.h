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
 * 03/11/22 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/

#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <stdlib.h>
#include <stdint.h>

#include "misc.h"

typedef struct {
    uint64_t v; // node id << 1 | rev
    uint64_t w; // node id << 1 | rev
    uint64_t ln; // v->w overlapping syncmer number
    uint64_t ls; // v->w overlapping consensus sequence length
    uint32_t cov:30, del:1, comp:1; // arc coverage, deleted, reverse complement
    uint64_t link_id; // arc and comp share the same link_id
} asmg_arc_t;

typedef struct {
    uint64_t n; // number syncmer
    uint64_t *a; // syncmer list
    char *seq; // concensus sequence
    uint64_t len; // sequence length
    uint32_t cov:30, del:1, circ:1;
} asmg_vtx_t;

typedef struct {
    uint64_t n_vtx, m_vtx; // number nodes
    asmg_vtx_t *vtx;
    uint64_t n_arc, m_arc;
    asmg_arc_t *arc;
    uint64_t *idx_p; // node starting positions
    uint64_t *idx_n; // node number arcs
} asmg_t;

#define asmg_arc_head(a) ((a).v)
#define asmg_arc_tail(a) ((a).w)
#define asmg_arc_n(g, v) ((g)->idx_n[(v)])
#define asmg_arc_a(g, v) (&(g)->arc[(g)->idx_p[(v)]])
#define asmg_arc_id(a) ((a).link_id << 1 | (a).comp)
#define asmg_vtx_n(g) ((g)->n_vtx << 1)

#ifdef __cplusplus
extern "C" {
#endif

void asmg_destroy(asmg_t *g);
void asmg_arc_index(asmg_t *g);
void asmg_arc_sort(asmg_t *g);
int asmg_arc_is_sorted(asmg_t *g);
uint64_t asmg_max_link_id(asmg_t *g);
void asmg_shrink_link_id(asmg_t *g);
void asmg_finalize(asmg_t *g, int do_cleanup);
void asmg_arc_fix_cov(asmg_t *g);
uint64_t *asmg_vtx_list(asmg_t *g, uint64_t *_n);
void asmg_print(asmg_t *g, FILE *fo, int no_seq);
uint32_t *asmg_uext_arc_group(asmg_t *g, uint32_t *n);
asmg_t *asmg_unitigging(asmg_t *g);
uint64_t asmg_drop_tip(asmg_t *g, int32_t tip_cnt, uint64_t tip_len, int protect_super_tip, int do_cleanup, int VERBOSE);
uint64_t asmg_pop_bubble(asmg_t *g, uint64_t radius, uint64_t max_del, int protect_tip, int protect_super_bubble, int do_cleanup, int VERBOSE);
uint64_t asmg_remove_weak_crosslink(asmg_t *g, double c_thresh, double m_cov, int do_cleanup, int VERBOSE);
uint32_t *asmg_subgraph(asmg_t *g, uint32_t *seeds, uint32_t n, uint32_t step, uint64_t dist, uint32_t *_nv, int modify_graph);
int asmg_tarjans_scc(asmg_t *g, int *scc);
int asmg_path_exists(asmg_t *g, uint32_t source, uint32_t sink, uint32_t step, uint64_t dist, uint32_t *_step, uint64_t *_dist);
#ifdef __cplusplus
}
#endif

static inline void asmg_arc_del(asmg_t *g, uint64_t v, uint64_t w, uint32_t del)
{
    uint64_t i, nv = asmg_arc_n(g, v);
    asmg_arc_t *av = asmg_arc_a(g, v);
    for (i = 0; i < nv; ++i)
        if (av[i].w == w)
            av[i].del = del;
}

static inline void asmg_arc_del_v(asmg_t *g, uint64_t v, uint32_t del)
{
    uint64_t i, nv = asmg_arc_n(g, v);
    asmg_arc_t *av = asmg_arc_a(g, v);
    for (i = 0; i < nv; ++i) {
        av[i].del = del;
        asmg_arc_del(g, av[i].w^1, v^1, del);
    }
}

static inline void asmg_vtx_del(asmg_t *g, uint64_t s, uint32_t del)
{
    g->vtx[s].del = del;
    asmg_arc_del_v(g, s<<1, del);
    asmg_arc_del_v(g, s<<1|1, del);
}

static inline void asmg_vtx_add(asmg_t *g)
{
    if (g->n_vtx == g->m_vtx)
        MYEXPAND(g->vtx, g->m_vtx);
    g->n_vtx++;
}

static inline void asmg_vtx_addp(asmg_t *g, asmg_vtx_t **p)
{
    if (g->n_vtx == g->m_vtx)
        MYEXPAND(g->vtx, g->m_vtx);
    *p = &g->vtx[g->n_vtx++];
}

static inline uint64_t asmg_vtx_n1(asmg_t *g)
{
    uint64_t i, n, n1;
    
    n1 = g->n_vtx;
    for (i = 0, n = g->n_vtx; i < n; ++i)
        if (g->vtx[i].del)
            --n1;

    return n1;
}

static inline uint64_t asmg_arc_n1(asmg_t *g, uint64_t v)
{
    uint64_t i, n, n1;
    asmg_arc_t *a;

    a = asmg_arc_a(g, v);
    n = asmg_arc_n(g, v);
    n1 = n;
    for (i = 0; i < n; ++i)
        if (a[i].del)
            --n1;

    return n1;
}

static inline asmg_arc_t *asmg_arc_a1(asmg_t *g, uint64_t v)
{
    uint64_t i, n;
    asmg_arc_t *a;

    a = asmg_arc_a(g, v);
    n = asmg_arc_n(g, v);
    for (i = 0; i < n; ++i)
        if (!a[i].del)
            return &a[i];

    return 0;
}

static inline asmg_arc_t *asmg_arc(asmg_t *g, uint64_t v, uint64_t w)
{
    uint64_t i, n;
    asmg_arc_t *a;

    a = asmg_arc_a(g, v);
    n = asmg_arc_n(g, v);
    for (i = 0; i < n; ++i)
        if (a[i].w == w)
            return &a[i];

    return 0;
}

static inline asmg_arc_t *asmg_arc1(asmg_t *g, uint64_t v, uint64_t w)
{
    uint64_t i, n;
    asmg_arc_t *a;

    a = asmg_arc_a(g, v);
    n = asmg_arc_n(g, v);
    for (i = 0; i < n; ++i)
        if (a[i].w == w && !a[i].del)
            return &a[i];

    return 0;
}

static inline asmg_arc_t *asmg_comp_arc(asmg_t *g, asmg_arc_t *a)
{
    return asmg_arc(g, a->w^1, a->v^1);
}

static inline asmg_arc_t *asmg_comp_arc1(asmg_t *g, asmg_arc_t *a)
{
    return asmg_arc1(g, a->w^1, a->v^1);
}

static inline int asmg_arc_exist(asmg_t *g, uint64_t v, uint64_t w)
{
    uint64_t i, n;
    asmg_arc_t *a;

    a = asmg_arc_a(g, v);
    n = asmg_arc_n(g, v);
    for (i = 0; i < n; ++i)
        if (a[i].w == w)
            return 1;

    return 0;
}

static inline int asmg_arc_exist1(asmg_t *g, uint64_t v, uint64_t w)
{
    uint64_t i, n;
    asmg_arc_t *a;

    a = asmg_arc_a(g, v);
    n = asmg_arc_n(g, v);
    for (i = 0; i < n; ++i)
        if (a[i].w == w && !a[i].del)
            return 1;

    return 0;
}

static inline asmg_arc_t *asmg_arc_add(asmg_t *g, uint64_t v, uint64_t w, uint64_t ln, uint64_t ls, uint64_t link_id, uint32_t cov, uint32_t comp)
{
    asmg_arc_t *a;
    if (g->n_arc == g->m_arc)
        MYEXPAND(g->arc, g->m_arc);
    a = &g->arc[g->n_arc++];
    a->v = v;
    a->w = w;
    a->ln = ln;
    a->ls = ls;
    a->link_id = link_id;
    a->cov = cov;
    a->del = 0;
    a->comp = comp;
    return a;
}

static inline void asmg_arc_add2(asmg_t *g, uint64_t v, uint64_t w, uint64_t ln, uint64_t ls, uint64_t link_id, uint32_t cov, uint32_t comp)
{
    // add both arc and its comp
    asmg_arc_add(g, v, w, ln, ls, link_id, cov, comp);
    if (v != (w^1) || w != (v^1))
        asmg_arc_add(g, w^1, v^1, ln, ls, link_id, cov, comp^1);
}

static inline uint64_t asmg_comp_arc_id(asmg_arc_t *a)
{
    return ((a->v^1) != a->w || (a->w^1) != a->v)?
        (asmg_arc_id(*a)^1) : asmg_arc_id(*a);
}

/***
static inline uint64_t asmg_comp_arc_id1(uint64_t v, uint64_t w, uint64_t arc_id)
{
    return ((v^1) != w || (w^1) != v)? (arc_id^1) : arc_id;
}
**/

static inline void asmg_clean_consensus(asmg_t *g)
{
    uint64_t i, n;
    for (i = 0, n = g->n_arc; i < n; ++i)
        g->arc[i].ls = 0;
    for (i = 0, n = g->n_vtx; i < n; ++i) {
        if (g->vtx[i].seq) {
            free(g->vtx[i].seq);
            g->vtx[i].seq = 0;
        }
        g->vtx[i].len = 0;
    }
}

static inline uint64_t asmg_arc_head_e(asmg_t *g, asmg_arc_t *a)
{
    uint64_t v = a->v;
    asmg_vtx_t *vtx = &g->vtx[v>>1];
    return (v&1)? (vtx->a[0]^1) : (vtx->a[vtx->n-1]);
}

static inline uint64_t asmg_arc_tail_e(asmg_t *g, asmg_arc_t *a)
{
    uint64_t w = a->w;
    asmg_vtx_t *vtx = &g->vtx[w>>1];
    return (w&1)? (vtx->a[vtx->n-1]^1) : (vtx->a[0]);
}

#endif


