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
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "kvec.h"

#include "misc.h"
#include "graph.h"

static void asmg_vtx_destroy(asmg_vtx_t *v)
{
    if (v->seq) free(v->seq);
    if (v->a) free(v->a);
}

void asmg_destroy(asmg_t *g)
{
    uint64_t i;
    for (i = 0; i < g->n_vtx; ++i)
        asmg_vtx_destroy(&g->vtx[i]);
    if (g->vtx) free(g->vtx);
    if (g->arc) free(g->arc);
    if (g->idx_p) free(g->idx_p);
    if (g->idx_n) free(g->idx_n);
    free(g);
}

int asmg_arc_is_sorted(asmg_t *g)
{
    uint64_t e;
    for (e = 1; e < g->n_arc; ++e)
        if (g->arc[e-1].v > g->arc[e].v ||
                (g->arc[e-1].v == g->arc[e].v && g->arc[e-1].w > g->arc[e].w))
            break;
    return (e == g->n_arc);
}

static int asmg_arc_cmpfunc(const void *a, const void *b)
{
    uint64_t x, y, s, t;
    x = ((asmg_arc_t *) a)->v;
    y = ((asmg_arc_t *) b)->v;
    s = ((asmg_arc_t *) a)->w;
    t = ((asmg_arc_t *) b)->w;
    return x == y? ((s > t) - (s < t)) : ((x > y) - (x < y));
}

void asmg_arc_sort(asmg_t *g)
{
    qsort(g->arc, g->n_arc, sizeof(asmg_arc_t), asmg_arc_cmpfunc);
}

void asmg_arc_index(asmg_t *g)
{
    uint64_t i, last, n, v, *idx_p, *idx_n;
    asmg_arc_t *a;

    if (g->idx_p) {
        free(g->idx_p);
        free(g->idx_n);
    }

    MYCALLOC(idx_p, g->n_vtx * 2);
    MYCALLOC(idx_n, g->n_vtx * 2);
    g->idx_p = idx_p;
    g->idx_n = idx_n;

    if (!g->n_arc) return;

    a = g->arc;
    n = g->n_arc;
    v = asmg_arc_head(a[0]);
    for (i = 1, last = 0; i < n; ++i) {
        if (v != asmg_arc_head(a[i])) {
            idx_p[v] = last;
            idx_n[v] = i - last;
            last = i;
            v = asmg_arc_head(a[i]);
        }
    }
    idx_p[v] = last;
    idx_n[v] = i - last;
}

uint64_t asmg_max_link_id(asmg_t *g)
{
    uint64_t i, n, link_id;

    link_id = 0;
    for (i = 0, n = g->n_arc; i < n; ++i)
        link_id = MAX(link_id, g->arc[i].link_id);

    return link_id;
}

void asmg_shrink_link_id(asmg_t *g)
{
    uint64_t i, n, link_id;
    asmg_arc_t *a;
    
    for (i = 0, n = g->n_arc; i < n; ++i)
        g->arc[i].link_id |= 0x8000000000000000;

    link_id = 0;
    for (i = 0, n = g->n_arc; i < n; ++i) {
        a = &g->arc[i];
        if (a->link_id & 0x8000000000000000) {
            a->link_id = link_id;
            a = asmg_arc(g, a->w^1, a->v^1);
            a->link_id = link_id;
            ++link_id;
        }
    }
}

static void asmg_cleanup(asmg_t *g) {
    uint64_t i, j, n, *v_idx;
    asmg_arc_t *a;

    // delete vtx
    n = g->n_vtx;
    MYMALLOC(v_idx, n);
    MYBONE(v_idx, n);
    for (i = 0; i < n; ++i)
        if (g->vtx[i].del)
            asmg_vtx_destroy(&g->vtx[i]);
    for (i = j = 0; i < n; ++i) {
        if (g->vtx[i].del) {
            // update vtx number
            --g->n_vtx;
            continue;
        }
        if (j < i)
            // copy vtx_i to vtx_j position
            g->vtx[j] = g->vtx[i];
        // recode new vtx index
        v_idx[i] = j++;
    }
    MYREALLOC(g->vtx, g->n_vtx);
    g->m_vtx = g->n_vtx;

    // delete arc
    n = g->n_arc;
    for (i = j = 0; i < n; ++i) {
        a = &g->arc[i];
        if (g->arc[i].del ||
                v_idx[a->v>>1] == UINT64_MAX ||
                v_idx[a->w>>1] == UINT64_MAX) {
            // update arc number
            --g->n_arc;
            continue;
        }
        if (j < i)
            // copy arc_i to arc_j position
            g->arc[j] = g->arc[i];
        ++j;
    }
    MYREALLOC(g->arc, g->n_arc);
    g->m_arc = g->n_arc;

    // update arc vtx
    n = g->n_arc;
    for (i = 0; i < n; ++i) {
        a = &g->arc[i];
        a->v = v_idx[a->v>>1] << 1 | (a->v&1);
        a->w = v_idx[a->w>>1] << 1 | (a->w&1);
    }

    free(v_idx);
}

static uint32_t asmg_arc_fix_symm(asmg_t *g) {
    uint32_t symm_fix;
    uint64_t i, v, w, n;
    asmg_arc_t *a, *a1;

    symm_fix = 0;
    for (i = 0, n = g->n_arc; i < n; ++i) {
        a = &g->arc[i];
        if (a->del)
            continue;
        v = a->v;
        w = a->w;
        a1 = asmg_arc1(g, w^1, v^1);
        if (a1 == 0) {
            // add symm arc
            a1 = asmg_arc_add(g, w^1, v^1, a->lo, a->link_id, a->cov, a->comp^1);
            ++symm_fix;
        } else {
            // fix symm flag
            a1->comp = a->comp^1;
        }
        // FIXME
        if (a->lo != a1->lo)
            a->lo = a1->lo = MIN(a->lo, a1->lo);
    }

    return symm_fix;
}

void asmg_finalize(asmg_t *g, int do_cleanup) {
    // TODO fix multi-arc
    uint32_t symm_fix;
    if (do_cleanup) asmg_cleanup(g);
    asmg_arc_sort(g);
    asmg_arc_index(g);
    symm_fix = asmg_arc_fix_symm(g);
    if (symm_fix > 0) {
        asmg_arc_sort(g);
        asmg_arc_index(g);
    }
    asmg_shrink_link_id(g);
}

typedef struct { size_t n, m; uint64_t *a; } u64_v_t;

/*********************
 * Probe unitig ends *
 *********************/

#define ASMG_VT_MERGEABLE 0
#define ASMG_VT_TIP       1
#define ASMG_VT_MULTI_OUT 2
#define ASMG_VT_MULTI_NEI 3
/*
 * ASMG_VT_TIP       : no outgoing arc
 * ASMG_VT_MERGEABLE : single outgoing arc v->w, w single incoming arc
 * ASMG_VT_MULTI_NEI : single outgoing arc v->w, w multiple incoming arc
 * ASMG_VT_MULTI_OUT : multiple outgoing arcs
 */
static inline uint64_t asmg_arc_n2(asmg_t *g, uint64_t v, uint64_t *w, uint64_t *l)
{
    uint64_t i, nv, nv0, k, lo, min_l;
    asmg_arc_t *av;
    *l = 0, *w = UINT64_MAX;
    if (g->vtx[v>>1].del)
        return 0;
    nv0 = asmg_arc_n(g, v);
    av = asmg_arc_a(g, v);
    lo = k = 0;
    for (i = nv = 0; i < nv0; ++i) {
        if (!av[i].del) {
            ++nv;
            k = i;
            lo = av[i].lo > lo? av[i].lo : lo;
        }
    }
    min_l = g->vtx[v>>1].len - lo;
    *l = min_l;
    *w = nv == 1? av[k].w : UINT64_MAX;
    return nv;
}

static inline int32_t asmg_uext(asmg_t *g, uint64_t v, int32_t max_ext, uint64_t *ne, uint64_t *le, u64_v_t *a)
{
    // the extension of circular unitigs will be bounded by max_ext
    // and the return value will be ASMG_VT_MERGEABLE
    int32_t vt;
    uint64_t nv, nw, l, w, n_ext, l_ext;
    n_ext = l_ext = 0;
    if (a) a->n = 0;
    if (a) kv_push(uint64_t, *a, v);
    do {
        nv = asmg_arc_n2(g, v, &w, &l);
        if (nv == 0) {
            vt = ASMG_VT_TIP;
        } else if (nv > 1) {
            vt = ASMG_VT_MULTI_OUT;
        } else {
            nw = asmg_arc_n1(g, w^1);
            vt = nw == 1? ASMG_VT_MERGEABLE : ASMG_VT_MULTI_NEI;
        }
        l_ext += l;
        if (vt != ASMG_VT_MERGEABLE) break;
        ++n_ext;
        if (a) kv_push(uint64_t, *a, w);
        v = w;
    } while (--max_ext > 0);
    if (ne) *ne = n_ext;
    if (le) *le = l_ext;
    return vt;
}

uint32_t *asmg_uext_arc_group(asmg_t *g, uint32_t *n_group)
{
    uint64_t i, j, k, v, n_vtx, n_arc;
    uint32_t *arc_group, group, na;
    uint8_t *visited;
    int32_t vt;
    u64_v_t a = {0, 0, 0};

    n_vtx = g->n_vtx;
    n_arc = asmg_max_link_id(g); // number unique arcs
    
    MYMALLOC(arc_group, n_arc);
    MYBONE(arc_group, n_arc);
    MYCALLOC(visited, n_vtx);

    group = 0;
    for (i = 0; i < n_vtx; ++i) {
        if (visited[i] || g->vtx[i].del)
            continue;
        na = 0;
        for (k = 0; k < 2; ++k) {
            v = i << 1 | k;
            vt = asmg_uext(g, v, n_vtx*2+1, 0, 0, &a);
            for (j = 1; j < a.n; ++j) {
                ++na;
                arc_group[asmg_arc1(g, a.a[j-1], a.a[j])->link_id] = group;
                visited[a.a[j]>>1] = 1;
            }
            if (vt == ASMG_VT_MULTI_NEI) {
                arc_group[asmg_arc_a1(g, a.a[a.n-1])->link_id] = group;
                ++na;
            }
        }
        if (na > 0) ++group;
        visited[i] = 1;
    }
    
    for (i = 0; i < g->n_arc; ++i) {
        if (g->arc[i].del || 
                arc_group[g->arc[i].link_id] != (uint32_t) -1)
            continue;
        arc_group[g->arc[i].link_id] = group++;
    }

    if (n_group) *n_group = group;

    free(a.a);
    free(visited);

    return arc_group;
}

/*************************
 * Topological extension *
 *************************/

typedef struct {
    uint64_t p; // the optimal parent vertex
    uint64_t d; // the shortest distance from the initial vertex
    uint64_t c; // coverage weight path length
    uint64_t r:63, s:1; // r: the number of remaining incoming arc; s: whether the vertex has been visited before
} asmg_tinfo_t;

typedef struct {
    asmg_tinfo_t *a;
    kvec_t(uint64_t) S; // set of vertices without parents
    kvec_t(uint64_t) b; // visited vertices
    kvec_t(uint64_t) e; // visited edges/arcs
    // n_short_tip: number of branching short tips
    // n_sink: number of sinks
    // dist: min distance from the input vertex
    // v_sink: end vertex; dist and v_sink only set when n_sink>0
    uint64_t n_short_tip, n_sink, dist, v_sink;
    uint8_t self_cycle;
} asmg_tbuf_t;

#define ASMG_TE_THRU_SHORT_TIP 0x1
#define ASMG_TE_THRU_BUBBLE    0x2

static asmg_tbuf_t *asmg_tbuf_init(asmg_t *g)
{
    uint64_t v, n_vtx;
    n_vtx = asmg_vtx_n(g);
    asmg_tbuf_t *b;
    MYCALLOC(b, 1);
    MYCALLOC(b->a, n_vtx);
    for (v = 0; v < n_vtx; ++v)
        b->a[v].p = UINT64_MAX;
    return b;
}

static void asmg_tbuf_destroy(asmg_tbuf_t *b)
{
    if (b == 0)
        return;
    free(b->a);
    free(b->S.a);
    free(b->b.a);
    free(b->e.a);
    free(b);
}

static void asmg_tbuf_reset(asmg_tbuf_t *b)
{
    uint64_t i;
    for (i = 0; i < b->b.n; ++i) { // clear the states of visited vertices
        asmg_tinfo_t *t = &b->a[b->b.a[i]];
        t->d = t->c = t->r = t->s = 0;
        t->p = UINT64_MAX;
    }
}

static uint64_t asmg_topo_ext(asmg_t *g, uint64_t v0, uint64_t max_dist, int32_t thru_flag, asmg_tbuf_t *b)
{
    if (g->vtx[v0>>1].del) return 0;

    uint64_t i, v, w, l, nv, d, c, a, n_pending, max_d;
    asmg_tinfo_t *t;
    asmg_arc_t *av;

    // n_pending: number of visited vertices that are not sorted
    // max_d: max asmg_tinfo_t::d of visited vertices
    n_pending = max_d = 0;
    b->S.n = b->b.n = b->e.n = 0;
    b->n_short_tip = b->n_sink = b->dist = 0;
    b->self_cycle = 0;
    b->v_sink = UINT64_MAX;
    t = &b->a[v0];
    t->d = t->c = t->r = t->s = 0;
    t->p = UINT64_MAX;

    kv_push(uint64_t, b->S, v0);

    while (b->S.n > 0 && max_d <= max_dist) {
        v = kv_pop(b->S);
        nv = asmg_arc_n(g, v);
        av = asmg_arc_a(g, v);
        d = b->a[v].d;
        c = b->a[v].c;
        if (b->S.n == 0 && n_pending == 0) { // a sink vertex
            b->dist = d;
            b->v_sink = v;
            if (v != v0) { // exclude the input vertex
                ++b->n_sink;
                if (!(thru_flag & ASMG_TE_THRU_BUBBLE))
                    break;
            }
        }

        if (asmg_arc_n1(g, v) == 0) { // a tip
            if (d + g->vtx[v>>1].len < max_dist) {
                // a tip shorter than max_dist
                if (b->S.n || n_pending)
                    // don't count a tip if it is the end of a bubble chain
                    ++b->n_short_tip;
                if (thru_flag & ASMG_TE_THRU_SHORT_TIP) continue;
                else break;
            } else break; // if we come here, we have a tip beyond max_dist; we stop
        }

        for (i = 0; i < nv; ++i) { // loop through v's neighbors
            if (av[i].del)
                continue;
            w = av[i].w;
            l = g->vtx[v>>1].len - av[i].lo;
            a = g->vtx[v>>1].cov * l;
            t = &b->a[w];
            if (w>>1 == v0>>1) {
                // cycle or bidirected cycle involving the input vertex
                b->self_cycle |= w == v0? 1 : 2;
                break;
            }
            kv_push(uint64_t, b->e, g->idx_p[v] + i); // TODO save the edge
            if (t->s == 0) { // this vertex has never been visited
                kv_push(uint64_t, b->b, w); // save it for asmg_tbuf_reset()
                t->p = v, t->s = 1, t->d = d + l, t->c = c + a;
                t->r = asmg_arc_n1(g, w^1);
                ++n_pending;
            } else { // visited before
                if (c + a > t->c || (c + a == t->c && d + l > t->d)) t->p = v;
                if (c + a > t->c) t->c = c + a;
                if (d + l < t->d) t->d = d + l; // update dist
            }
            max_d = max_d > t->d? max_d : t->d;
            assert(t->r > 0);
            assert(n_pending > 0);
            if (--(t->r) == 0) {
                kv_push(uint64_t, b->S, w);
                --n_pending;
            }
        }
        if (i < nv) break;
    }
    return b->n_sink;
}

/******************
 * graph cleaning *
 * ****************/

uint64_t asmg_drop_tip(asmg_t *g, int32_t tip_cnt, uint64_t tip_len, int do_cleanup, int VERBOSE)
{
    uint64_t i, v, n_vtx, cnt, l_ext;
    int32_t vt;
    u64_v_t a = {0, 0, 0};
    n_vtx = asmg_vtx_n(g);
    cnt = 0;
    for (v = 0; v < n_vtx; ++v) {
        if (g->vtx[v>>1].del) continue;
        if (asmg_arc_n1(g, v^1) != 0) continue; // not a tip
        vt = asmg_uext(g, v, tip_cnt, 0, &l_ext, &a);
        if (vt == ASMG_VT_MERGEABLE) continue; // circular unitig
        if (l_ext > tip_len) continue; // tip too long
        for (i = 0; i < a.n; ++i)
            asmg_vtx_del(g, a.a[i]>>1, 1);
        ++cnt;
    }
    free(a.a);
    if (do_cleanup && cnt > 0) asmg_finalize(g, 1);
    if (VERBOSE)
        fprintf(stderr, "[M::%s] dropped %lu tips\n", __func__, cnt);
    return cnt;
}

static inline uint64_t asmg_arc_n3(asmg_t *g, uint64_t v, uint64_t w, uint32_t *c)
{
    uint64_t i, nv, nv0, k;
    asmg_arc_t *av;
    nv0 = asmg_arc_n(g, v);
    av = asmg_arc_a(g, v);
    for (i = k = nv = 0; i < nv0; ++i) {
        if (!av[i].del) {
            ++nv;
            if (av[i].w != w)
                k = i;
        }
    }
    *c = nv == 2? av[k].cov : 0;
    return nv;
}

uint64_t asmg_remove_weak_crosslink(asmg_t *g, double c_thresh, int do_cleanup, int VERBOSE)
{
    uint64_t i, n, v, w, cnt = 0;
    uint32_t cv, cw;
    asmg_arc_t *a;

    for (i = 0, n = g->n_arc; i < n; ++i) {
        a = &g->arc[i];
        if (a->del) continue;
        v = a->v;
        w = a->w;
        if (asmg_arc_n1(g, v^1) != 1 || asmg_arc_n1(g, w) != 1)
            continue;
        if (asmg_arc_n3(g, v, w, &cv) != 2 || 
                asmg_arc_n3(g, w^1, v^1, &cw) != 2)
            continue;
        cv = MIN(cv, cw);
        if (cv > 0 && (double) a->cov / cv < c_thresh) {
            a->del = 1;
            asmg_arc_del(g, w^1, v^1, 1);
            ++cnt;
        }
    }
    if (do_cleanup && cnt > 0) asmg_finalize(g, 1);
    if (VERBOSE)    
        fprintf(stderr, "[M::%s] dropped %lu weak cross links\n", __func__, cnt);
    return cnt;
}

/******************
 * Bubble popping *
 ******************/

// in a resolved bubble, mark unused vertices and arcs as "reduced"
static void asmg_bub_backtrack(asmg_t *g, uint64_t v0, uint64_t max_del, int protect_super_bubble, asmg_tbuf_t *b)
{
    uint64_t i, v, w;
    asmg_arc_t *a;

    assert(b->S.n == 0);
    if (max_del > 0) {
        uint64_t n_kept = 0;
        v = b->v_sink;
        do { ++n_kept, v = b->a[v].p; } while (v != v0);
        if (b->b.n > n_kept + max_del) return;
    }
    if (protect_super_bubble) {
        uint64_t n_kept, b_kept, b_tot;
        n_kept = b_kept = 0;
        v = b->v_sink;
        do { ++n_kept, b_kept += g->vtx[v>>1].len, v = b->a[v].p; } while (v != v0);
        b_tot = 0;
        for (i = 0; i < b->b.n; ++i) b_tot += g->vtx[b->b.a[i]>>1].len;
        if ((b_tot - b_kept) * 2 > (g->vtx[v0>>1].len + g->vtx[b->v_sink>>1].len) * (b->b.n - n_kept)) return;
    }
    for (i = 0; i < b->b.n; ++i)
        g->vtx[b->b.a[i]>>1].del = 1;
    for (i = 0; i < b->e.n; ++i) {
        a = &g->arc[b->e.a[i]];
        a->del = 1;
        asmg_arc_del(g, a->w^1, a->v^1, 1);
    }
    v = b->v_sink;
    do {
        w = b->a[v].p; // w->v
        g->vtx[v>>1].del = 0;
        asmg_arc_del(g, w, v, 0);
        asmg_arc_del(g, v^1, w^1, 0);
        v = w;
    } while (v != v0);
}

// pop bubbles from vertex v0; the graph MUST BE symmetric: if u->v present, v'->u' must be present as well
// radius is calculated from the end of v0, not the start
static uint64_t asmg_bub_pop1(asmg_t *g, uint64_t v0, uint64_t radius, uint64_t max_del, int protect_tip, int protect_super_bubble, asmg_tbuf_t *b)
{
    uint64_t ret = 0;
    if (asmg_arc_n1(g, v0) < 2) return 0; // no bubbles
    asmg_topo_ext(g, v0, g->vtx[v0>>1].len + radius, protect_tip? 0 : ASMG_TE_THRU_SHORT_TIP, b);
    if (b->n_sink) {
        asmg_bub_backtrack(g, v0, max_del, protect_super_bubble, b);
        ret = 1 | (uint64_t) b->n_short_tip << 32;
    }
    asmg_tbuf_reset(b);
    return ret;
}

// pop bubbles
// super bubble: the average length of sequences to be removed is greater than the average length of source and sink sequence
// super buubles are protected with protect_super_bubble set
uint64_t asmg_pop_bubble(asmg_t *g, uint64_t radius, uint64_t max_del, int protect_tip, int protect_super_bubble, int do_cleanup, int VERBOSE)
{
    uint64_t v, n_vtx, n_pop;
    asmg_tbuf_t *b;

    n_vtx = asmg_vtx_n(g);
    b = asmg_tbuf_init(g);
    n_pop = 0;
    for (v = 0; v < n_vtx; ++v)
        if (!g->vtx[v>>1].del && asmg_arc_n1(g, v) >= 2)
            n_pop += asmg_bub_pop1(g, v, radius, max_del, protect_tip, protect_super_bubble, b);
    asmg_tbuf_destroy(b);
    if (do_cleanup && n_pop > 0) asmg_finalize(g, 1);
    if (VERBOSE)
        fprintf(stderr, "[M::%s] popped %u bubbles and trimmed %u short tips\n", __func__, (uint32_t) n_pop, (uint32_t) (n_pop>>32));
    return n_pop;
}

/****************
 * Unitig graph *
 ****************/

static inline int is_junction(asmg_t *g, uint64_t s)
{
    return asmg_arc_n1(g, s << 1) > 1 || asmg_arc_n1(g, s << 1 | 1) > 1;
}

static inline void vec_add(u64_v_t *v_vec, uint64_t *v, size_t n, int r)
{
    size_t i, j;
    if (r) {
        for (i = 0, j = n-1; i < n; ++i, --j)
            kv_push(uint64_t, *v_vec, v[j]^1);
    } else {
        for (i = 0; i < n; ++i)
            kv_push(uint64_t, *v_vec, v[i]);
    }
}

asmg_t *asmg_unitigging(asmg_t *g)
{
    uint64_t i, j, k, n_vtx, n_arc, n_arc1, v, w, n, m, *vtx_p;
    asmg_arc_t *arc, *arc1;
    asmg_vtx_t *vtx, *vtx1;
    uint8_t *vtx_visited;
    asmg_t *utg_asmg;
    kvec_t(asmg_vtx_t) utgs = {0, 0 ,0};
    kvec_t(asmg_arc_t) arcs = {0, 0, 0};
    vtx = g->vtx;
    n_vtx = g->n_vtx;
    MYCALLOC(vtx_visited, n_vtx);

    /* first visit for all untigs connected to a junction */
    for (i = 0; i < n_vtx; ++i) {
        if (vtx[i].del || !is_junction(g, i))
            continue;
        for (k = 0; k < 2; ++k) {
            v = i << 1 | k;
            n_arc = asmg_arc_n(g, v);
            n_arc1 = asmg_arc_n1(g, v);
            arc = asmg_arc_a(g, v);
            for (j = 0; j < n_arc; ++j) {
                if (arc[j].del)
                    continue;
                kvec_t(uint64_t) vec = {0, 0, 0};
                if (!vtx_visited[v>>1] && n_arc1 == 1)
                    kv_push(uint64_t, vec, v); /* add to utg vec */
                v = asmg_arc_tail(arc[j]);
                /* v can be the start point of a unitig */
                while (!vtx_visited[v>>1] && asmg_arc_n1(g, v^1) == 1) {
                    kv_push(uint64_t, vec, v);
                    vtx_visited[v>>1] = 1;
                    if (asmg_arc_n1(g, v) == 1) /* at most one outgoing edge */
                        v = asmg_arc_a1(g, v)->w;
                    else
                        break;
                }
                if (vec.n > 1) {
                    kv_pushp(asmg_vtx_t, utgs, &vtx1);
                    vtx1->n = vec.n;
                    vtx1->a = vec.a;
                    vtx1->seq = 0;
                    vtx1->len = 0;
                    vtx1->cov = vtx1->del = vtx1->circ = 0;
                    MYREALLOC(vtx1->a, vtx1->n);
                }
            }
        }
        vtx_visited[i] = 1;
    }

    /* second visit to get linear path */
    /* TODO */
    /* consider this special case as linear */
    /* v(-)->v(+) v(+)->w(+) w(+)->w(-) */
    for (i = 0; i < n_vtx; ++i) {
        if (vtx[i].del || vtx_visited[i] ||
                (asmg_arc_n1(g, i << 1) > 0 &&
                 asmg_arc_n1(g, i << 1 | 1) > 0))
            continue;
        v = asmg_arc_n1(g, i << 1) > 0? i << 1 : i << 1 | 1;
        kvec_t(uint64_t) vec = {0, 0, 0};
        kv_push(uint64_t, vec, v);
        vtx_visited[v>>1] = 1;
        while (asmg_arc_n1(g, v) == 1) { /* should always be 1 */
            v = asmg_arc_a1(g, v)->w;
            if (!vtx_visited[v>>1])
                kv_push(uint64_t, vec, v);
            else
                break;
            vtx_visited[v>>1] = 1;
        }
        if (vec.n > 1) {
            kv_pushp(asmg_vtx_t, utgs, &vtx1);
            vtx1->n = vec.n;
            vtx1->a = vec.a;
            vtx1->seq = 0;
            vtx1->len = 0;
            vtx1->cov = vtx1->del = vtx1->circ = 0;
            MYREALLOC(vtx1->a, vtx1->n);
        }
    }

    /* last visit for the remaining unvisited non-juntions on circles */
    for (i = 0; i < n_vtx; ++i) {
        if (vtx[i].del || vtx_visited[i]) continue;
        v = i << 1; /* either direction would work */
        kvec_t(uint64_t) vec = {0, 0, 0};
        kv_push(uint64_t, vec, v);
        vtx_visited[v>>1] = 1;
        while (asmg_arc_n1(g, v) > 0) {
            v = asmg_arc_a1(g, v)->w;
            if (!vtx_visited[v>>1])
                kv_push(uint64_t, vec, v);
            else
                break;
            vtx_visited[v>>1] = 1;
        }
        if (vec.n > 1) {
            kv_pushp(asmg_vtx_t, utgs, &vtx1);
            vtx1->n = vec.n;
            vtx1->a = vec.a;
            vtx1->seq = 0;
            vtx1->len = 0;
            vtx1->cov = vtx1->del = 0;
            vtx1->circ = 1;
            MYREALLOC(vtx1->a, vtx1->n);
        }
    }
    free(vtx_visited);

    /* singleton: UINT64_MAX; start: u << 1 | 0; end: u << 1 | 1; mid: UINT64_MAX - 1 */
    MYMALLOC(vtx_p, n_vtx);
    MYBONE(vtx_p, n_vtx);
    for (i = 0; i < utgs.n; ++i) {
        vtx1 = &utgs.a[i];
        vtx_p[vtx1->a[0]>>1] = i << 1;
        vtx_p[vtx1->a[vtx1->n-1]>>1] = i << 1 | 1;
        for (j = 1; j < vtx1->n-1; ++j) { /* vtx1->n > 1 */
            vtx_p[vtx1->a[j]>>1] = UINT64_MAX - 1;
            asmg_arc_del(g, vtx1->a[j-1], vtx1->a[j], 1);
            asmg_arc_del(g, vtx1->a[j]^1, vtx1->a[j-1]^1, 1);
        }
        asmg_arc_del(g, vtx1->a[vtx1->n-2], vtx1->a[vtx1->n-1], 1);
        asmg_arc_del(g, vtx1->a[vtx1->n-1]^1, vtx1->a[vtx1->n-2]^1, 1);
    }
    
    /* add singletons to utg set */
    for (i = 0; i < n_vtx; ++i) {
        if (vtx_p[i] == UINT64_MAX && !vtx[i].del) {
            vtx_p[i] = utgs.n << 1; /*update index*/
            kv_pushp(asmg_vtx_t, utgs, &vtx1);
            v = i << 1;
            vtx1->n = 1;
            MYMALLOC(vtx1->a, 1);
            vtx1->a[0] = v;
            vtx1->seq = 0;
            vtx1->len = 0;
            vtx1->cov = vtx1->del = 0;
            vtx1->circ = asmg_arc_exist1(g, v, v);
        }
    }

    /* add arcs */
    for (i = 0, n_arc = g->n_arc; i < n_arc; ++i) {
        arc = &g->arc[i];
        if (arc->del) continue;
        v = vtx_p[arc->v>>1];
        w = vtx_p[arc->w>>1];
        if (v == UINT64_MAX - 1 || w == UINT64_MAX - 1) continue;
        kv_pushp(asmg_arc_t, arcs, &arc1);
        arc1->v = utgs.a[v>>1].n>1? v^1 : v|(arc->v&1);
        arc1->w = utgs.a[w>>1].n>1?  w  : w|(arc->w&1);
        arc1->lo = arc->lo;
        arc1->link_id = arc->link_id;
        arc1->cov = arc->cov;
        arc1->del = 0;
        arc1->comp = arc->comp;
    }
    free(vtx_p);

    /* expand utgs */
    for (i = 0, n = utgs.n; i < n; ++i) {
        vtx_p = utgs.a[i].a;
        m = utgs.a[i].n;
        u64_v_t v_vec = {0, 0, 0};
        for (j = 0; j < m; ++j) {
            if (j > 0)
                v_vec.n -= asmg_arc(g, vtx_p[j-1], vtx_p[j])->lo;
            vtx1 = &vtx[vtx_p[j]>>1];
            vec_add(&v_vec, vtx1->a, vtx1->n, !!(vtx_p[j]&1));
        }
        vtx1 = &utgs.a[i];
        free(vtx1->a);
        vtx1->n = v_vec.n;
        vtx1->a = v_vec.a;
        vtx1->cov = 0; // need to redo alignment and coverage estimation
        MYREALLOC(vtx1->a, vtx1->n);
    }

    /* finalize */
    MYCALLOC(utg_asmg, 1);
    utg_asmg->n_vtx = utgs.n;
    utg_asmg->m_vtx = utgs.m;
    utg_asmg->vtx = utgs.a;
    utg_asmg->n_arc = arcs.n;
    utg_asmg->m_arc = arcs.m;
    utg_asmg->arc = arcs.a;

    asmg_finalize(utg_asmg, 1);

    return utg_asmg;
}

