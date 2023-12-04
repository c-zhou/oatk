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
#include "kdq.h"

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

static void asmg_cleanup(asmg_t *g) 
{
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

static uint32_t asmg_arc_fix_symm(asmg_t *g) 
{
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
            a1 = asmg_arc_add(g, w^1, v^1, a->ln, a->ls, a->link_id, a->cov, a->comp^1);
            ++symm_fix;
        } else {
            // fix symm flag
            a1->comp = a->comp^1;
            // FIXME overlap length difference
            if (a->ln != a1->ln)
                a->ln = a1->ln = MIN(a->ln, a1->ln);
            if (a->ls != a1->ls)
                a->ls = a1->ls = MIN(a->ls, a1->ls);
        }
    }

    return symm_fix;
}

void asmg_arc_fix_cov(asmg_t *g)
{
    uint64_t i, n;
    uint32_t c;
    asmg_arc_t *arc;
    for (i = 0, n = g->n_arc; i < n; ++i) {
        arc = &g->arc[i];
        if (arc->del) continue;
        c = MIN(g->vtx[arc->v >> 1].cov, g->vtx[arc->w >> 1].cov);
        if (c < arc->cov) arc->cov = c;
    }
}

void asmg_finalize(asmg_t *g, int do_cleanup) 
{
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

uint64_t *asmg_vtx_list(asmg_t *g, uint64_t *_n)
{
    uint64_t i;
    u64_v_t vlist = {0, 0, 0};
    for (i = 0; i < g->n_vtx; ++i)
        if (!g->vtx[i].del)
            kv_push(uint64_t, vlist, i);
    MYREALLOC(vlist.a, vlist.n);
    
    if (_n) *_n = vlist.n;

    return vlist.a;
}

void asmg_print(asmg_t *g, FILE *fo, int no_seq)
{
    uint64_t i, n;
    int l, c;
    char *s;
    asmg_vtx_t *vtx;
    asmg_arc_t *arc;

    fprintf(fo, "H\tVN:Z:1.0\n");
    for (i = 0, n = g->n_vtx; i < n; ++i) {
        vtx = &g->vtx[i];
        if (vtx->del) continue;
        l = vtx->len;
        s = vtx->seq;
        c = vtx->cov;
        if (s && (!no_seq))
            fprintf(fo, "S\tu%lu\t%.*s\tLN:i:%d\tKC:i:%ld\tSC:f:%.3f\n", i, l, s, l, (int64_t) l * c, (double) c);
        else
            fprintf(fo, "S\tu%lu\t*\tLN:i:%d\tKC:i:%ld\tSC:f:%.3f\n", i, l, (int64_t) l * c, (double) c);
    }

    for (i = 0, n = g->n_arc; i < n; ++i) {
        arc = &g->arc[i];
        if (arc->del || arc->comp) continue;
        fprintf(fo, "L\tu%lu\t%c\tu%lu\t%c\t%ldM\tEC:i:%u\n", arc->v>>1, "+-"[arc->v&1], arc->w>>1, "+-"[arc->w&1], arc->ls, arc->cov);
        fprintf(fo, "L\tu%lu\t%c\tu%lu\t%c\t%ldM\tEC:i:%u\n", arc->w>>1, "-+"[arc->w&1], arc->v>>1, "-+"[arc->v&1], arc->ls, arc->cov);
    }
}

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
    uint64_t i, nv, nv0, k, ls, min_l;
    asmg_arc_t *av;
    *l = 0, *w = UINT64_MAX;
    if (g->vtx[v>>1].del)
        return 0;
    nv0 = asmg_arc_n(g, v);
    av = asmg_arc_a(g, v);
    ls = k = 0;
    for (i = nv = 0; i < nv0; ++i) {
        if (!av[i].del) {
            ++nv;
            k = i;
            ls = av[i].ls > ls? av[i].ls : ls;
        }
    }
    min_l = g->vtx[v>>1].len - ls;
    *l = min_l;
    *w = nv == 1? av[k].w : UINT64_MAX;
    return nv;
}

static inline int32_t asmg_uext(asmg_t *g, uint64_t v, int32_t max_ext, uint64_t *ne, uint64_t *le, u64_v_t *a, int tip_only)
{
    // the extension of circular unitigs will be bounded by max_ext
    // and the return value will be ASMG_VT_MERGEABLE
    int32_t vt;
    uint64_t nv, nw, l, w, n_ext, l_ext;
    l = n_ext = l_ext = 0;
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
    if (tip_only && vt == ASMG_VT_MULTI_OUT) {
        // for tip extenstion
        l_ext -= l;
        if (a) (*a).n--;
    }
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
    n_arc = asmg_max_link_id(g) + 1; // number unique arcs
    
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
            vt = asmg_uext(g, v, n_vtx*2+1, 0, 0, &a, 0);
            for (j = 1; j < a.n; ++j) {
                arc_group[asmg_arc1(g, a.a[j-1], a.a[j])->link_id] = group;
                visited[a.a[j]>>1] = 1;
                ++na;
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
            l = g->vtx[v>>1].len - av[i].ls;
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
#undef DEBUG_EXEC_ORDER

#ifdef DEBUG_EXEC_ORDER
#include <time.h>
void shuffle(uint64_t *a, size_t n)
{
    while (n > 1) {
        size_t k = (size_t) ((double)(n--) * (rand() / (RAND_MAX+1.0)));
        SWAP(a[n], a[k]);
    }
}
#endif

static uint64_t asmg_cwt_len(asmg_t *g, uint64_t *v, uint64_t nv)
{
    if (nv == 0) return 0;
    uint64_t i, wt_l, ov_l;
    wt_l = g->vtx[v[0]>>1].len * g->vtx[v[0]>>1].cov;
    for (i = 1; i < nv; ++i) {
        ov_l = asmg_arc(g, v[i-1], v[i])->ls;
        wt_l += (g->vtx[v[i]>>1].len - ov_l) * g->vtx[v[i]>>1].cov;
    }
    return wt_l;
}

// super tip: a tip with a coverage larger than [half] the coverage of the vertex it attached to
uint64_t asmg_drop_tip(asmg_t *g, int32_t tip_cnt, uint64_t tip_len, int protect_super_tip, int do_cleanup, int VERBOSE)
{
    uint64_t i, v, w, w1, n1, n_vtx, cnt, b_tip, c_tip, l_ext;
    int32_t vt;
    int is_tip;
    asmg_arc_t *a1;
    u64_v_t a = {0, 0, 0};
    u64_v_t b = {0, 0, 0};
    u64_v_t d = {0, 0, 0};
    n_vtx = asmg_vtx_n(g);
    if ((uint64_t) tip_cnt > n_vtx)
        tip_cnt = n_vtx;
    cnt = 0;
#ifdef DEBUG_EXEC_ORDER
    srand((unsigned) time(NULL));
    uint64_t j, idx[n_vtx];
    for (v = 0; v < n_vtx; ++v) idx[v] = v;
    shuffle(idx, n_vtx);
    for (j = 0; j < n_vtx; ++j) {
        v = idx[j];
#else
    for (v = 0; v < n_vtx; ++v) {
#endif
        if (g->vtx[v>>1].del) continue;
        if (asmg_arc_n1(g, v^1) != 0) continue; // not a tip
        vt = asmg_uext(g, v, tip_cnt, 0, &l_ext, &a, 1);
        if (a.n == 0) continue; // v is ASMG_VT_MULTI_OUT
        if (vt == ASMG_VT_MERGEABLE) continue; // circular unitig
        if (l_ext > tip_len) continue; // tip too long
        if (vt != ASMG_VT_TIP && protect_super_tip) {
            w = a.a[a.n - 1]; // the last vertex
            // asmg_arc_n1(g, w) is 0 or 1 by definition of tips
            // zero for a tip subgraph
            // assert(asmg_arc_n1(g, w) == 1);
            b_tip = l_ext;
            c_tip = asmg_cwt_len(g, a.a, a.n);
            /***
            for (i = 0; i < a.n; ++i) {
                b_tip += g->vtx[a.a[i]>>1].len;
                c_tip += g->vtx[a.a[i]>>1].len * g->vtx[a.a[i]>>1].cov;
            }
            **/
            // a tip with high coverage
            // if (c_tip * 2 > b_tip * asmg_arc_a1(g, w)->cov) continue;
            // compare max coverage neighbour vtx
            w1 = asmg_arc_a1(g, w)->w ^ 1; // this is the only vtx by definition
            a1 = asmg_arc_a(g, w1);
            n1 = asmg_arc_n(g, w1);
            is_tip = 0;
            for (i = 0; i < n1; ++i) {
                if ((a1[i].del || a1[i].w ^ 1) == w) 
                    continue;
                asmg_uext(g, a1[i].w, n_vtx + 1, 0, &l_ext, &b, 0);
                if (b_tip <= l_ext || c_tip * 2 <= asmg_cwt_len(g, b.a, b.n)) {
                    is_tip = 1;
                    break;
                }
            }
            if (!is_tip) continue;
        }
        kv_pushn(uint64_t, d, a.a, a.n);
        ++cnt;
    }
    for (i = 0; i < d.n; ++i)
        asmg_vtx_del(g, d.a[i]>>1, 1);
    free(a.a);
    free(b.a);
    free(d.a);
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

uint64_t asmg_remove_weak_crosslink(asmg_t *g, double c_thresh, double m_cov, int do_cleanup, int VERBOSE)
{
    uint64_t i, k, n, n1, v, w, cnt = 0;
    // uint32_t cv, cw;
    int weak;
    asmg_arc_t *a, *a1;
    u64_v_t d = {0, 0, 0};

#ifdef DEBUG_EXEC_ORDER
    srand((unsigned) time(NULL));
    n = g->n_arc;
    uint64_t j, idx[n];
    for (j = 0; j < n; ++j) idx[j] = j;
    shuffle(idx, n);
    for (j = 0; j < n; ++j) {
        i = idx[j];
#else
    for (i = 0, n = g->n_arc; i < n; ++i) {
#endif
        a = &g->arc[i];
        if (a->del || a->comp) continue;
        v = a->v;
        w = a->w;
        /***
        if (asmg_arc_n1(g, v^1) != 1 || asmg_arc_n1(g, w) != 1)
            continue;
        if (asmg_arc_n3(g, v, w, &cv) != 2 || 
                asmg_arc_n3(g, w^1, v^1, &cw) != 2)
            continue;
        cv = MIN(cv, cw);
        if (cv > 0 && (double) a->cov / cv < c_thresh) {
            kv_push(uint64_t, d, i);
            ++cnt;
        }
        **/
        
        // check if dominating outgoing arc exists
        n1 = asmg_arc_n(g, v);
        a1 = asmg_arc_a(g, v);
        weak = 0;
        for (k = 0; k < n1; ++k) {
            if (a1[k].del || a1[k].cov < m_cov)
                continue;
            if ((double) a->cov / a1[k].cov < c_thresh) {
                weak = 1;
                break;
            }
        }
        if (!weak) continue;
        
        // check if dominating incoming arc exists
        n1 = asmg_arc_n(g, w^1);
        a1 = asmg_arc_a(g, w^1);
        weak = 0;
        for (k = 0; k < n1; ++k) {
            if (a1[k].del || a1[k].cov < m_cov)
                continue;
            if ((double) a->cov / a1[k].cov < c_thresh) {
                weak = 1;
                break;
            }
        }
        if (!weak) continue;

        kv_push(uint64_t, d, i);
        ++cnt;
    }
    for (i = 0, n = d.n; i < n; ++i) {
        a = &g->arc[d.a[i]];
        a->del = 1;
        asmg_arc_del(g, a->w^1, a->v^1, 1);
    }
    free(d.a);
    if (do_cleanup && cnt > 0) asmg_finalize(g, 1);
    if (VERBOSE)    
        fprintf(stderr, "[M::%s] dropped %lu weak cross links\n", __func__, cnt);
    return cnt;
}

/******************
 * Bubble popping *
 ******************/

// in a resolved bubble, mark unused vertices and arcs as "reduced"
static int asmg_bub_backtrack(asmg_t *g, uint64_t v0, uint64_t max_del, int protect_super_bubble, asmg_tbuf_t *b)
{
    uint64_t i, v, w;
    asmg_arc_t *a;

    assert(b->S.n == 0);
    if (max_del > 0) {
        uint64_t n_kept = 0;
        v = b->v_sink;
        do { ++n_kept, v = b->a[v].p; } while (v != v0);
        if (b->b.n > n_kept + max_del) return 0;
    }
    if (protect_super_bubble) {
        uint64_t n_kept, b_kept, c_kept, b_tot, c_tot;
        n_kept = b_kept = c_kept = 0;
        v = b->v_sink;
        do { ++n_kept, b_kept += g->vtx[v>>1].len, c_kept += g->vtx[v>>1].len * g->vtx[v>>1].cov, v = b->a[v].p; } while (v != v0);
        b_tot = c_tot = 0;
        for (i = 0; i < b->b.n; ++i) b_tot += g->vtx[b->b.a[i]>>1].len, c_tot += g->vtx[b->b.a[i]>>1].len * g->vtx[b->b.a[i]>>1].cov;
        // check length for super bubble protection: b_delt / n_delt * 2 > b_source + b_sink
        uint64_t le, re, le_wt, re_wt; // left and right extension 
        u64_v_t a = {0, 0, 0};
        le = re = le_wt = re_wt = 0;
        asmg_uext(g, v0^1,      g->n_vtx*2+1, 0, &le, &a, 0);
        le_wt = asmg_cwt_len(g, a.a, a.n);
        asmg_uext(g, b->v_sink, g->n_vtx*2+1, 0, &re, &a, 0);
        re_wt = asmg_cwt_len(g, a.a, a.n);
        free(a.a);
        //if ((b_tot - b_kept) * 2 > (g->vtx[v0>>1].len + g->vtx[b->v_sink>>1].len) * (b->b.n - n_kept)) return 0;
        // if ((b_tot - b_kept) * 2 > (le + re) * (b->b.n - n_kept)) return 0;
        if ((c_tot - c_kept) * (le + re) * 2 > (le_wt + re_wt) * (b_tot - b_kept)) return 0;
        // check coverage for super bubbble protection: c_delt / b_delt * 2 > c_kept / b_kept
        if ((c_tot - c_kept) * b_kept * 2 > c_kept * (b_tot - b_kept)) return 0;
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
    return 1;
}

// pop bubbles from vertex v0; the graph MUST BE symmetric: if u->v present, v'->u' must be present as well
// radius is calculated from the end of v0, not the start
static uint64_t asmg_bub_pop1(asmg_t *g, uint64_t v0, uint64_t radius, uint64_t max_del, int protect_tip, int protect_super_bubble, asmg_tbuf_t *b)
{
    uint64_t ret = 0;
    if (asmg_arc_n1(g, v0) < 2) return 0; // no bubbles
    asmg_topo_ext(g, v0, g->vtx[v0>>1].len + radius, protect_tip? 0 : ASMG_TE_THRU_SHORT_TIP, b);
    if (b->n_sink) {
        ret = asmg_bub_backtrack(g, v0, max_del, protect_super_bubble, b);
        if (ret) ret |= (uint64_t) b->n_short_tip << 32;
    }
    asmg_tbuf_reset(b);
    return ret;
}

// pop bubbles
// super bubble: the average length of sequences to be removed is greater than the average length of source and sink sequence?
// super bubble: the average coverage of sequences to be removed is greater than [half] the average coverage of source and sink sequence
// super bubble: the average coverage of sequences to be removed is greater than [half] the average coverage of sequences to be kept
// super buubles are protected with protect_super_bubble set
// TODO does vertex traversal order matters? in theory the number of bubbles being processed could be different but the resulted graph should be same
uint64_t asmg_pop_bubble(asmg_t *g, uint64_t radius, uint64_t max_del, int protect_tip, int protect_super_bubble, int do_cleanup, int VERBOSE)
{
    uint64_t v, n_vtx, n_pop;
    asmg_tbuf_t *b;

    n_vtx = asmg_vtx_n(g);
    b = asmg_tbuf_init(g);
    n_pop = 0;

#ifdef DEBUG_EXEC_ORDER
    srand((unsigned) time(NULL));
    uint64_t j, idx[n_vtx];
    for (j = 0; j < n_vtx; ++j) idx[j] = j;
    shuffle(idx, n_vtx);
    for (j = 0; j < n_vtx; ++j) {
        v = idx[j];
#else
    for (v = 0; v < n_vtx; ++v) {
#endif
        if (!g->vtx[v>>1].del && asmg_arc_n1(g, v) >= 2)
            n_pop += asmg_bub_pop1(g, v, radius, max_del, protect_tip, protect_super_bubble, b);
    }
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
    kvec_t(asmg_vtx_t) utgs = {0, 0, 0};
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
                } else {
                    kv_destroy(vec);
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
        } else {
            kv_destroy(vec);
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
        } else {
            kv_destroy(vec);
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
        arc1->ln = arc->ln;
        arc1->ls = arc->ls;
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
                v_vec.n -= asmg_arc(g, vtx_p[j-1], vtx_p[j])->ln;
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

KDQ_INIT(uint64_t)

typedef kdq_t(uint64_t) kdq_u64_t;

uint32_t *asmg_subgraph(asmg_t *g, uint32_t *seeds, uint32_t n, uint32_t step, 
        uint64_t dist, uint32_t *_nv, int modify_graph)
{
    // make a subgraph from given a vertex set and radius
    // if modify_graph
    //     all vertices and arcs in other components will be marked as deleted
    //     return 0
    // else
    //     keep graph unchanged
    //     return the list of vertices in the subgraph
    if (n == 0) {
        if (_nv) *_nv = 0;
        return 0;
    }

    uint32_t r;
    uint64_t i, v, x, rd, nv;
    int8_t *flag;
    kdq_u64_t *q;
    kdq_u64_t *d;
    asmg_arc_t *av;

    if (step == 0)
        step = UINT32_MAX;
    if (dist == 0)
        dist = UINT64_MAX;

    MYCALLOC(flag, asmg_vtx_n(g));
    // flag deleted segments
    for (i = 0; i < g->n_vtx; ++i) {
        if (g->vtx[i].del) {
            flag[i<<1|0] = -1;
            flag[i<<1|1] = -1;
        }
    }

    q = kdq_init(uint64_t, 0);
    d = kdq_init(uint64_t, 0);
    for (i = 0; i < n; ++i) {
        if (seeds[i] < g->n_vtx) {
            kdq_push(uint64_t, q, ((uint64_t)seeds[i]<<1|0)<<32);
            kdq_push(uint64_t, d, 0);
            kdq_push(uint64_t, q, ((uint64_t)seeds[i]<<1|1)<<32);
            kdq_push(uint64_t, d, 0);
        }
    }

    if (modify_graph) {
        for (i = 0; i < g->n_vtx; ++i) // mark all segments to be deleted
            g->vtx[i].del = 1;
        // for (i = 0; i < g->n_arc; ++i) // mark all arcs to be deleted
        //     g->arc[i].del = 1;
    }

    while (kdq_size(q) > 0) {
        x = *kdq_shift(uint64_t, q);
        v = x>>32;
        r = (uint32_t) x;
        rd = *kdq_shift(uint64_t, d);
        if (flag[v] != 0) continue; // already visited or deleted
        flag[v] = 1;
        if (modify_graph) g->vtx[v>>1].del = 0;
        if (r < step && rd < dist) {
            nv = asmg_arc_n(g, v);
            av = asmg_arc_a(g, v);
            for (i = 0; i < nv; ++i) {
                if (av[i].del) continue;
                // if (modify_graph) av[i].del = 0;
                if (flag[av[i].w] == 0) {
                    kdq_push(uint64_t, q, (uint64_t)av[i].w<<32 | (r + 1));
                    kdq_push(uint64_t, d, rd + g->vtx[av[i].w>>1].len - av[i].ls);
                }
                if (flag[av[i].w^1] == 0) {
                    kdq_push(uint64_t, q, (uint64_t)(av[i].w^1)<<32 | (r + 1));
                    kdq_push(uint64_t, d, rd + g->vtx[av[i].w>>1].len - av[i].ls);
                }
            }
        }
    }
    assert(kdq_size(d) == 0);

    // change flag index from vtx<<1 to vtx
    for (i = 0; i < g->n_vtx; ++i)
        flag[i] = (flag[i<<1|0] > 0 || flag[i<<1|1] > 0);

    uint32_t *vs;
    nv = 0;
    vs = 0;
    if (!modify_graph) {
        // make vlist
        kvec_t(uint32_t) vlist;
        kv_init(vlist);
        for (i = 0; i < g->n_vtx; ++i)
            if (flag[i])
                kv_push(uint32_t, vlist, i);
        MYREALLOC(vlist.a, vlist.n);
        vs = vlist.a;
        nv = vlist.n;
    } else {
        // delete dirty arcs
        for (i = 0; i < g->n_arc; ++i) {
            if (!flag[g->arc[i].v>>1] || !flag[g->arc[i].w>>1])
                g->arc[i].del = 1;
        }
        if (_nv) {
            for (i = 0; i < g->n_vtx; ++i)
                nv += flag[i];
        }
    }
    
    if (_nv) *_nv = nv;
    
    kdq_destroy(uint64_t, q);
    kdq_destroy(uint64_t, d);
    free(flag);

    return vs;
}

int asmg_path_exists(asmg_t *g, uint32_t source, uint32_t sink, uint32_t step, uint64_t dist, uint32_t *_step, uint64_t *_dist)
{
    // TODO dealing with graphs with deleted vertices and arcs
    // detect path between two vertices
    if (source >= asmg_vtx_n(g) || sink >= asmg_vtx_n(g))
        return 0;
    
    if (_step) *_step = 0;
    if (_dist) *_dist = 0;

    uint32_t i, v, r, nv;
    uint64_t x, rd;
    int8_t *flag;
    kdq_u64_t *q;
    kdq_u64_t *d;
    asmg_arc_t *av;

    int exists = 0;

    if (step == 0)
        step = UINT32_MAX;
    if (dist == 0)
        dist = UINT64_MAX;
    q = kdq_init(uint64_t, 0);
    d = kdq_init(uint64_t, 0);
    MYCALLOC(flag, asmg_vtx_n(g));
    kdq_push(uint64_t, q, (uint64_t)source<<32);
    kdq_push(uint64_t, d, 0);
    while (kdq_size(q) > 0) {
        x = *kdq_shift(uint64_t, q);
        v = x>>32;
        r = (uint32_t) x;
        rd = *kdq_shift(uint64_t, d);
        if (flag[v]) continue; // already visited
        flag[v] = 1;
        if (r < step && rd < dist) {
            nv = asmg_arc_n(g, v);
            av = asmg_arc_a(g, v);
            for (i = 0; i < nv; ++i) {
                if (av[i].w == sink) {
                    exists = 1;
                    if (_step) *_step = r;
                    if (_dist) *_dist = rd;
                    goto do_clean;
                }
                if (flag[av[i].w] == 0) {
                    kdq_push(uint64_t, q, (uint64_t)av[i].w<<32 | (r + 1));
                    kdq_push(uint64_t, d, rd + g->vtx[av[i].w>>1].len - av[i].ls);
                }
            }
        }
    }
    assert(kdq_size(d) == 0);

do_clean:
    kdq_destroy(uint64_t, q);
    kdq_destroy(uint64_t, d);
    free(flag);

    return exists;
}

static void scc_util(asmg_t *g, uint64_t v, int *disc, int *low, kdq_u64_t *st, int *stb, int *depth, int *scc, int *n_scc)
{
    uint64_t i, w, n;
    asmg_arc_t *a;

    disc[v] = low[v] = ++*depth;
    kdq_push(uint64_t, st, v);
    stb[v] = 1;

    a = asmg_arc_a(g, v);
    n = asmg_arc_n(g, v);
    for (i = 0; i < n; ++i) {
        if (a[i].del) continue;
        w = a[i].w;
        if (g->vtx[w>>1].del) continue;
        if (disc[w] == -1) {
            scc_util(g, w, disc, low, st, stb, depth, scc, n_scc);
            low[v] = MIN(low[v], low[w]);
        } else if (stb[w] == 1) {
            low[v] = MIN(low[v], disc[w]);
        }
    }
    
    if (low[v] == disc[v]) {
        do {
            w = *kdq_pop(uint64_t, st);
            stb[w] = 0;
            scc[w] = *n_scc;
        } while (w != v);
        ++*n_scc;
    }
}

int asmg_tarjans_scc(asmg_t *g, int *scc)
{
    uint64_t i, n_seg;
    int n_scc, depth, *low, *disc, *stb;
    kdq_u64_t *st;

    n_seg = asmg_vtx_n(g);
    MYMALLOC(low, n_seg);
    MYMALLOC(disc, n_seg);
    MYMALLOC(stb, n_seg);
    st = kdq_init(uint64_t, 0);
    for (i = 0; i < n_seg; ++i) {
        scc[i] = -1;
        low[i] = -1;
        disc[i] = -1;
        stb[i] = 0;
    }
    
    n_scc = depth = 0;
    for (i = 0; i < n_seg; ++i)
        if (disc[i] == -1 && !g->vtx[i>>1].del)
            scc_util(g, i, disc, low, st, stb, &depth, scc, &n_scc);

    free(low);
    free(disc);
    free(stb);
    kdq_destroy(uint64_t, st);
    
    return n_scc;
}

