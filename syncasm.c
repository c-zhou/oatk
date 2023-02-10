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
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <pthread.h>
#include <assert.h>

#include "kvec.h"
#include "kdq.h"
#include "khashl.h"
#include "kstring.h"

#include "syncmer.h"
#include "syncasm.h"
#include "graph.h"

#undef DEBUG_KMER_OVERLAP
#undef DEBUG_CONSENSUS
#undef DEBUG_SCM_UTG_INDEX
#undef DEBUG_UTG_MULTIPLEX
#undef DEBUG_UTG_DEMULTIPLEX

static kh_inline khint_t kh_hash_uint128(uint128_t key)
{
    khint_t k1 = kh_hash_uint64((khint64_t) key);
    khint_t k2 = kh_hash_uint64((khint64_t) (key >> 64));
    return kh_hash_uint64((khint64_t) ((uint64_t) k1 << 32 | k2));
}

static kh_inline khint_t kh_hash_generic(int key)
{
    return (khint_t) key;
}

KHASHL_MAP_INIT(KH_LOCAL, kh_scm_t, kh_scm, uint128_t, uint64_t, kh_hash_uint128, kh_eq_generic)
KHASHL_MAP_INIT(KH_LOCAL, kh_generic_t, kh_generic, int, int, kh_hash_generic, kh_eq_generic)
KHASHL_MAP_INIT(KH_LOCAL, kh_u64_t, kh_u64, uint64_t, uint64_t, kh_hash_uint64, kh_eq_generic)
KHASHL_MAP_INIT(KH_LOCAL, kh_dbl_t, kh_dbl, uint128_t, double, kh_hash_uint128, kh_eq_generic)
KHASHL_SET_INIT(KH_LOCAL, ks_u128_t, ks_u128, uint128_t, kh_hash_uint128, kh_eq_generic)

static void scg_destroy_scm(scg_t *g)
{
    if (!g->scm) return;
    size_t i;
    for (i = 0; i < g->n_scm; ++i) {
        if (g->scm[i].m_pos)
            free(g->scm[i].m_pos);
    }
    free(g->scm);
}

void scg_destroy(scg_t *g)
{
    if (!g) return;
    scg_destroy_scm(g);
    if (g->h_scm) kh_scm_destroy((kh_scm_t *) g->h_scm);
    if (g->utg_asmg) asmg_destroy(g->utg_asmg);
    if (g->scm_u) free(g->scm_u);
    if (g->idx_u) free(g->idx_u);
    free(g);
}

typedef struct {uint128_t h; uint64_t s; uint64_t m_pos;} syncmer1_t;

static int syncmer1_h_cmpfunc(const void *a, const void *b)
{
    uint128_t x, y;
    uint64_t s, t;
    x = ((syncmer1_t *) a)->h;
    y = ((syncmer1_t *) b)->h;
    s = ((syncmer1_t *) a)->s;
    t = ((syncmer1_t *) b)->s;
    return x == y? (s == t? 0 : (s > t? 1 : -1)) : (x > y? 1 : -1);
}

static int u128_cmpfunc(const void *a, const void *b)
{
    return (*(uint128_t *) a > *(uint128_t *) b) - (*(uint128_t *) a < *(uint128_t *) b);
}

static void scg_scm_utg_index(scg_t *g)
{
    if (g->idx_u) free(g->idx_u);
    if (g->scm_u) free(g->scm_u);

    kvec_t(uint128_t) scm_u;

    if (!g->utg_asmg) return;

    size_t i, j, m, n;
    uint64_t *a;
    scg_utg_t *utg;
    
    utg = g->utg_asmg->vtx;
    kv_init(scm_u);
    for (i = 0, n = g->utg_asmg->n_vtx; i < n; ++i) {
        if (utg[i].del)
            continue;
        a = utg[i].a;
        for (j = 0, m = utg[i].n; j < m; ++j)
            kv_push(uint128_t, scm_u, ((uint128_t) a[j] << 78) | ((uint128_t) i << 36) | j);
    }

    if (scm_u.n == 0)
        return;

    // scm position array sorting
    qsort(scm_u.a, scm_u.n, sizeof(uint128_t), u128_cmpfunc);

    // scm position array indexing
    uint128_t **idx_u, *p, *p_n;
    uint64_t s0, s;
    MYMALLOC(idx_u, g->n_scm + 1);
    p = scm_u.a;
    s = scm_utg_sid(p[0]);
    p_n = p + scm_u.n;
    i = 0;
    n = g->n_scm;
    while (i < n && p < p_n) {
        while (i <= s)
            idx_u[i++] = p;
        s0 = s;
        while (s0 == s && ++p < p_n)
            s = scm_utg_sid(p[0]);
    }
    while (i <= n)
        idx_u[i++] = p;

    assert(scm_u.a + scm_u.n == idx_u[g->n_scm]);

    g->idx_u = idx_u;
    g->scm_u = scm_u.a;

#ifdef DEBUG_SCM_UTG_INDEX
    for (i = 0; i < n; ++i) {
        m = scm_utg_n(g, i);
        if (m > 0) {
            p = idx_u[i];
            s = scm_utg_sid(p[0]);
            for (j = 1; j < m; ++j) {
                if (s != scm_utg_sid(p[j])) {
                    fprintf(stderr, "[DEBUG_SCM_UTG_INDEX::%s] syncmer unitigging index error!\n", __func__);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
#endif
}

static uint64_t scg_get_scm_id(scg_t *g, uint128_t key)
{
    kh_scm_t *h = (kh_scm_t *) g->h_scm;
    khint_t k = kh_scm_get(h, key);
    if (k == kh_end(h))
        return UINT64_MAX;
    return kh_val(h, k);
}

static inline void add_a_cov(kh_scm_t *h, uint128_t v_p, uint64_t c)
{
    int absent;
    khint_t k;
    k = kh_scm_put(h, v_p, &absent);
    if (absent)
        kh_val(h, k) = 0;
    kh_val(h, k) += c;
    return;
}

static int uint64_cmpfunc(const void *a, const void *b)
{
    uint64_t x, y;
    x = *(uint64_t *) a;
    y = *(uint64_t *) b;
    return (x > y) - (x < y);
}


scg_t *make_syncmer_graph(sr_v *sr, uint32_t min_k_cov, uint32_t min_a_cov)
{
    size_t i, j;
    uint32_t c;
    uint64_t n1, n2;
    sr_t *s;
    kvec_t(syncmer1_t) scm1;
    kv_init(scm1);
    n1 = 0;
    for (i = 0; i < sr->n; ++i) {
        s = &sr->a[i];
        assert(s->sid == i);
        n1 += s->n;
        for (j = 0; j < s->n; ++j) {
            syncmer1_t s1 = {s->k_mer_h[j], s->s_mer[j], (s->sid << 25) | (j << 1) | (s->m_pos[j] & 1)};
            kv_push(syncmer1_t, scm1, s1);
        }
    }
    if (scm1.n == 0) {
        kv_destroy(scm1);
        return 0;
    }

    scg_t *scg;
    MYCALLOC(scg, 1);

    qsort(scm1.a, scm1.n, sizeof(syncmer1_t), syncmer1_h_cmpfunc);

    // pack syncmers by kmer hash128
    kvec_t(syncmer_t) scm;
    kv_init(scm);
    uint128_t h128 = scm1.a[0].h;
    uint64_t s64 = scm1.a[0].s;
    c = 1;
    for (i = 1; i < scm1.n; ++i) {
        if (scm1.a[i].h == h128) {
            ++c;
            // syncmers with same hash values should have same smers
            if (scm1.a[i].s != s64) {
                fprintf(stderr, "[E::%s] identical kmers have different smers\n", __func__);
                fprintf(stderr, "[E::%s] kmer hash: %064lu%064lu\n", __func__, (uint64_t) (h128>>64), (uint64_t) h128);
                fprintf(stderr, "[E::%s] smer code: %lu\n", __func__, s64);
                for (j = i + 1 - c; j <= i; ++j) {
                    syncmer1_t *m = &scm1.a[j];
                    fprintf(stderr, "[E::%s] %064lu%064lu %lu %lu %llu %u\n", __func__,
                            (uint64_t) (m->h>>64), (uint64_t) m->h, m->s, m->m_pos>>25, m->m_pos>>1&MAX_RD_SCM,
                            sr->a[m->m_pos>>25].m_pos[m->m_pos>>1&MAX_RD_SCM]>>1);
                }
                exit(EXIT_FAILURE);
            }
        } else {
            uint64_t *m_pos;
            MYMALLOC(m_pos, c);
            for (j = i - c; j < i; ++j)
                m_pos[j+c-i] = scm1.a[j].m_pos;
            syncmer_t s1 = {h128, s64, m_pos, c, c < min_k_cov};
            kv_push(syncmer_t, scm, s1);

            h128 = scm1.a[i].h;
            s64 = scm1.a[i].s;
            c = 1;
        }
    }
    uint64_t *m_pos;
    MYMALLOC(m_pos, c);
    for (j = i - c; j < i; ++j)
        m_pos[j+c-i] = scm1.a[j].m_pos;
    syncmer_t s1 = {h128, s64, m_pos, c, c < min_k_cov};
    kv_push(syncmer_t, scm, s1);

    n2 = 0;
    for (i = 0; i < scm.n; ++i) n2 += scm.a[i].k_cov;

    assert(n1 == n2);

    kv_destroy(scm1);

    scg->scm = scm.a;
    scg->m_scm = scg->n_scm = scm.n;
    MYREALLOC(scg->scm, scg->n_scm);

    // add utg vtx
    // each single scm is a utg
    asmg_vtx_t *vtx, *vtx1;
    MYCALLOC(vtx, scg->n_scm);
    for (i = 0; i < scg->n_scm; ++i) {
        vtx1 = &vtx[i];
        MYMALLOC(vtx1->a, 1);
        vtx1->n = 1;
        vtx1->a[0] = i<<1;
        vtx1->cov = scg->scm[i].k_cov;
        vtx1->del = scg->scm[i].del;
    }

    // sort syncmer m_pos
    for (i = 0; i < scg->n_scm; ++i)
        qsort(scg->scm[i].m_pos, scg->scm[i].k_cov, sizeof(uint64_t), uint64_cmpfunc);

    // build syncmer hash h128 to id map
    int absent;
    khint_t k;
    kh_scm_t *h_scm;
    h_scm = kh_scm_init();
    for (i = 0; i < scm.n; ++i) {
        k = kh_scm_put(h_scm, scm.a[i].h, &absent);
        kh_val(h_scm, k) = i;
    }
    scg->h_scm = h_scm;

    // add arcs
    kvec_t(asmg_arc_t) arc;
    kv_init(arc);
    uint64_t v0, v1;
    uint32_t r0, r1, v_v;
    uint128_t v_p; // v0, v1 packed
    kh_scm_t *a_cov; // arc cov
    a_cov = kh_scm_init();
    for (i = 0; i < sr->n; ++i) {
        s = &sr->a[i];
        if (s->n == 0) continue;
        v0 = scg_get_scm_id(scg, s->k_mer_h[0]), assert(v0 != UINT64_MAX);
        r0 = s->m_pos[0] & 1;
        v0 = v0 << 1 | r0;

        for (j = 1; j < s->n; ++j) {
            v1 = scg_get_scm_id(scg, s->k_mer_h[j]), assert(v1 != UINT64_MAX);
            r1 = s->m_pos[j] & 1;
            v1 = v1 << 1 | r1;

            v0 <= v1? add_a_cov(a_cov, (uint128_t) v0 << 64 | v1, 1) :
                add_a_cov(a_cov, ((uint128_t) v1^1) << 64 | (v0^1), 1);

            v0 = v1;
        }
    }

    for (k = (khint_t) 0; k < kh_end(a_cov); ++k) {
        if (kh_exist(a_cov, k)) {
            v_p = kh_key(a_cov, k);
            v_v = kh_val(a_cov, k);
            v0 = (uint64_t) (v_p >> 64);
            v1 = (uint64_t) v_p;
            if (v_v < min_a_cov || vtx[v0>>1].del || vtx[v1>>1].del)
                continue;
            asmg_arc_t a = {v0, v1, 0, UINT64_MAX, v_v, 0, 0};
            kv_push(asmg_arc_t, arc, a);
            if ((v1^1) != v0 || (v0^1) != v1) {
                // to avoid multi-arcs
                // for arcs like (v+)->(v-)
                asmg_arc_t a_c = {v1^1, v0^1, 0, UINT64_MAX, v_v, 0, 1};
                kv_push(asmg_arc_t, arc, a_c);
            }
        }
    }
    kh_scm_destroy(a_cov);

    asmg_t *utg_asmg;
    MYCALLOC(utg_asmg, 1);
    utg_asmg->vtx = vtx;
    utg_asmg->m_vtx = utg_asmg->n_vtx = scg->n_scm;
    utg_asmg->arc = arc.a;
    utg_asmg->m_arc = utg_asmg->n_arc = arc.n;
    MYREALLOC(utg_asmg->arc, utg_asmg->n_arc);
    // sort arc and build index
    asmg_finalize(utg_asmg, 1);

    scg->utg_asmg = utg_asmg;
    scg_scm_utg_index(scg);

    return scg;
}

void scg_arc_coverage(scg_t *scg, sr_v *sr)
{
    // add arcs
    uint64_t i, j, n;
    uint64_t v0, v1;
    uint32_t r0, r1;
    khint_t k;
    kh_scm_t *a_cov; // arc cov
    asmg_t *g;
    asmg_arc_t *a;
    asmg_vtx_t *v;
    sr_t *s;

    g = scg->utg_asmg;
    a_cov = kh_scm_init();
    for (i = 0, n = g->n_arc; i < n; ++i) {
        a = &g->arc[i];
        if (a->del) continue;
        v = &g->vtx[a->v>>1];
        v0 = a->v&1? v->a[0]^1 : v->a[v->n-1];
        v = &g->vtx[a->w>>1];
        v1 = a->w&1? v->a[v->n-1]^1 : v->a[0];
        add_a_cov(a_cov, (uint128_t) v0 << 64 | v1, 0);
    }
    for (i = 0, n = sr->n; i < n; ++i) {
        s = &sr->a[i];
        if (s->n == 0) continue;
        v0 = scg_get_scm_id(scg, s->k_mer_h[0]), assert(v0 != UINT64_MAX);
        r0 = s->m_pos[0] & 1;
        v0 = v0 << 1 | r0;

        for (j = 1; j < s->n; ++j) {
            v1 = scg_get_scm_id(scg, s->k_mer_h[j]), assert(v1 != UINT64_MAX);
            r1 = s->m_pos[j] & 1;
            v1 = v1 << 1 | r1;

            k = kh_scm_get(a_cov, (uint128_t) v0 << 64 | v1);
            if (k < kh_end(a_cov)) {
                kh_val(a_cov, k) += 1;
                if ((v1^1) != v0 || (v0^1) != v1) {
                    k = kh_scm_get(a_cov, (uint128_t) (v1^1) << 64 | (v0^1));
                    kh_val(a_cov, k) += 1;
                }
            }

            v0 = v1;
        }
    }
    for (i = 0, n = g->n_arc; i < n; ++i) {
        a = &g->arc[i];
        if (a->del) continue;
        v = &g->vtx[a->v>>1];
        v0 = a->v&1? v->a[0]^1 : v->a[v->n-1];
        v = &g->vtx[a->w>>1];
        v1 = a->w&1? v->a[v->n-1]^1 : v->a[0];
        k = kh_scm_get(a_cov, (uint128_t) v0 << 64 | v1);
        a->cov = kh_val(a_cov, k);
    }
    kh_scm_destroy(a_cov);
}

int scg_is_empty(scg_t *scg)
{
    uint64_t i, n, n_scm;
    n_scm = scg->n_scm;
    for (i = 0, n = scg->n_scm; i < n; ++i)
        if (scg->scm[i].del) --n_scm;
    return n_scm == 0;
}

void scg_stat(scg_t *scg, FILE *fo, uint64_t *stats)
{
    uint64_t i, n, n_scm, u_scm, n_utg, n_arc;
    asmg_t *utg_asmg;
    utg_asmg = scg->utg_asmg;
    n_utg = utg_asmg->n_vtx;
    n_arc = utg_asmg->n_arc;
    n_scm = 0;
    u_scm = scg->n_scm;
    
    for (i = 0, n = scg->n_scm; i < n; ++i)
        if (scg->scm[i].del)
            --u_scm;
    for (i = 0, n = utg_asmg->n_arc; i < n; ++i)
        if (utg_asmg->arc[i].del)
            --n_arc;
    for (i = 0, n = utg_asmg->n_vtx; i < n; ++i) {
        if (utg_asmg->vtx[i].del)
            --n_utg;
        else
            n_scm += utg_asmg->vtx[i].n;
    }

    if (fo != 0) {
        fprintf(fo, "[M::%s] number unitigs  : %lu\n", __func__, n_utg);
        fprintf(fo, "[M::%s] number syncmers : %lu\n", __func__, n_scm);
        fprintf(fo, "[M::%s] number arcs     : %lu\n", __func__, n_arc);
    }

    if (stats != 0) {
        stats[0] = n_scm;
        stats[1] = u_scm;
        stats[2] = n_utg;
        stats[3] = n_arc;
    }

    return;
}

static inline void add_ovl_count(kh_generic_t *h, int key)
{
    int absent;
    khint_t k;
    k = kh_generic_put(h, key, &absent);
    if (absent)
        kh_val(h, k) = 1;
    else
        ++kh_val(h, k);
    return;
}

static int calc_syncmer_overlap(sr_v *sr, syncmer_t *m1, uint64_t rc1, syncmer_t *m2, uint64_t rc2, void *hm)
{
    // more precisely the distance between adjacent syncmers
    // always m1 -> m2
    // rc1 and rc2 specify the strand of m1 and m2
    uint64_t *pos1, *pos2, n1, n2, r1, r2, i1, i2, c1, c2;
    uint32_t i, p1, p2, c;

    pos1 = m1->m_pos;
    pos2 = m2->m_pos;
    n1 = m1->k_cov;
    n2 = m2->k_cov;

    assert(n1>0 && n2>0);

    // compare each position pair to calculate overlap
    kh_generic_t *h;
    h = hm? (kh_generic_t *) hm : kh_generic_init();
    kh_generic_m_clear(h);
    for (p1 = p2 = 0; p1 < n1; ++p1) {
        r1 = pos1[p1]>>25;
        i1 = pos1[p1]>>1&MAX_RD_SCM;
        c1 = pos1[p1]&1;
        while (p2 < n2 && (r2 = pos2[p2]>>25) < r1) ++p2;
        if (r1 != r2) continue;
        for (i = p2; i < n2; ++i) {
            c = 0;

            r2 = pos2[i]>>25;
            if (r1 != r2) break;
            i2 = pos2[i]>>1&MAX_RD_SCM;
            c2 = pos2[i]&1;

            // rc1 == 0 && rc2 == 0
            //   c1 == 0 && c2 == 0 && i1+1 == i2
            //   c1 == 1 && c2 == 1 && i1 == i2+1
            //
            // rc1 == 0 && rc2 == 1
            //   c1 == 0 && c2 == 1 && i1+1 == i2
            //   c1 == 1 && c2 == 0 && i1 == i2+1
            //
            // rc1 == 1 && rc2 == 0
            //   c1 == 1 && c2 == 0 && i1+1 == i2
            //   c1 == 0 && c2 == 1 && i1 == i2+1
            //
            // rc1 == 1 && rc2 == 1
            //   c1 == 1 && c2 == 1 && i1+1 == i2
            //   c1 == 0 && c2 == 0 && i1 == i2+1

            if (i1 == i2+1 && c1 != rc1 && c2 != rc2) {
                add_ovl_count(h, (sr->a[r1].m_pos[i1]>>1)-(sr->a[r2].m_pos[i2]>>1));
                ++c;
            } else if (i1+1 == i2 && c1 == rc1 && c2 == rc2) {
                add_ovl_count(h, (sr->a[r2].m_pos[i2]>>1)-(sr->a[r1].m_pos[i1]>>1));
                ++c;
            }

            if (c > 1) {
                // this should never happen
                // for checking only
                fprintf(stderr, "[W::%s] multiple syncmer order found: %s ", __func__, sr->a[r1].sname);
                print_seq(&sr->a[r1], stderr);
                fprintf(stderr, "[W::%s] == %lu %lu %u\n", __func__, r1, i1, sr->a[r1].m_pos[i1]>>1);
                for (i = p2; i < n2; ++i) {
                    r2 = pos2[i]>>25;
                    if (r1 != r2) break;
                    i2 = pos2[i]>>1&MAX_RD_SCM;
                    fprintf(stderr, "[W::%s] ## %lu %lu %u\n", __func__, r2, i2, sr->a[r2].m_pos[i2]>>1);
                }
            }
        }
    }

    int movl, mcnt, smovl, smcnt, cnt;
    khint_t k;
    movl = mcnt = smovl = smcnt = 0;
    for (k = (khint_t) 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            cnt = kh_val(h, k);
            if (cnt > mcnt) {
                smcnt = mcnt;
                smovl = movl;
                mcnt = cnt;
                movl = kh_key(h, k);
            } else if (cnt > smcnt) {
                smcnt = cnt;
                smovl = kh_key(h, k);
            }
        }
    }

#ifdef DEBUG_KMER_OVERLAP
    fprintf(stderr, "[DEBUG_KMER_OVERLAP::%s] overlap count: %d[%d], %d[%d]\n", __func__, movl, mcnt, smovl, smcnt);
#endif

    // TODO overlap/gap length could be a indication of haplotypic differences
    // return smovl? 0 : movl;
    return movl;
}

static inline double utg_avg_cov(scg_t *scg, const scg_utg_t *s)
{
    if (s->del) return .0;

    uint64_t i, cov;
    cov = 0;
    for (i = 0; i < s->n; ++i)
        cov += scg->scm[s->a[i]>>1].k_cov;
    return (double) cov / s->n;
}

void scg_consensus(sr_v *sr, scg_t *scg, int w, int save_seq, FILE *fo)
{
    uint64_t i, v, t, z, n;
    int64_t l;
    double cov;
    asmg_t *utg_asmg;
    asmg_arc_t *a;
    scg_utg_t *s;
    kstring_t c_seq = {0, 0, 0};

    utg_asmg = scg->utg_asmg;
    if (fo) fprintf(fo, "H\tVN:Z:1.0\n");
    for (i = 0, n = utg_asmg->n_vtx; i < n; ++i) {
        s = &utg_asmg->vtx[i];
        if (s->del) continue;
        c_seq.l = 0;
        l = scg_unitig_consensus(sr, s, scg->scm, w, &c_seq);
        cov = s->cov? s->cov : utg_avg_cov(scg, s);
        s->cov = cov;
        s->len = l;
        if (save_seq) {
            if (s->seq) free(s->seq);
            MYMALLOC(s->seq, l+1);
            strncpy(s->seq, c_seq.s, l);
            s->seq[l] = 0;
        }
        if (fo) fprintf(fo, "S\tu%lu\t%.*s\tLN:i:%ld\tKC:i:%ld\tSC:f:%.3f\n", i, (int) l, c_seq.s, l, (int64_t) (l*cov), cov);
    }

    for (i = 0, n = utg_asmg->n_arc; i < n; ++i) {
        a = &utg_asmg->arc[i];
        if (a->del || a->comp) continue;
        s = &utg_asmg->vtx[a->v>>1];
        z = a->v&1;
        v = s->a[(s->n-1) * (!z)]^z;
        s = &utg_asmg->vtx[a->w>>1];
        z = a->w&1;
        t = s->a[(s->n-1) * z]^z;

        l = calc_syncmer_overlap(sr, &scg->scm[v>>1], v&1, &scg->scm[t>>1], t&1, 0);
        if (l < w) {
            c_seq.l = 0;
            l = scg_syncmer_consensus(sr, &scg->scm[v>>1], v&1, l, w, &c_seq);
        } else {
            l = 0;
        }
        // FIXME: consensus problem seql < lo
        l = MIN((uint64_t) l, utg_asmg->vtx[a->v>>1].len);
        l = MIN((uint64_t) l, utg_asmg->vtx[a->w>>1].len);
        a->lo = l;
        asmg_arc(utg_asmg, a->w^1, a->v^1)->lo = l;
        if (fo) {
            fprintf(fo, "L\tu%lu\t%c\tu%lu\t%c\t%ldM\tEC:i:%u\n", a->v>>1, "+-"[a->v&1], a->w>>1, "+-"[a->w&1], l, a->cov);
            fprintf(fo, "L\tu%lu\t%c\tu%lu\t%c\t%ldM\tEC:i:%u\n", a->w>>1, "-+"[a->w&1], a->v>>1, "-+"[a->v&1], l, a->cov);
        }
    }
    
    free(c_seq.s);
}

void scg_print_unitig_syncmer_list(scg_t *g, FILE *fo)
{
    uint64_t i, j, n, m, *a;
    asmg_t *utg;
    utg = g->utg_asmg;
    for (i = 0, n = utg->n_vtx; i < n; ++i) {
        if (utg->vtx[i].del)
            continue;
        a = utg->vtx[i].a;
        m = utg->vtx[i].n;
        fprintf(fo, "u%lu syncmer list:", i);
        for (j = 0; j < m; ++j)
            fprintf(fo, " %lu%c", a[j]>>1, "+-"[a[j]&1]);
        fprintf(fo, "\n");
    }
}


static void get_kmer_sequence(uint8_t *hoco_s, uint32_t pos, uint32_t rev, uint8_t *kmer_s, int w)
{
    int i, j;
    uint8_t t;
    uint32_t p;
    for (i = 0; i < w; ++i) {
        p = pos + i;
        kmer_s[i] = (hoco_s[p/4]>>(((p&3)^3)<<1))&3;
    }
    if (rev) {
        for (i = 0, j = w - 1; i < j; ++i, --j) {
            t = kmer_s[i];
            kmer_s[i] = (kmer_s[j]^3)&3;
            kmer_s[j] = (t^3)&3;
        }
        if (i == j) kmer_s[i] = (kmer_s[i]^3)&3;
    }
}

int64_t scg_syncmer_consensus(sr_v *sr, syncmer_t *scm, int rev, int64_t beg, int w, kstring_t *c_seq)
{
    assert(beg < w);

    uint32_t i, j, k, n_seq, rl;
    uint64_t *m_pos, p, r, l;
    uint8_t *ho_rl;
    uint32_t *ho_l_rl;
    int64_t bl;

    bl = beg<0? -beg : 0;
    while (beg < 0) {kputc_('N', c_seq), ++beg;}

    n_seq = scm->k_cov;
    m_pos = scm->m_pos;
    l = w - beg;
    bl += l;

    uint64_t tot_rl[l];
    uint8_t kmer_s[l];
#ifdef DEBUG_CONSENSUS
    uint8_t kmer_s1[l];
#endif

    MYBZERO(tot_rl, l);
    for (i = 0; i < n_seq; ++i) {
        const sr_t *s = &sr->a[m_pos[i]>>25];
        p = s->m_pos[m_pos[i]>>1&MAX_RD_SCM];
        r = (p&1)^rev;
        p >>= 1;
        if (!r) p += beg;

        if (i == 0) get_kmer_sequence(s->hoco_s, p, r, kmer_s, l);
#ifdef DEBUG_CONSENSUS
        get_kmer_sequence(s->hoco_s, p, r, kmer_s1, l);
        for (j = 0; j < l; ++j) assert(kmer_s[j] == kmer_s1[j]);
#endif
        ho_rl = s->ho_rl;
        ho_l_rl = s->ho_l_rl;
        k = 0;
        if (ho_l_rl)
            for (j = 0; j < p; ++j)
                if (ho_rl[j] == 255) ++k;
        if (r) {
            for (j = 0; j < l; ++j) {
                rl = ho_rl[p+j];
                if (rl == 255)
                    rl = ho_l_rl[k++];
                tot_rl[l-1-j] += rl;
            }
        } else {
            for (j = 0; j < l; ++j) {
                rl = ho_rl[p+j];
                if (rl == 255)
                    rl = ho_l_rl[k++];
                tot_rl[j] += rl;
            }
        }
    }

    for (i = 0; i < l; ++i) {
        kputc_(char_nt4_table[kmer_s[i]], c_seq);
        long b = lround((double) tot_rl[i] / n_seq);
        for (j = 0; j < b; ++j)
            kputc_(char_nt4_table[kmer_s[i]], c_seq);
        bl += b;
    }

    return bl;
}

int64_t scg_unitig_consensus(sr_v *sr, scg_utg_t *utg, syncmer_t *scm, int w, kstring_t *c_seq)
{
    uint64_t i;
    int64_t beg_pos, end_pos, l;
    kh_generic_t *h;

    h = kh_generic_init();
    uint64_t *a = utg->a;
    int64_t *pos;
    MYMALLOC(pos, utg->n);
    pos[0] = 0;
    for (i = 1; i < utg->n; ++i) {
        syncmer_t *m1 = &scm[a[i-1]>>1];
        syncmer_t *m2 = &scm[a[i]>>1];
        pos[i] = pos[i-1] + calc_syncmer_overlap(sr, m1, a[i-1]&1, m2, a[i]&1, h);
    }

#ifdef DEBUG_KMER_OVERLAP
    for (i = 0; i < utg->n; ++i)
        fprintf(stderr, "[DEBUG_KMER_OVERLAP::%s] overlap [%lu %lu]\n", __func__, i, pos[i]);
#endif

    beg_pos = end_pos = l = 0;
    for (i = 0; i < utg->n; ++i) {
        while (i+1 < utg->n && pos[i+1] <= end_pos) ++i;
        beg_pos = pos[i];
        l += scg_syncmer_consensus(sr, &scm[a[i]>>1], a[i]&1, end_pos - beg_pos, w, c_seq);
#ifdef DEBUG_KMER_OVERLAP
        fprintf(stderr, "[DEBUG_KMER_OVERLAP::%s] consensus start position [%lu %lu%c %ld %064lu%064lu]\n", __func__, i, a[i]>>1, "+-"[a[i]&1], end_pos-beg_pos, (uint64_t)(scm[a[i]>>1].h>>64), (uint64_t)(scm[a[i]>>1].h));
#endif
        end_pos = beg_pos + w;
    }
    free(pos);
    kh_generic_destroy(h);

    assert(l >= 0 && (uint64_t) l == c_seq->l);

    return l;
}

void process_mergeable_unitigs(scg_t *g)
{
    // postprocessing for mergeable utgs
    asmg_t *utg_asmg;

    // make utgs from the unitig graph
    utg_asmg = asmg_unitigging(g->utg_asmg);

    // update utg
    asmg_destroy(g->utg_asmg);
    g->utg_asmg = utg_asmg;

    scg_scm_utg_index(g);
}

static inline uint64_t scg_add_utg_p(asmg_t *g, scg_utg_t **p) {
    if (g->n_vtx == g->m_vtx)
        MYEXPAND(g->vtx, g->m_vtx);
    ++g->n_vtx;
    *p = &g->vtx[g->n_vtx-1];
    return g->n_vtx-1;
}

typedef struct { size_t n, m; uint64_t *a; } u64_v_t;

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

typedef struct {
    uint64_t vtx_new;
    u64_v_t arc_next;
} multi_arc_t;

// static uint64_t LONG_UTG_SCM_NUM = 20;
static int DOM_PATH_MIN_FC = 10;

int scg_multiplex(scg_t *g, scg_ra_v *ra_v, double r_thresh)
{
    uint64_t i, j, l0, l1, c0, c1, v, v1, w, av, aw, s, t, n, m, n_out, n_in, n_out1, n_in1, n_vtx, n_arc, max_l_id;
    uint64_t *a, stats0[4], stats1[4];
    double score, m_score, intpart;
    multi_arc_t *multi_arc;
    uint8_t *multi_vtx, multi;
    kvec_t(uint8_t) uniq;
    khint_t k;
    kh_dbl_t *tri_s; // spanning triplets
    ks_u128_t *arc_s;
    uint128_t k1;
    asmg_t *utg_g;
    asmg_vtx_t *utg, *vtx;
    asmg_arc_t *a_out, *a_in, *arc, *arc1;
    ra_frg_t *ra;
    int absent, updated;
    
    utg_g = g->utg_asmg;
    scg_stat(g, 0, stats0);

    // build hash table for spanning triplets represented by a pair of links
    tri_s = kh_dbl_init();
    kv_init(uniq);
    for (i = 0, n = ra_v->n; i < n; ++i) {
        m = ra_v->a[i].n;
        if (m < 3)
            continue;
        ra = ra_v->a[i].a;
        score = modf(ra_v->a[i].s, &intpart);
        if (score < DBL_EPSILON)
            score = 1.0;

        if (m > uniq.m)
            kv_resize(uint8_t, uniq, m);
        if (score < .99) { // score < 1.0
            // read was not uniquely mapped
            MYBZERO(uniq.a, m);
            for (j = 0; j < m; ++j) {
                a = utg_g->vtx[ra[j].uid>>1].a;
                for (s = ra[j].u_beg; s <= ra[j].u_end; ++s) {
                    if (scm_utg_n(g, a[s]>>1) == 1) {
                        // find a unique syncmer in the alignment
                        uniq.a[j] = 1;
                        break;
                    }
                }
            }
        } else {
            // no need to check one by one
            MYBONE(uniq.a, m);
        }
    
        arc = asmg_arc(utg_g, ra[0].uid, ra[1].uid);
        l0 = asmg_arc_id(*arc);
        c0 = asmg_comp_arc_id(arc);
        for (j = 2; j < m; ++j) {
            arc = asmg_arc(utg_g, ra[j-1].uid, ra[j].uid);
            l1 = asmg_arc_id(*arc);
            c1 = asmg_comp_arc_id(arc);
            if (uniq.a[j-2] && uniq.a[j-1] && uniq.a[j]) {
                // all three utgs are uniquely mapped to
                k = kh_dbl_put(tri_s, (uint128_t) l0 << 64 | l1, &absent);
                if (absent) {
                    kh_val(tri_s, k) = score;
                    k = kh_dbl_put(tri_s, (uint128_t) c1 << 64 | c0, &absent);
                    kh_val(tri_s, k) = score;
                } else {
                    kh_val(tri_s, k) += score;
                    k = kh_dbl_put(tri_s, (uint128_t) c1 << 64 | c0, &absent);
                    kh_val(tri_s, k) += score;
                }
            }
            l0 = l1;
            c0 = c1;
        }
    }

    // at most (max_l_id * 2 + 2) arcs
    max_l_id = asmg_max_link_id(utg_g);
    n_arc = utg_g->n_arc;
    n_vtx = utg_g->n_vtx;
    MYCALLOC(multi_arc, max_l_id * 2 + 2);
    MYCALLOC(multi_vtx, n_vtx);
    
    for (i = 0, n = max_l_id * 2 + 2; i < n; ++i)
        multi_arc[i].vtx_new = UINT64_MAX;

    // multiplex threadables
    for (i = 0; i < n_vtx; ++i) {
        if (utg_g->vtx[i].del)
            continue;
        
        v1 = i << 1;
        n_in = asmg_arc_n(utg_g, v1^1);
        n_in1 = asmg_arc_n1(utg_g, v1^1);
        a_in = asmg_arc_a(utg_g, v1^1);
        n_out = asmg_arc_n(utg_g, v1);
        n_out1 = asmg_arc_n1(utg_g, v1);
        a_out = asmg_arc_a(utg_g, v1);

        // singletons
        if (n_in1 == 0 && n_out1 == 0) {
            multi_vtx[i] = 2;
            continue;
        }
       
        m_score = 0;
        for (s = 0; s < n_in; ++s) {
            if (a_in[s].del)
                continue;
            // v0 -> v1 arc link_id
            l0 = asmg_comp_arc_id(&a_in[s]);

            for (t = 0; t < n_out; ++t) {
                if (a_out[t].del)
                    continue;
                // v1 -> v2 arc link_id
                l1 = asmg_arc_id(a_out[t]);

                k = kh_dbl_get(tri_s, (uint128_t) l0 << 64 | l1);
                score = k < kh_end(tri_s)? kh_val(tri_s, k) : .0;

                m_score = MAX(m_score, score);
            
                if (score >= r_thresh) {
                    // v0->v1->v2 is potentially multiplexable
                    // mark l0 l1
                    kv_push(uint64_t, multi_arc[l0].arc_next, a_out[t].w);
                    kv_push(uint64_t, multi_arc[l1^1].arc_next, a_in[s].w);
                }
            }
        }

        // mark as unmultiplexable if exist orphan arcs
        multi = 1;
        for (s = 0; s < n_in; ++s) {
            if (a_in[s].del)
                continue;
            // v0 -> v1 arc link_id
            l0 = asmg_comp_arc_id(&a_in[s]);
            if (multi_arc[l0].arc_next.n == 0) {
                multi = 0;
                break;
            }
        }

        if (multi == 1) {
            for (t = 0; t < n_out; ++t) {
                if (a_out[t].del)
                    continue;
                // v1 -> v2 arc link_id
                l1 = asmg_arc_id(a_out[t]);
                if (multi_arc[l1^1].arc_next.n == 0) {
                    multi = 0;
                    break;
                }
            }
        }

        // mark as multiplexable if exists a dominating arc link
        if (multi == 0)
            // check if the max_score is dominating
            multi = m_score >= r_thresh * DOM_PATH_MIN_FC;
        
        // retain all arc pairs if unmultiplexable
        if (multi == 0) {
            // clean arc pairs
            for (s = 0; s < n_in; ++s) {
                if (a_in[s].del)
                    continue;
                // v0 -> v1 arc link_id
                l0 = asmg_comp_arc_id(&a_in[s]);
                multi_arc[l0].arc_next.n = 0;
            }

            for (t = 0; t < n_out; ++t) {
                if (a_out[t].del)
                    continue;
                // v1 -> v2 arc link_id
                l1 = asmg_arc_id(a_out[t]);
                multi_arc[l1^1].arc_next.n = 0;
            }
                
            // mark all arc pairs
            for (s = 0; s < n_in; ++s) {
                if (a_in[s].del)
                    continue;
                // v0 -> v1 arc link_id
                l0 = asmg_comp_arc_id(&a_in[s]);

                for (t = 0; t < n_out; ++t) {
                    if (a_out[t].del)
                        continue;
                    // v1 -> v2 arc link_id
                    l1 = asmg_arc_id(a_out[t]);
                        
                    kv_push(uint64_t, multi_arc[l0].arc_next, a_out[t].w);
                    kv_push(uint64_t, multi_arc[l1^1].arc_next, a_in[s].w);
                }
            }
        }

        multi_vtx[i] = multi;
    }

    kh_dbl_destroy(tri_s);
    
    updated = 0;
    for (i = 0; i < n_arc; ++i) {
        arc = &utg_g->arc[i];
        if (arc->del || arc->comp)
            continue;
        if (multi_vtx[arc->v>>1] != 1 &&
                multi_vtx[arc->w>>1] != 1)
            continue;
        l0 = asmg_arc_id(*arc);
        v = scg_add_utg_p(utg_g, &utg);
        // record vtx id for the old arc
        multi_arc[l0].vtx_new = v << 1;
        multi_arc[l0^1].vtx_new = v << 1 | 1;
        // update a compound vtx
        u64_v_t v_vec = {0, 0, 0};
        vtx = &utg_g->vtx[arc->v>>1];
        vec_add(&v_vec, vtx->a, vtx->n, !!(arc->v&1));
        v_vec.n -= arc->lo;
        vtx = &utg_g->vtx[arc->w>>1];
        vec_add(&v_vec, vtx->a, vtx->n, !!(arc->w&1));
        MYREALLOC(v_vec.a, v_vec.n);
        utg->n = v_vec.n;
        utg->a = v_vec.a;
        utg->seq = 0;
        utg->len = 0;
        utg->cov = 0;
        utg->del = 0;
        utg->circ = 0;
        ++updated;
#ifdef DEBUG_UTG_MULTIPLEX
        fprintf(stderr, "[DEBUG_UTG_MULTIPLEX::%s] add vtx [u%lu]: u%lu%c u%lu%c\n", __func__,
                v,
                arc->v>>1, "+-"[arc->v&1],
                arc->w>>1, "+-"[arc->w&1]);
#endif
    }

    if (updated == 0)
        goto no_updates;

    // make new arcs
    arc_s = ks_u128_init();
    for (i = 0; i < n_arc; ++i) {
        arc = &utg_g->arc[i];
        if (arc->del)
            continue;

        av = arc->v;
        aw = arc->w;
        l0 = asmg_arc_id(*arc);
        v = multi_arc[l0].vtx_new;
        s = v == UINT64_MAX? aw : v;

        a = multi_arc[l0].arc_next.a;
        for (j = 0, m = multi_arc[l0].arc_next.n; j < m; ++j) {
            arc1 = asmg_arc(utg_g, aw, a[j]);
            l1 = asmg_arc_id(*arc1);
            w = multi_arc[l1].vtx_new;
            t = w == UINT64_MAX? aw : w;

            if (v != UINT64_MAX || w != UINT64_MAX) {
                // check if arc has been added
                k1 = (uint128_t) s << 64 | t;
                if (ks_u128_get(arc_s, k1) < kh_end(arc_s))
                    continue;
                ks_u128_put(arc_s, k1, &absent);

                // add arc v1 -> w
                // fix comp flag in asmg_finalize
                asmg_arc_add(utg_g, 
                        s,
                        t,
                        utg_g->vtx[aw>>1].n, 
                        UINT64_MAX, 0, 0);
#ifdef DEBUG_UTG_MULTIPLEX
                fprintf(stderr, "[DEBUG_UTG_MULTIPLEX::%s] add arc: u%lu%c u%lu%c u%lu%c %s[u%lu%c]->%s[u%lu%c]\n", __func__, 
                        av>>1, "+-"[av&1],
                        aw>>1, "+-"[aw&1],
                        a[j]>>1, "+-"[a[j]&1],
                        v == UINT64_MAX? "OLD" : "NEW",
                        s>>1, "+-"[s&1],
                        w == UINT64_MAX? "OLD" : "NEW",
                        t>>1, "+-"[t&1]);
#endif
            }
        }
    }
    ks_u128_destroy(arc_s);

    // mark arcs to delete
    for (i = 0; i < n_arc; ++i) {
        arc = &utg_g->arc[i];
        if (arc->del)
            continue;
        l0 = asmg_arc_id(*arc);
        if (multi_arc[l0].vtx_new != UINT64_MAX) {
            arc->del = 1;
#ifdef DEBUG_UTG_MULTIPLEX
            fprintf(stderr, "[DEBUG_UTG_MULTIPLEX::%s] del arc: u%lu%c u%lu%c\n", __func__,
                    arc->v>>1, "+-"[arc->v&1],
                    arc->w>>1, "+-"[arc->w&1]);
#endif
        }
    }

    // mark vertices to delete
    for (i = 0; i < n_vtx; ++i) {
        if (utg_g->vtx[i].del || multi_vtx[i] == 2)
            continue;

        v1 = i << 1;
        n_in1 = asmg_arc_n1(utg_g, v1^1);
        n_out1 = asmg_arc_n1(utg_g, v1);
        if (n_in1 == 0 && n_out1 == 0) {
            utg_g->vtx[i].del = 1;
#ifdef DEBUG_UTG_MULTIPLEX
            fprintf(stderr, "[DEBUG_UTG_MULTIPLEX::%s] del vtx: u%lu\n", __func__, i);
#endif
        }
    }

#ifdef DEBUG_UTG_MULTIPLEX
    for (i = 0, n = utg_g->n_arc; i < n; ++i) {
        arc = &utg_g->arc[i];
        fprintf(stderr, "[DEBUG_UTG_MULTIPLEX::%s] arc list: u%lu%c u%lu%c%s\n", __func__,
                arc->v>>1, "+-"[arc->v&1],
                arc->w>>1, "+-"[arc->w&1],
                arc->del? " [DEL]" : "");
    }
#endif

    // sort and index arcs
    asmg_finalize(utg_g, 1);

    // postprocessing for mergeable utgs after multiplexing
    process_mergeable_unitigs(g);

no_updates:
    
    for (i = 0, n = max_l_id * 2; i < n; ++i)
        kv_destroy(multi_arc[i].arc_next);
    free(multi_arc);
    free(multi_vtx);

    updated = 0;
    scg_stat(g, 0, stats1);
    for (i = 0; i < 4; ++i) {
        if (stats0[i] != stats1[i]) {
            updated = 1;
            break;
        }
    }

    return updated;
}

KDQ_INIT(uint64_t)

void scg_demultiplex(scg_t *g)
{
    uint64_t i, j, k, n, m, v, w, nv, pv, *a;
    int8_t *flag;
    kdq_t(uint64_t) *q;
    kvec_t(uint64_t) sub_g;
    asmg_arc_t *av;
    asmg_vtx_t *vtx;
    asmg_t *utg_g, *de_g;
    khint_t k1;
    kh_u64_t *h_scm;
    ks_u128_t *arc_s;
    int absent;

    utg_g = g->utg_asmg;
    q = kdq_init(uint64_t, 0);
    MYCALLOC(flag, utg_g->n_vtx * 2);
    MYCALLOC(de_g, 1);
    h_scm = kh_u64_init();
    arc_s = ks_u128_init();
    kv_init(sub_g);
    for (i = 0, n = utg_g->n_vtx * 2; i < n; ++i) {
        if (flag[i] || utg_g->vtx[i>>1].del)
            continue;

        // get a sub graph
        sub_g.n = 0;
        kdq_push(uint64_t, q, i);
        kdq_push(uint64_t, q, i^1);
        while (kdq_size(q) > 0) {
            v = *kdq_shift(uint64_t, q);
            if (flag[v]) continue;
            if (v & 1) kv_push(uint64_t, sub_g, v>>1);
            nv = asmg_arc_n(utg_g, v);
            av = asmg_arc_a(utg_g, v);
            for (j = 0; j < nv; ++j) {
                if (av[j].del) continue;
                if (flag[av[j].w] == 0)
                    kdq_push(uint64_t, q, av[j].w);
                if (flag[av[j].w^1] == 0)
                    kdq_push(uint64_t, q, av[j].w^1);
            }
            flag[v] = 1;
        }

#ifdef DEBUG_UTG_DEMULTIPLEX
        fprintf(stderr, "[DEBUG_UTG_DEMULTIPLEX::%s] processing subgraph:", __func__);
        for (j = 0; j < sub_g.n; ++j) fprintf(stderr, " u%lu", sub_g.a[j]);
        fprintf(stderr, "\n");
#endif
        // processing individual unitigs in sub_g
        for (j = 0; j < sub_g.n; ++j) {
#ifdef DEBUG_UTG_DEMULTIPLEX
            fprintf(stderr, "[DEBUG_UTG_DEMULTIPLEX::%s] process utg: u%lu\n", __func__, sub_g.a[j]);
#endif
            a = utg_g->vtx[sub_g.a[j]].a;
            m = utg_g->vtx[sub_g.a[j]].n;
            nv = 0;
            for (k = 0; k < m; ++k) {
                pv = nv;
                v = a[k] >> 1;
                k1 = kh_u64_put(h_scm, v, &absent);
                if (absent) {
                    // add v
                    nv = scg_add_utg_p(de_g, &vtx);
                    MYBZERO(vtx, 1);
                    MYMALLOC(vtx->a, 1);
                    vtx->n = 1;
                    vtx->a[0] = v << 1;
                    kh_val(h_scm, k1) = nv;
#ifdef DEBUG_UTG_DEMULTIPLEX
                    fprintf(stderr, "[DEBUG_UTG_DEMULTIPLEX::%s] add vtx: %lu [%lu%c]\n", 
                            __func__, nv, v, "+-"[a[k]&1]);
#endif
                } else {
                    // v has been added
                    nv = kh_val(h_scm, k1);
#ifdef DEBUG_UTG_DEMULTIPLEX
                    fprintf(stderr, "[DEBUG_UTG_DEMULTIPLEX::%s] existed vtx: %lu [%lu%c]\n", 
                            __func__, nv, v, "+-"[a[k]&1]);
#endif
                }
                if (k > 0) {
                    // add arc pv->nv
                    v = pv << 1 | (a[k-1]&1);
                    w = nv << 1 | (a[k]&1);
                    k1 = ks_u128_get(arc_s, (uint128_t) v<<64 | w);
                    if (k1 >= kh_end(arc_s)) {
                        // arc does not exist
                        asmg_arc_add2(de_g, v, w, 0, 0, 0, 0);
                        ks_u128_put(arc_s, (uint128_t) v<<64 | w, &absent);
                        ks_u128_put(arc_s, (uint128_t) (w^1)<<64 | (v^1), &absent);
#ifdef DEBUG_UTG_DEMULTIPLEX
                        fprintf(stderr, "[DEBUG_UTG_DEMULTIPLEX::%s] add arc: %lu%c -> %lu%c\n", 
                                __func__, v>>1, "+-"[v&1], w>>1, "+-"[w&1]);
#endif
                    }
#ifdef DEBUG_UTG_DEMULTIPLEX
                    else {
                        fprintf(stderr, "[DEBUG_UTG_DEMULTIPLEX::%s] existed arc: %lu%c -> %lu%c\n",
                                __func__, v>>1, "+-"[v&1], w>>1, "+-"[w&1]);
                    }
#endif
                }
            }
        }
        
        // add arcs between unitigs
        m = sub_g.n * 2;
        for (j = 0; j < m; ++j) {
            v = sub_g.a[j>>1];
            a = utg_g->vtx[v].a;
            pv = (j&1)? a[0]^1 : a[utg_g->vtx[v].n-1];
            k1 = kh_u64_get(h_scm, pv>>1);
            pv = kh_val(h_scm, k1)<<1 | (pv&1);
            for (k = 0; k < m; ++k) {
                w = sub_g.a[k>>1];
                av = asmg_arc1(utg_g, v<<1|(j&1), w<<1|(k&1));
                if (av == 0 || av->lo > 0)
                    continue;
                a = utg_g->vtx[w].a;
                nv = (k&1)? a[utg_g->vtx[w].n-1]^1 : a[0];
                k1 = kh_u64_get(h_scm, nv>>1);
                nv = kh_val(h_scm, k1)<<1 | (nv&1);
                if (ks_u128_get(arc_s, (uint128_t) pv<<64 | nv) >= kh_end(arc_s)) {
                    // add vtx pv->nv
                    asmg_arc_add(de_g, pv, nv, 0, 0, 0, 0);
                    ks_u128_put(arc_s, (uint128_t) pv<<64 | nv, &absent);
#ifdef DEBUG_UTG_DEMULTIPLEX
                    fprintf(stderr, "[DEBUG_UTG_DEMULTIPLEX::%s] add inter-utg arc: %lu%c -> %lu%c\n",
                            __func__, pv>>1, "+-"[pv&1], nv>>1, "+-"[nv&1]);
#endif
                }
            }
        }
        
        kh_u64_m_clear(h_scm);
        ks_u128_s_clear(arc_s);
    }
    
    kdq_destroy(uint64_t, q);
    kv_destroy(sub_g);
    kh_u64_destroy(h_scm);
    free(flag);
    
    // sort arc and build index
    asmg_finalize(de_g, 1);
    
    // update utg graph
    asmg_destroy(utg_g);
    g->utg_asmg = de_g;

    // postprocessing for mergeable utgs after purging
    process_mergeable_unitigs(g);
}



