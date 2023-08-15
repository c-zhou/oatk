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

#include "syncasm.h"

#undef DEBUG_KMER_OVERLAP
#undef DEBUG_CONSENSUS
#undef DEBUG_SCM_UTG_INDEX
#undef DEBUG_UTG_MULTIPLEX
#undef DEBUG_UTG_DEMULTIPLEX
#undef DEBUG_UTG_COVERAGE
#undef DEBUG_LCS_BLOCK
#undef DEBUG_LCS_MALIGN
#undef DEBUG_UTG_RA_COV

static kh_inline khint_t kh_hash_uint128(uint128_t key)
{
    khint_t k1 = kh_hash_uint64((khint64_t) key);
    khint_t k2 = kh_hash_uint64((khint64_t) (key >> 64));
    return kh_hash_uint64((khint64_t) ((uint64_t) k1 << 32 | k2));
}

KHASHL_MAP_INIT(KH_LOCAL, kh_generic_t, kh_generic, int, int, kh_hash_dummy, kh_eq_generic)
KHASHL_MAP_INIT(KH_LOCAL, kh_dbl_t, kh_dbl, uint128_t, double, kh_hash_uint128, kh_eq_generic)
KHASHL_MAP_INIT(KH_LOCAL, kh_u64_t, kh_u64, uint64_t, uint64_t, kh_hash_uint64, kh_eq_generic)
KHASHL_MAP_INIT(KH_LOCAL, kh_u128_t, kh_u128, uint128_t, uint64_t, kh_hash_uint128, kh_eq_generic)
KHASHL_SET_INIT(KH_LOCAL, ks_u128_t, ks_u128, uint128_t, kh_hash_uint128, kh_eq_generic)

scg_t *scg_clone(scg_t *g)
{
    if (!g) return 0;
    scg_t *g1;
    MYCALLOC(g1, 1);
    g1->n_scm = g->n_scm;
    g1->m_scm = g->m_scm;
    g1->scm = g->scm;
    g1->h_scm = g->h_scm;
    g1->clone = true;
    return g1;
}

static void scg_destroy_scm(scg_t *g)
{
    if (!g) return;
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
    if (!g->clone) scg_destroy_scm(g);
    if (g->h_scm && !g->clone) kh_u128_destroy((kh_u128_t *) g->h_scm);
    if (g->utg_asmg) asmg_destroy(g->utg_asmg);
    if (g->scm_u) free(g->scm_u);
    if (g->idx_u) free(g->idx_u);
    free(g);
}

void scg_meta_clean(scg_meta_t *meta)
{
    if (!meta) return;
    scg_destroy(meta->scg);
    sr_v_destroy(meta->sr);
    scg_ra_v_destroy(meta->ra);
}

void scg_meta_destroy(scg_meta_t *meta)
{
    if (!meta) return;
    scg_meta_clean(meta);
    free(meta);
}

static int uint128_cmpfunc(const void *a, const void *b)
{
    return (*(uint128_t *) a > *(uint128_t *) b) - (*(uint128_t *) a < *(uint128_t *) b);
}

static int uint64_cmpfunc(const void *a, const void *b)
{
    return (*(uint64_t *) a > *(uint64_t *) b) - (*(uint64_t *) a < *(uint64_t *) b);
}

static int dbl_cmpfunc(const void *a, const void *b)
{
    return (*(double *) a > *(double *) b) - (*(double *) a < *(double *) b);
}

/***
static int uint32_cmpfunc(const void *a, const void *b)
{
    return (*(uint32_t *) a > *(uint32_t *) b) - (*(uint32_t *) a < *(uint32_t *) b);
}
**/

static void scg_scm_utg_index(scg_t *g)
{
    if (g->idx_u) {
        free(g->idx_u);
        g->idx_u = 0;
    }
    if (g->scm_u) {
        free(g->scm_u);
        g->scm_u = 0;
    }

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
    qsort(scm_u.a, scm_u.n, sizeof(uint128_t), uint128_cmpfunc);

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

uint64_t scg_get_scm_id(scg_t *g, uint128_t key)
{
    kh_u128_t *h = (kh_u128_t *) g->h_scm;
    khint_t k = kh_u128_get(h, key);
    if (k == kh_end(h))
        return UINT64_MAX;
    return kh_val(h, k);
}

static inline void add_a_cov(kh_u128_t *h, uint128_t v_p, uint64_t c)
{
    int absent;
    khint_t k;
    k = kh_u128_put(h, v_p, &absent);
    if (absent)
        kh_val(h, k) = 0;
    kh_val(h, k) += c;
    return;
}

static kh_u128_t *build_scm_hash_table(syncmer_t *scm, uint64_t n) {
    // build syncmer hash h128 to id map
    uint64_t i;
    int absent;
    khint_t k;
    kh_u128_t *h_scm;
    h_scm = kh_u128_init();
    for (i = 0; i < n; ++i) {
        k = kh_u128_put(h_scm, scm[i].h, &absent);
        kh_val(h_scm, k) = i;
    }
    return h_scm;
}

scg_t *make_syncmer_graph(sr_v *sr, uint32_t min_k_cov, double min_a_cov_f)
{
    scg_t *scg;
    MYCALLOC(scg, 1);

    scg->scm = collect_syncmer_from_reads(sr, &scg->n_scm);
    scg->m_scm = scg->n_scm;
    if (scg->n_scm == 0) {
        scg_destroy(scg);
        return 0;
    }
    scg->h_scm = build_scm_hash_table(scg->scm, scg->n_scm);

    size_t i, j;
    sr_t *s;
    // add utg vtx
    // each single scm is a utg
    asmg_vtx_t *vtx, *vtx1;
    MYCALLOC(vtx, scg->n_scm);
    for (i = 0; i < scg->n_scm; ++i) {
        scg->scm[i].del = scg->scm[i].k_cov < min_k_cov; // filter by k_cov
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

    // add arcs
    kvec_t(asmg_arc_t) arc;
    kv_init(arc);
    uint64_t v0, v1;
    uint32_t r0, r1, v_v;
    uint128_t v_p; // v0, v1 packed
    kh_u128_t *a_cov; // arc cov
    a_cov = kh_u128_init();
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

    khint_t k;
    for (k = (khint_t) 0; k < kh_end(a_cov); ++k) {
        if (kh_exist(a_cov, k)) {
            v_p = kh_key(a_cov, k);
            v_v = kh_val(a_cov, k);
            v0 = (uint64_t) (v_p >> 64);
            v1 = (uint64_t) v_p;
            if (v_v < min_a_cov_f * MIN(scg->scm[v0>>1].k_cov, scg->scm[v1>>1].k_cov)
                    || vtx[v0>>1].del || vtx[v1>>1].del)
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
    kh_u128_destroy(a_cov);

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
    kh_u128_t *a_cov; // arc cov
    asmg_t *g;
    asmg_arc_t *a;
    asmg_vtx_t *v;
    sr_t *s;

    g = scg->utg_asmg;
    a_cov = kh_u128_init();
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

            k = kh_u128_get(a_cov, (uint128_t) v0 << 64 | v1);
            if (k < kh_end(a_cov)) {
                kh_val(a_cov, k) += 1;
                if ((v1^1) != v0 || (v0^1) != v1) {
                    k = kh_u128_get(a_cov, (uint128_t) (v1^1) << 64 | (v0^1));
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
        k = kh_u128_get(a_cov, (uint128_t) v0 << 64 | v1);
        a->cov = kh_val(a_cov, k);
    }
    kh_u128_destroy(a_cov);
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
    uint64_t *pos1, *pos2, n1, n2, r1, r2, i1, i2, l1, l2, c1, c2;
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
        r1 = pos1[p1]>>32;
        i1 = pos1[p1]>>1&MAX_RD_SCM;
        l1 = sr->a[r1].m_pos[i1]>>1;
        if (l1 == MAX_RD_LEN) continue; // this mer is error-corrected
        c1 = pos1[p1]&1;
        while (p2 < n2 && (r2 = pos2[p2]>>32) < r1) ++p2;
        if (r1 != r2) continue;
        for (i = p2; i < n2; ++i) {
            c = 0;

            r2 = pos2[i]>>32;
            if (r1 != r2) break;
            i2 = pos2[i]>>1&MAX_RD_SCM;
            l2 = sr->a[r2].m_pos[i2]>>1;
            if (l2 == MAX_RD_LEN) continue; // this mer is error-corrected
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
                add_ovl_count(h, l1 - l2);
                ++c;
            } else if (i1+1 == i2 && c1 == rc1 && c2 == rc2) {
                add_ovl_count(h, l2 - l1);
                ++c;
            }

            if (c > 1) {
                // this should never happen
                // for checking only
                fprintf(stderr, "[W::%s] multiple syncmer order found: %s ", __func__, sr->a[r1].sname);
                print_hoco_seq(&sr->a[r1], stderr);
                fprintf(stderr, "[W::%s] == %lu %lu %lu\n", __func__, r1, i1, l1);
                for (i = p2; i < n2; ++i) {
                    r2 = pos2[i]>>32;
                    if (r1 != r2) break;
                    i2 = pos2[i]>>1&MAX_RD_SCM;
                    l2 = sr->a[r2].m_pos[i2]>>1;
                    fprintf(stderr, "[W::%s] ## %lu %lu %lu\n", __func__, r2, i2, l2);
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
    
    if (!hm) kh_generic_destroy(h);

    // TODO overlap/gap length could be a indication of haplotypic differences
    // return smovl? 0 : movl;
    return movl;
}

static double quantile(double *a, int n, double q, int sorted)
{
    if (q < 0 || q > 1.0 || n <= 0) {
        fprintf(stderr, "[E::%s] numeric error", __func__);
        exit(EXIT_FAILURE);
    }
    if (n == 1) return a[0];
    if (!sorted) qsort(a, n, sizeof(double), dbl_cmpfunc);
    double fractpart, intpart;
    fractpart = modf(q * (n - 1), &intpart);
    int i = lround(intpart);
    if (i == n - 1) return a[i];
    return a[i] + (a[i + 1] - a[i]) * fractpart;
}

static double average_IQR(double *a, int n, int sorted)
{
    if (n == 0) return 0.;

    double Q1, Q3, IQR, s;
    int i, n0;
    if (!sorted) qsort(a, n, sizeof(double), dbl_cmpfunc);
    Q1 = quantile(a, n, 0.25, 1);
    Q3 = quantile(a, n, 0.75, 1);
    IQR = Q3 - Q1;
    Q1 -= 1.5 * IQR;
    Q3 += 1.5 * IQR;
    n0 = 0;
    s = 0.;
    for (i = 0; i < n; ++i) {
        if (a[i] >= Q1 && a[i] <= Q3) {
            ++n0;
            s += a[i];
        }
    }
#ifdef DEBUG_UTG_RA_COV
    fprintf(stderr, "[DEBUG_UTG_RA_COV::%s] %d single copy syncmers:",
            __func__, n);
    for (i = 0; i < n; ++i) fprintf(stderr, " %.0f", a[i]);
    fputc('\n', stderr);
    fprintf(stderr, "[DEBUG_UTG_RA_COV::%s] %d outliers out of range [%.3f %.3f] removed\n",
            __func__, n - n0, Q1, Q3);
#endif
    return n0? s / n0 : 0.;
}

static inline double utg_avg_cov(scg_t *scg, const scg_utg_t *s)
{
    if (s->del) return .0;

    uint64_t i, u;
    double *cov, avg_cov;

    MYCALLOC(cov, s->n);
    // try to use single copy syncmers
    for (i = 0; i < s->n; ++i) {
        u = s->a[i]>>1;
        if (scm_utg_n(scg, u) == 1)
            cov[i] = scg->scm[u].k_cov;
    }
    qsort(cov, s->n, sizeof(double), dbl_cmpfunc);
    i = 0;
    while (i < s->n && cov[i] < DBL_EPSILON) ++i;

    if (i == s->n) {
        // using all syncmers
        MYBZERO(cov, s->n);
        for (i = 0; i < s->n; ++i)
            cov[i] = scg->scm[s->a[i]>>1].k_cov;
        qsort(cov, s->n, sizeof(double), dbl_cmpfunc);
        i = 0;
    }
    // calculate average using [Q1-IQR, Q3+IQR]
    avg_cov = average_IQR(&cov[i], s->n - i, 1);
    
    free(cov);
    
    return avg_cov;
}

/***
static inline double utg_med_cov(scg_t *scg, const scg_utg_t *s)
{
    if (s->del) return .0;
    uint64_t i;
    uint32_t *cov;
    MYMALLOC(cov, s->n);
    for (i = 0; i < s->n; ++i)
        cov[i] = scg->scm[s->a[i]>>1].k_cov;
    qsort(cov, s->n, sizeof(uint32_t), uint32_cmpfunc);
    double med = (cov[(s->n-1)/2] + cov[s->n/2]) / 2.0;
    free(cov);
    return med;
}
**/

void scg_update_utg_cov(scg_t *scg)
{
    uint64_t i, n;
    scg_utg_t *s;
    asmg_t *utg_asmg = scg->utg_asmg;
    for (i = 0, n = utg_asmg->n_vtx; i < n; ++i) {
        s = &utg_asmg->vtx[i];
        s->cov = utg_avg_cov(scg, s);
    }
}

void scg_consensus(sr_v *sr, scg_t *scg, int w, int hoco_seq, int save_seq, FILE *fo)
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
        l = scg_unitig_consensus(sr, s, scg->scm, w, &c_seq, hoco_seq);
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

#ifdef DEBUG_UTG_COVERAGE
        uint64_t j;
        fprintf(stderr, "[DEBUG_UTG_COVERAGE::%s] u%lu [N=%lu]:", __func__, i, s->n);
        for (j = 0; j < s->n; ++j) fprintf(stderr, " s%lu%c", s->a[j]>>1, "+-"[s->a[j]&1]);
        fputc('\n', stderr);
        fprintf(stderr, "[DEBUG_UTG_COVERAGE::%s] u%lu [N=%lu]:", __func__, i, s->n);
        for (j = 0; j < s->n; ++j) fprintf(stderr, " %u", scg->scm[s->a[j]>>1].k_cov);
        fputc('\n', stderr);
#endif
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
            l = scg_syncmer_consensus(sr, &scg->scm[v>>1], v&1, l, w, &c_seq, hoco_seq);
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

int64_t scg_syncmer_consensus(sr_v *sr, syncmer_t *scm, int rev, int64_t beg, int w, kstring_t *c_seq, int hoco_seq)
{
    assert(beg < w);
    
    uint32_t i, n_seq;
    uint64_t *m_pos, p, r, l;
    int64_t bl;
    sr_t *s;

    bl = beg<0? -beg : 0;
    while (beg < 0) {kputc_('N', c_seq), ++beg;}

    n_seq = scm->k_cov;
    m_pos = scm->m_pos;
    l = w - beg;
    bl += l;

    // get the kmer sequence
    uint8_t kmer_s[l];
    i = 0;
    while (1) {
        s = &sr->a[m_pos[i]>>32];
        p = s->m_pos[m_pos[i]>>1&MAX_RD_SCM];
        r = (p&1)^rev;
        p >>= 1;
        if (p != MAX_RD_LEN)
            break;
        ++i;
    }
    if (!r) p += beg;
    get_kmer_seq(s->hoco_s, p, l, r, kmer_s);

    if (hoco_seq) {
        // do not need to do run length consensus
        for (i = 0; i < l; ++i)
            kputc_(char_nt4_table[kmer_s[i]], c_seq);
        return bl;
    }

    uint32_t j, k, rl;
    uint8_t *ho_rl;
    uint32_t *ho_l_rl;

    uint64_t tot_rl[l];
#ifdef DEBUG_CONSENSUS
    uint8_t kmer_s1[l];
#endif

    MYBZERO(tot_rl, l);
    for (i = 0; i < n_seq; ++i) {
        s = &sr->a[m_pos[i]>>32];
        p = s->m_pos[m_pos[i]>>1&MAX_RD_SCM];
        r = (p&1)^rev;
        p >>= 1;
        if (p == MAX_RD_LEN)
            continue;
        if (!r) p += beg;

#ifdef DEBUG_CONSENSUS
        get_kmer_seq(s->hoco_s, p, l, r, kmer_s1);
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

int64_t scg_unitig_consensus(sr_v *sr, scg_utg_t *utg, syncmer_t *scm, int w, kstring_t *c_seq, int hoco_seq)
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
        l += scg_syncmer_consensus(sr, &scm[a[i]>>1], a[i]&1, end_pos - beg_pos, w, c_seq, hoco_seq);
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

typedef struct {
    uint64_t l0, l1, w0, w1;
    double score, s_score;
    int rank;
} p_info_t;

static int pinfo_s_cmpfunc(const void *a, const void *b)
{
    return (((p_info_t *) a)->s_score < ((p_info_t *) b)->s_score) - (((p_info_t *) a)->s_score > ((p_info_t *) b)->s_score);
}

int scg_multiplex(scg_t *g, scg_ra_v *ra_v, uint32_t max_n_scm, double min_n_r, double min_d_f)
{
    uint64_t i, j, l0, l1, c0, c1, v, v1, w, aw, s, t, n, m, n_out, n_in, n_out1, n_in1, n_vtx, n_arc, max_l_id;
    uint64_t *a, stats0[4], stats1[4];
    double score, s_score, m_score, intpart;
    multi_arc_t *multi_arc;
    uint8_t *multi_vtx, multi;
    kvec_t(uint8_t) uniq;
    kvec_t(p_info_t) p_info;
    khint_t k;
    kh_dbl_t *tri_s; // spanning triplets
    ks_u128_t *arc_s;
    uint128_t k1;
    asmg_t *utg_g;
    asmg_vtx_t *utg, *vtx;
    asmg_arc_t *a_out, *a_in, *arc, *arc1;
    ra_frg_t *ra;
    int absent, updated, max_rank, rank;
    
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
    kv_destroy(uniq);

    // at most (max_l_id * 2 + 2) arcs
    max_l_id = asmg_max_link_id(utg_g);
    n_arc = utg_g->n_arc;
    n_vtx = utg_g->n_vtx;
    MYCALLOC(multi_arc, max_l_id * 2 + 2);
    MYCALLOC(multi_vtx, n_vtx);
    kv_init(p_info);

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

        p_info.n = 0;
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

                if (score >= min_n_r) {
                    // v0->v1->v2 is potentially multiplexable
                    // mark l0 l1
                    // scaled score by vtx coverage
                    s_score = score / MAX(1, MIN(utg_g->vtx[a_in[s].w>>1].cov, utg_g->vtx[a_out[t].w>>1].cov));
                    p_info_t p = {l0, l1^1, a_out[t].w, a_in[s].w, score, s_score, INT_MAX};
                    kv_push(p_info_t, p_info, p);

#ifdef DEBUG_UTG_MULTIPLEX
                    fprintf(stderr, "[DEBUG_UTG_MULTIPLEX::%s] triplex score u%lu%c u%lu+ u%lu%c: %.3f\n", __func__,
                            a_in[s].w>>1, "-+"[a_in[s].w&1],
                            v1>>1,
                            a_out[t].w>>1, "+-"[a_out[t].w&1],
                            score);
#endif

                }
            }
        }

        multi = 0;
        if (p_info.n > 0) {
            // rank all paths
            // sort by s_score in descending order
            qsort(p_info.a, p_info.n, sizeof(p_info_t), pinfo_s_cmpfunc);
            p_info.a[0].rank = 0;
            max_rank = 0;
            for (s = 1; s < p_info.n; ++s) {
                // find lowest possible rank
                p_info_t *p1 = &p_info.a[s];
                rank = 0;
                for (t = 0; t < s; ++t) {
                    p_info_t *p0 = &p_info.a[t];
                    if (p0->l0 == p1->l0 || p0->l1 == p1->l1) {
                        // comparable
                        if (p0->s_score * min_d_f > p1->s_score)
                            ++rank;
                    }
                }
                p1->rank = rank;
                max_rank = MAX(rank, max_rank);
            }

#ifdef DEBUG_UTG_MULTIPLEX
            for (s = 0; s < p_info.n; ++s) {
                p_info_t *p = &p_info.a[s];
                fprintf(stderr, "[DEBUG_UTG_MULTIPLEX::%s] triplex score u%lu%c u%lu+ u%lu%c: %.3f %.3f %d\n", __func__,
                        p->w1>>1, "-+"[p->w1&1],
                        v1>>1,
                        p->w0>>1, "+-"[p->w0&1],
                        p->score, p->s_score,
                        p->rank);
            }
#endif

            // add paths to threadable
            // paths with the same rank should be added together
            m_score = 0;
            for (rank = 0; rank <= max_rank; ++rank) {
                for (t = 0; t < p_info.n; ++t) {
                    if (p_info.a[t].rank != rank)
                        continue;
                    p_info_t *p = &p_info.a[t];
                    kv_push(uint64_t, multi_arc[p->l0].arc_next, p->w0);
                    kv_push(uint64_t, multi_arc[p->l1].arc_next, p->w1);
                    m_score = MAX(m_score, p->score);
                }

#ifdef DEBUG_UTG_MULTIPLEX
                fprintf(stderr, "[DEBUG_UTG_MULTIPLEX::%s] paths ranking %d added\n", __func__, rank);
#endif
                // check if need to add more
                // mark as unmultiplexable if exist orphan arcs
                uint8_t multi1 = 1;
                for (s = 0; s < n_in; ++s) {
                    if (a_in[s].del)
                        continue;
                    // v0 -> v1 arc link_id
                    l0 = asmg_comp_arc_id(&a_in[s]);
                    if (multi_arc[l0].arc_next.n == 0) {
                        multi1 = 0;
                        break;
                    }
                }

                if (multi1 == 1) {
                    for (t = 0; t < n_out; ++t) {
                        if (a_out[t].del)
                            continue;
                        // v1 -> v2 arc link_id
                        l1 = asmg_arc_id(a_out[t]);
                        if (multi_arc[l1^1].arc_next.n == 0) {
                            multi1 = 0;
                            break;
                        }
                    }
                }

                multi = multi1;
                if (multi1 == 1) break;
            }
            
            // mark as multiplexable if exists a dominating arc link
            if (multi == 0)
                // check if the max_score is dominating
                multi = m_score * min_d_f >= min_n_r;
        }

        /***
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
            
                if (score >= min_n_r) {
                    // v0->v1->v2 is potentially multiplexable
                    // mark l0 l1
                    kv_push(uint64_t, multi_arc[l0].arc_next, a_out[t].w);
                    kv_push(uint64_t, multi_arc[l1^1].arc_next, a_in[s].w);
#ifdef DEBUG_UTG_MULTIPLEX
                    fprintf(stderr, "[DEBUG_UTG_MULTIPLEX::%s] triplex score u%lu%c u%lu+ u%lu%c: %.3f\n", __func__,
                            a_in[s].w>>1, "-+"[a_in[s].w&1],
                            v1>>1,
                            a_out[t].w>>1, "+-"[a_out[t].w&1],
                            score);
#endif
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
            multi = m_score * min_d_f >= min_n_r;
        **/

        // do not do multiplexing for large utgs
        // the read spanning might be not confident
        if (utg_g->vtx[i].n > max_n_scm)
            multi = 0;
        
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
    kv_destroy(p_info);

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

#ifdef DEBUG_UTG_MULTIPLEX
    uint64_t av;
#endif
    // make new arcs
    arc_s = ks_u128_init();
    for (i = 0; i < n_arc; ++i) {
        arc = &utg_g->arc[i];
        if (arc->del)
            continue;
#ifdef DEBUG_UTG_MULTIPLEX
        av = arc->v;
#endif
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
    
    for (i = 0, n = max_l_id * 2 + 2; i < n; ++i)
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
    ks_u128_destroy(arc_s);
    free(flag);
    
    // sort arc and build index
    asmg_finalize(de_g, 1);
    
    // update utg graph
    asmg_destroy(utg_g);
    g->utg_asmg = de_g;

    // postprocessing for mergeable utgs after purging
    process_mergeable_unitigs(g);
}

// multiple alignment block
typedef struct {
    uint32_t a, b; // number alignments, blocks
    uint32_t *n; // matched syncmers in each block
    uint64_t *u; // unitig list, size a * b
} ma_t;

typedef struct { size_t n, m; uint64_t *a; } lcs_block_t;

static void lcs_backtrace(int **L, uint64_t *s_scm, int s_n, uint64_t *u_scm, int u_n, int offset, lcs_block_t *t)
{
    if (s_n == 0 || u_n == 0)
        return;
    if ((s_scm[s_n] >> 1) == (u_scm[u_n] >> 1)) {
        kv_push(uint64_t, *t, (uint64_t) (s_n - 1 + offset) << 32 | 1);
        return lcs_backtrace(L, s_scm, s_n - 1, u_scm, u_n - 1, offset, t);
    }
    if (L[s_n][u_n - 1] > L[s_n - 1][u_n])
        return lcs_backtrace(L, s_scm, s_n, u_scm, u_n - 1, offset, t);
    return lcs_backtrace(L, s_scm, s_n - 1, u_scm, u_n, offset, t);
}

static void lcs_block_merge(lcs_block_t *blocks)
{
    // merge blocks
    size_t k, p;
    if (blocks->n > 1) {
        for (p = 0, k = 1; k < blocks->n; ++k) {
            if ((blocks->a[p] >> 32) + (uint32_t) blocks->a[p] == (blocks->a[k] >> 32))
                blocks->a[p] += (uint32_t) blocks->a[k];
            else
                blocks->a[++p] = blocks->a[k];
        }
        blocks->n = p + 1;
    }
}

// TODO use max_d for banded LCS
static uint64_t *find_lcs(uint64_t *s_scm, int s_n, uint64_t *u_scm, int u_n, int offset, int max_d, int *n_block)
{
    int i, j, start, s_end, u_end;
    lcs_block_t blocks;

    kv_init(blocks);
    *n_block = 0;

    start = 0;
    s_end = s_n - 1;
    u_end = u_n - 1;

    // trim off the matching syncmers at the beginning and end
    while (start < s_n && start < u_n && (s_scm[start] >> 1) == (u_scm[start] >> 1)) { ++start; }
    while (start <= s_end && start <= u_end && (s_scm[s_end] >> 1) == (u_scm[u_end] >> 1)) { --s_end; --u_end; }

    if (start > 0)
        kv_push(uint64_t, blocks, (uint64_t) offset << 32 | start);

    s_scm = &s_scm[start];
    u_scm = &u_scm[start];
    s_end = s_end - start + 1;
    u_end = u_end - start + 1;

    // core LCS
    int *D, **L;
    MYCALLOC(D, (s_end + 1) * (u_end + 1));
    MYMALLOC(L, s_end + 1);
    L[0] = D;
    for (i = 1; i <= s_end; ++i) {
        L[i] = &L[i - 1][u_end + 1];
        for (j = 1; j <= u_end; ++j) {
            if ((s_scm[i - 1] >> 1) == (u_scm[j - 1] >> 1))
                L[i][j] = L[i - 1][j - 1] + 1;
            else
                L[i][j] = MAX(L[i - 1][j], L[i][j - 1]);
        }
    }

    // LCS backtrace
    uint32_t n_b = blocks.n;
    lcs_backtrace(L, s_scm - 1, s_end, u_scm - 1, u_end, offset + start, &blocks);
    free(D);
    free(L);

    // reverse backtrace array
    array_reverse(&blocks.a[n_b], blocks.n - n_b);

    if (start + s_end < s_n)
        kv_push(uint64_t, blocks, (uint64_t) (offset + start + s_end) << 32 | (s_n - start - s_end));

    // merge LCS blocks
    lcs_block_merge(&blocks);

#ifdef DEBUG_LCS_BLOCK
    size_t k;
    for (k = 0; k < blocks.n; ++k)
        fprintf(stderr, "[DEBUG_LCS_BLOCK::%s] LCS_MATCH_BLOCK_%lu: %lu %u\n",
                __func__, k, blocks.a[k]>>32, (uint32_t) blocks.a[k]);
#endif

    *n_block = blocks.n;

    return blocks.a;
}

// since LCS blocks were generated for each fragment
// a LCS block will always covered by a fragment
#define shift_lcs_block(i) do { \
    begs[(i)] = lcs_blocks[(i)].a[lcsb[(i)]] >> 32; \
    lens[(i)] = (uint32_t) lcs_blocks[(i)].a[lcsb[(i)]]; \
    while (ra[(i)].a[frgs[(i)]].s_end < begs[(i)]) ++frgs[(i)]; \
    uids[(i)] = ra[(i)].a[frgs[(i)]].uid >> 1; \
} while (0)

static void make_ma_block(scg_t *g, sr_t *sr, scg_ra_t *ra, uint32_t n, ma_t *ma)
{
    uint32_t i, j, n_scm, n_frg;
    uint64_t uid, *scm, *u_scm, *lcs_b, s_beg, u_n, s_n;
    double score, intpart;
    int max_d, n_b;
    ra_frg_t *frg;
    kh_u128_t *h;

    // collect syncmers on read
    h = (kh_u128_t *) g->h_scm;
    n_scm = sr->n;
    MYMALLOC(scm, n_scm);
    for (i = 0; i < n_scm; ++i)
        scm[i] = kh_val(h, kh_u128_get(h, sr->k_mer_h[i])) << 1;

    score = modf(ra[0].s, &intpart);
    if (score < FLT_EPSILON) score = 1.0;
    max_d = n_scm - lround(ra[0].s - score);
    
    // find LCS blocks for each alignment record
    lcs_block_t lcs_blocks[n];
    for (i = 0; i < n; ++i) {
        n_frg = ra[i].n;
        kv_init(lcs_blocks[i]);

        for (j = 0; j < n_frg; ++j) {
            frg = &ra[i].a[j];
            uid = frg->uid >> 1;
            u_n = frg->u_end - frg->u_beg + 1;
            s_n = frg->s_end - frg->s_beg + 1;
            if (frg->uid & 1) {
                // mapped to reverse strand
                MYMALLOC(u_scm, u_n);
                memcpy(u_scm, &g->utg_asmg->vtx[uid].a[frg->u_beg], sizeof(uint64_t) * u_n);
                array_reverse(u_scm, u_n);
            } else {
                u_scm = &g->utg_asmg->vtx[uid].a[frg->u_beg];
            }
            // find LCS
            lcs_b = find_lcs(&scm[frg->s_beg], s_n, u_scm, u_n, frg->s_beg, max_d, &n_b);
            // add LCS blocks
            kv_pushn(uint64_t, lcs_blocks[i], lcs_b, n_b);
            free(lcs_b);
            if (frg->uid & 1) free(u_scm);
        }
        // merge blocks
        // lcs_block_merge(&lcs_blocks[i]);
    }

    // make multiple alignment
    uint64_t uids[n], frgs[n], lcsb[n], begs[n], lens[n];
    int m_ext, ext;
    kvec_t(uint32_t) n_match;
    kvec_t(uint64_t) u_match;
    kv_init(n_match);
    kv_init(u_match);
    MYBZERO(lcsb, n);
    MYBZERO(frgs, n);

    for (i = 0; i < n; ++i) {
        if (lcs_blocks[i].n == 0)
            goto ma_done;
        shift_lcs_block(i);
    }

    while (1) {
        // find left alignment position
        s_beg = 0;
        for (i = 0; i < n; ++i)
            s_beg = MAX(s_beg, begs[i]);
        // find max extension
        m_ext = INT_MAX;
        for (i = 0; i < n; ++i) {
            ext = (int) lens[i] - s_beg + begs[i];
            m_ext = MIN(m_ext, ext);
        }
        if (m_ext > 0) {
            // add ma block
            kv_push(uint32_t, n_match, m_ext);
            for (i = 0; i < n; ++i)
                kv_push(uint64_t, u_match, uids[i]);
            // shift LCS blocks
            for (i = 0; i < n; ++i) {
                ext = (int) lens[i] - s_beg + begs[i];
                if (ext == m_ext) {
                    // shift the block
                    ++lcsb[i];
                    if (lcsb[i] == lcs_blocks[i].n)
                        goto ma_done;
                    shift_lcs_block(i);
                } else {
                    begs[i] = s_beg + m_ext;
                    lens[i] = ext - m_ext;
                }
            }
        } else {
            // shift the block with smallest begs
            for (i = 0, j = 1; j < n; ++j)
                if (begs[j] < begs[i])
                    i = j;
            ++lcsb[i];
            if (lcsb[i] == lcs_blocks[i].n)
                goto ma_done;
            shift_lcs_block(i);
        }
    }

ma_done:
    for (i = 0; i < n; ++i)
        kv_destroy(lcs_blocks[i]);
    free(scm);

    ma->a = n;
    ma->b = n_match.n;
    ma->n = n_match.a;
    ma->u = u_match.a;

#ifdef DEBUG_LCS_MALIGN
    fprintf(stderr, "[DEBUG_LCS_MALIGN::%s] MA BLOCKs of read %lu - %s\n", __func__, sr->sid, sr->sname);
    for (i = 0; i < ma->a; ++i) {
        fprintf(stderr, "[DEBUG_LCS_MALIGN::%s] MA_RECORD %u -", __func__, i);
        for (j = 0; j < ma->b; ++j)
            fprintf(stderr, " BLOCK %u: %lu %u;", j, ma->u[j * (ma->a) + i], ma->n[j]);
        fputc('\n', stderr);
    }
#endif
}

#define EM_MAX_ITER 1000

void scg_ra_utg_coverage(scg_t *g, sr_v *sr_v, scg_ra_v *ra_v, int verbose)
{
    if (ra_v->n == 0) {
        fprintf(stderr, "[W::%s] no read alignment, unitig coverage estimation skipped\n", __func__);
        return;
    }

    uint64_t i, j, k, s, n, m, n_vtx, n_scm, uid, sid, *a;
    asmg_t *utg_g;

    utg_g = g->utg_asmg;
    n_vtx = utg_g->n_vtx;

    /***
    // first round estimation with single copy syncmers
    double *avg_covs;
    kvec_t(double) kcovs;
    
    // FIXME utgs with no single copy mers
    MYCALLOC(avg_covs, n_vtx);
    kv_init(kcovs);
    for (i = 0; i < n_vtx; ++i) {
        a = utg_g->vtx[i].a;
        n = utg_g->vtx[i].n;
        for (j = 0; j < n; ++j) {
            s = a[j] >> 1;
            if (scm_utg_n(g, s) == 1)
                kv_push(double, kcovs, (double) g->scm[s].k_cov);
        }
        avg_covs[i] = average_IQR(kcovs.a, kcovs.n, 0);
        avg_covs[i] = MAX(1., avg_covs[i]);
        kcovs.n = 0;
    }
    kv_destroy(kcovs);
    **/
    
    // first round estimation with uniquely mapped reads
    double *avg_covs, *covs, **C, intpart;
    ra_frg_t *ra;
    n_scm = 0;
    for (i = 0; i < n_vtx; ++i)
        n_scm += utg_g->vtx[i].n;
    MYCALLOC(covs, n_scm);
    MYMALLOC(C, n_vtx);
    C[0] = covs;
    for (i = 1; i < n_vtx; ++i)
        C[i] = C[i - 1] + utg_g->vtx[i - 1].n;
    for (i = 0, n = ra_v->n; i < n; ++i) {
        if (modf(ra_v->a[i].s, &intpart) > DBL_EPSILON)
            continue;
        ra = ra_v->a[i].a;
        m = ra_v->a[i].n;
        for (j = 0; j < m; ++j) {
            uid = ra[j].uid >> 1;
            for (k = ra[j].u_beg, s = ra[j].u_end; k <= s; ++k)
                C[uid][k] += 1.;
        }
    }
    MYCALLOC(avg_covs, n_vtx);
    for (i = 0; i < n_vtx; ++i) {
        n = utg_g->vtx[i].n;
        qsort(C[i], n, sizeof(double), dbl_cmpfunc);
        // do not use syncmers not covered
        m = 0;
        while (m < n && C[i][m] < DBL_EPSILON) ++m;
        avg_covs[i] = average_IQR(&C[i][m], n - m, 1);
        avg_covs[i] = MAX(1., avg_covs[i]);
    }
    free(C);
    free(covs);

#ifdef DEBUG_UTG_RA_COV
    fprintf(stderr, "[DEBUG_UTG_RA_COV::%s] unitig coverage estimation first round -\n", __func__);
    for (i = 0; i < n_vtx; ++i)
        fprintf(stderr, "[DEBUG_UTG_RA_COV::%s] u%lu: %u\n", __func__, i, (uint32_t) avg_covs[i]);
#endif

    // find alignment blocks for multiple alignments
    // here assume alignment records were sorted by reads
    kvec_t(ma_t) mas = {0, 0, 0};
    ma_t *ma;
    sid = ra_v->a[0].sid;
    for (i = j = 0; i < ra_v->n; ++i) {
        if (sid != ra_v->a[i].sid) {
            // new read
            kv_pushp(ma_t, mas, &ma);
            make_ma_block(g, &sr_v->a[ra_v->a[j].sid], &ra_v->a[j], i - j, ma);
            j = i;
            sid = ra_v->a[j].sid;
        }
    }
    kv_pushp(ma_t, mas, &ma);
    make_ma_block(g, &sr_v->a[ra_v->a[j].sid], &ra_v->a[j], i - j, ma);

    // second round estimation with mapped reads
    double diff, covt;
    // this is an interface for better EM algorithm
    // currently it converges quickly
    MYMALLOC(covs, n_vtx);
    for (i = 0; i < EM_MAX_ITER; ++i) {
        MYBZERO(covs, n_vtx);
        for (j = 0; j < mas.n; ++j) {
            ma = &mas.a[j];
            for (k = 0; k < ma->b; ++k) {
                a = &ma->u[k * (ma->a)];
                covt = 0.;
                for (n = 0; n < ma->a; ++n)
                    covt += avg_covs[a[n]];
                if (covt == 0.) continue;
                for (n = 0; n < ma->a; ++n)
                    covs[a[n]] += avg_covs[a[n]] / covt * ma->n[k];
            }
        }
        diff = 0.;
        for (j = 0; j < n_vtx; ++j) {
            covt = covs[j] / utg_g->vtx[j].n;
            diff += fabs(covt - avg_covs[j]);
            avg_covs[j] = covt;
        }
        if (verbose > 2)
            fprintf(stderr, "[M::%s] unitig coverage estimation iteration %lu: diff = %.6f\n",
                    __func__, i, diff);
        if (diff < DBL_EPSILON)
            break;
    }
    if (verbose > 2)
        fprintf(stderr, "[M::%s] unitig coverage estimation ended at iteration %lu\n", __func__, i);
    free(covs);

#ifdef DEBUG_UTG_RA_COV
    fprintf(stderr, "[DEBUG_UTG_RA_COV::%s] unitig coverage estimation second round -\n", __func__);
    for (i = 0; i < n_vtx; ++i)
        fprintf(stderr, "[DEBUG_UTG_RA_COV::%s] u%lu: %u\n", __func__, i, (uint32_t) avg_covs[i]);
#endif

    // third round estimation with all syncmers weighted by utg coverage estimated in last round
    uint128_t **idx_u;
    idx_u = g->idx_u;
    MYCALLOC(covs, n_scm);
    MYMALLOC(C, n_vtx);
    C[0] = covs;
    for (i = 1; i < n_vtx; ++i)
        C[i] = C[i - 1] + utg_g->vtx[i - 1].n;
    for (i = 0; i < g->n_scm; ++i) {
        m = scm_utg_n(g, i);
        if (m == 0) continue;
        covt = 0.;
        for (j = 0; j < m; ++j)
            covt += avg_covs[scm_utg_uid(idx_u[i][j])];
        if (covt < DBL_EPSILON) continue;
        for (j = 0; j < m; ++j) {
            uid = scm_utg_uid(idx_u[i][j]);
            C[uid][scm_utg_pos(idx_u[i][j])] = avg_covs[uid] / covt * g->scm[i].k_cov;
        }
    }
    for (i = 0; i < n_vtx; ++i) {
        n = utg_g->vtx[i].n;
        avg_covs[i] = average_IQR(C[i], n, 0);
        avg_covs[i] = MAX(1., avg_covs[i]);
    }
    free(C);
    free(covs);

#ifdef DEBUG_UTG_RA_COV
    fprintf(stderr, "[DEBUG_UTG_RA_COV::%s] unitig coverage estimation third round -\n", __func__);
    for (i = 0; i < n_vtx; ++i)
        fprintf(stderr, "[DEBUG_UTG_RA_COV::%s] u%lu: %u\n", __func__, i, (uint32_t) avg_covs[i]);
#endif

    // update unitig coverage on graph
    for (i = 0; i < n_vtx; ++i)
        utg_g->vtx[i].cov = (uint32_t) avg_covs[i]; 

    for (i = 0; i < mas.n; ++i) {
        free(mas.a[i].n);
        free(mas.a[i].u);
    }
    kv_destroy(mas);

    free(avg_covs);
}

void scg_ra_arc_coverage(scg_t *g, int fix, int verbose)
{
    // update arc coverage on graph
    // the arc coverage on the input graph is the number of syncmer pair links in the read set
    // need to check if there are other pair links on the graph
    // if yes, calculate weighted arc coverage by utg coverge
    
    uint64_t i, j, n, m, p, v, w, c, nl, ld, *a;
    asmg_t *utg_g;
    asmg_arc_t *arc;
    u64_v_t *link_pairs, *pair;
    khint_t k;
    kh_u128_t *h_arc;
    int absent;

    utg_g = g->utg_asmg;
    h_arc = kh_u128_init();
    nl = asmg_max_link_id(utg_g) + 1; // number of links
    MYCALLOC(link_pairs, nl);
    // add arcs
    for (i = 0, n = utg_g->n_arc; i < n; ++i) {
        arc = &utg_g->arc[i];
        if (arc->del || arc->comp) continue;
        v = asmg_arc_head_e(utg_g, arc);
        w = asmg_arc_tail_e(utg_g, arc);
        if (v > w) {
            v ^= 1;
            w ^= 1;
            SWAP(v, w);
        }
        ld = arc->link_id;
        // add arc
        k = kh_u128_put(h_arc, (uint128_t) v << 64 | w, &absent);
        // only update if absent
        // for (v,w) with multiple link ids
        // all pack in the first link
        if (absent)
            kh_val(h_arc, k) = ld;
        else
            ld = kh_val(h_arc, k);
        kv_push(uint64_t, link_pairs[ld], arc->link_id);
        kv_push(uint64_t, link_pairs[ld], ((uint64_t) 
                    utg_g->vtx[arc->v >> 1].cov + 
                    utg_g->vtx[arc->w >> 1].cov) / 2);
    }
    // check utgs
    for (i = 0, n = utg_g->n_vtx; i < n; ++i) {
        a = utg_g->vtx[i].a;
        for (j = 1, m = utg_g->vtx[i].n; j < m; ++j) {
            v = a[j - 1];
            w = a[j];
            if (v > w) {
                v ^= 1;
                w ^= 1;
                SWAP(v, w);
            }
            k = kh_u128_get(h_arc, (uint128_t) v << 64 | w);
            if (k == kh_end(h_arc)) continue;
            ld = kh_val(h_arc, k);
            kv_push(uint64_t, link_pairs[ld], UINT64_MAX);
            kv_push(uint64_t, link_pairs[ld], utg_g->vtx[i].cov);
        } 
    }
    // update arc covs on graph
    for (i = 0, n = utg_g->n_arc; i < n; ++i) {
        arc = &utg_g->arc[i];
        if (arc->del || arc->comp) continue;
        v = asmg_arc_head_e(utg_g, arc);
        w = asmg_arc_tail_e(utg_g, arc);
        if (v > w) {
            v ^= 1;
            w ^= 1;
            SWAP(v, w);
        }
        k = kh_u128_get(h_arc, (uint128_t) v << 64 | w);
        ld = kh_val(h_arc, k);
        pair = &link_pairs[ld];
        if (pair->n == 2) {
            assert(pair->a[0] == ld);
            continue; // no need to update
        }
        // find the pair and total cov
        c = 0;
        p = UINT64_MAX;
        for (j = 1, m = pair->n; j < m; j += 2) {
            if (pair->a[j - 1] == arc->link_id)
                p = j;
            c += pair->a[j];
        }
        if (c == 0) continue;
        assert(p < m);
        c = lround((double) arc->cov / c * pair->a[p]);
        // c = MIN(c, utg_g->vtx[arc->v >> 1].cov);
        // c = MIN(c, utg_g->vtx[arc->w >> 1].cov);
        if (verbose > 2)
            fprintf(stderr, "[M::%s] arc u%lu%c -> u%lu%c coverage updated: %u -> %lu\n", 
                    __func__, arc->v>>1, "+-"[arc->v&1], arc->w>>1, "+-"[arc->w&1], arc->cov, c);
        arc->cov = c;

        // for reverse complementary arcs
        arc = asmg_comp_arc1(utg_g, arc);
        if (arc) arc->cov = c;
    }

    // fix arc coverage
    // bounded by the minimum utg cov
    if (fix) {
        for (i = 0, n = utg_g->n_arc; i < n; ++i) {
            arc = &utg_g->arc[i];
            if (arc->del) continue;
            c = MIN(utg_g->vtx[arc->v >> 1].cov, utg_g->vtx[arc->w >> 1].cov);
            if (c >= arc->cov) continue;
            if (verbose > 2)
                fprintf(stderr, "[M::%s] arc u%lu%c -> u%lu%c coverage fixed: %u -> %lu\n",
                        __func__, arc->v>>1, "+-"[arc->v&1], arc->w>>1, "+-"[arc->w&1], arc->cov, c);
            arc->cov = c;
        }
    }

    for (i = 0; i < nl; ++i)
        free(link_pairs[i].a);
    free(link_pairs);
    kh_u128_destroy(h_arc);
}

