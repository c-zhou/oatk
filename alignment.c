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
 * 19/01/23 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <pthread.h>

#include "khashl.h"
#include "kvec.h"

#include "syncmer.h"
#include "syncasm.h"
#include "misc.h"

#undef DEBUG_READ_ALIGNMENT

static kh_inline khint_t kh_hash_uint128(uint128_t key)
{
    khint_t k1 = kh_hash_uint64((khint64_t) key);
    khint_t k2 = kh_hash_uint64((khint64_t) (key >> 64));
    return kh_hash_uint64((khint64_t) ((uint64_t) k1 << 32 | k2));
}

KHASHL_MAP_INIT(KH_LOCAL, kh_128_t, kh_128, uint128_t, uint64_t, kh_hash_uint128, kh_eq_generic)

void scg_ra_v_destroy(scg_ra_v *ra_v)
{
    if (!ra_v) return;
    size_t i;
    for (i = 0; i < ra_v->n; ++i)
        free(ra_v->a[i].a);
    kv_destroy(*ra_v);
}

typedef struct { size_t n, m; uint32_t *a; } u32_v;
typedef struct { size_t n, m; u32_v *a; } v32_v;

typedef struct {
    uint64_t uid;
    uint64_t u_pos;
    uint64_t s_pos;
} sr_scm_t;

typedef struct {
    uint64_t uid;
    uint64_t u_beg, u_end, u_gap;
    uint64_t s_beg, s_end, s_cnt;
    int64_t score0, score;
    u32_v prev;
} sr_frg_t;

typedef struct {
    int tid;
    sr_t *sr;
    scg_ra_v *ra_v;
    scg_t *g;
    int64_t *old_ra;
    uint64_t n;
    uint64_t m_stats[4];
} p_data_t;

pthread_mutex_t lock;

static int sr_scm_cmpfunc(const void *a, const void *b)
{
    uint64_t x, y;
    x = ((sr_scm_t *) a)->uid;
    y = ((sr_scm_t *) b)->uid;
    if (x != y)
        return (x > y) - (x < y);
    x = ((sr_scm_t *) a)->s_pos;
    y = ((sr_scm_t *) b)->s_pos;
    if (x != y)
        return (x > y) - (x < y);
    x = ((sr_scm_t *) a)->u_pos;
    y = ((sr_scm_t *) b)->u_pos;
    return (x > y) - (x < y);
}

static int sr_frg_cmpfunc(const void *a, const void *b)
{
    uint64_t x, y;
    x = ((sr_frg_t *) a)->s_beg;
    y = ((sr_frg_t *) b)->s_beg;
    if (x != y)
        return (x > y) - (x < y);
    x = ((sr_frg_t *) a)->s_end;
    y = ((sr_frg_t *) b)->s_end;
    return (x > y) - (x < y);
}

static inline void rev_array(uint32_t *arr, int n)
{
    int i, j;
    uint32_t tmp;
    for (i = 0, j = n - 1; i < j; ++i, --j) {
        tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
}

static void aln_frg_backtrace(uint32_t node, uint32_t len, sr_frg_t *frg_v, u32_v *utg_v, v32_v *aln_v)
{
    if (node == UINT32_MAX)
        return;

    kv_push(uint32_t, *utg_v, node);
    ++len;

    if (frg_v[node].prev.n == 0) {
        // no prev node
        // add aln path to aln_v
        u32_v *v;
        kv_pushp(u32_v, *aln_v, &v);
        v->n = v->m = len;
        MYMALLOC(v->a, len);
        memcpy(v->a, utg_v->a, sizeof(uint32_t) * len);
        rev_array(v->a, v->n);
    } else {
        // add paths for prev nodes
        size_t i, n;
        for (i = 0, n = frg_v[node].prev.n; i < n; ++i) {
            aln_frg_backtrace(frg_v[node].prev.a[i], len, frg_v, utg_v, aln_v);
            utg_v->n = len;
        }
    }
}

static int match_score = 1;
static int gap_penalty = 1;
// static int clip_penalty = 1;
static double min_a_frac = .9;

static void *scg_ra_analysis_thread(void *args)
{
    p_data_t *dat = (p_data_t *) args;

    if (dat->n == 0)
        return NULL;

    sr_t *sr;
    scg_ra_v *ra_v;
    scg_ra_t *ra;
    scg_t *g;
    scg_utg_t *utg;
    asmg_arc_t *arc;
    uint64_t i, k, n, s, u, p, p1;
    int64_t score, score1, max_score, u_gap, u_clip, u_ovl, *old_ra;
    uint32_t j, t, m, n_m, n_u, n_a;
    int incr, d;
    uint128_t x;
    kh_128_t *h;
    sr_scm_t *scm;
    sr_frg_t *frg, *frg1;
    kvec_t(sr_scm_t) scm_v; // scm position triple (utg_id, utg_pos, sr_pos)
    kvec_t(sr_frg_t) frg_v;
    u32_v utg_v, *v;
    v32_v aln_v;

    g = dat->g;
    ra_v = dat->ra_v;
    old_ra = dat->old_ra;
    h = (kh_128_t *) g->h_scm;
    utg = g->utg_asmg->vtx;

    kv_init(scm_v);
    kv_init(frg_v);
    kv_init(utg_v);
    kv_init(aln_v);
    n_m = n_u = 0;
    for (i = 0, n = dat->n; i < n; ++i) {
        if ((old_ra[i]&1) == 0)
            continue;

        sr = &dat->sr[i];

        if (sr->n == 0)
            continue;

        // find all syncmer positions
        scm_v.n = 0; // reset scm vector
        for (j = 0, m = sr->n; j < m; ++j) {
            s = kh_val(h, kh_128_get(h, sr->k_mer_h[j]));
            u = scm_utg_n(g, s);

            for (k = 0; k < u; ++k) {
                x = g->idx_u[s][k];
                kv_pushp(sr_scm_t, scm_v, &scm);
                scm->uid = scm_utg_uid(x)<<1 | (scm_utg_rev(x)^(sr->m_pos[j]&1));
                scm->u_pos = scm_utg_pos(x);
                scm->s_pos = j;
            }
        }

        // sort scm triple by uid - s_pos - u_pos
        qsort(scm_v.a, scm_v.n, sizeof(sr_scm_t), sr_scm_cmpfunc);

        // collect all alignment fragments
        frg_v.n = 0; // reset frg vector
        for (j = 0, m = scm_v.n; j < m; ) {
            s = j; // fragment start position
            u = scm_v.a[j].uid;
            p = scm_v.a[j].u_pos;

            // check alignment along a utg
            incr = 0;
            while (++j < m) {
                if (scm_v.a[j].uid != u) // reach a new utg
                    break;
                // mapping positions need to change monotonically
                p1 = scm_v.a[j].u_pos;
                d = (p<p1) - (p>p1);
                if (d == 0) // this is possible for singleton utgs
                    break;
                if (d * incr < 0) {
                    // mapping direction changed
                    // there is a special case when with only two syncmers
                    // need to decide if the second syncmer should be with
                    // the first one or the next one
                    // decision is made by comparing the gap size
                    if (j - s == 2) {
                        if (((int64_t) scm_v.a[s+1].u_pos - scm_v.a[s].u_pos) * incr >
                                ((int64_t) scm_v.a[j].u_pos - scm_v.a[j-1].u_pos) * d)
                            --j;
                    }
                    break;
                }
                p = p1;
                incr = d;
            }

            // calculate gaps and score
            u_gap = 0;
            for (k = s + 1; k < j; ++k)
                u_gap += labs(scm_v.a[k].u_pos - scm_v.a[k-1].u_pos) - 1;
            if (u_gap < 0) u_gap = 0;
            score = ((int64_t) j - s) * match_score - u_gap * gap_penalty;

            if (score >= 0 ) {
                // push fragment
                kv_pushp(sr_frg_t, frg_v, &frg);
                frg->uid = u;
                if (incr > 0) {
                    frg->u_beg = scm_v.a[s].u_pos;
                    frg->u_end = scm_v.a[j-1].u_pos;
                } else {
                    frg->u_beg = scm_v.a[j-1].u_pos;
                    frg->u_end = scm_v.a[s].u_pos;
                }
                frg->s_beg = scm_v.a[s].s_pos;
                frg->s_end = scm_v.a[j-1].s_pos;
                frg->u_gap = u_gap;
                frg->s_cnt = j - s;
                frg->score0 = score;
                frg->score = score;
                kv_init(frg->prev);
            }
        }

#ifdef DEBUG_READ_ALIGNMENT
        uint64_t dbg_nr = 303579; //866; //910584;
        pthread_mutex_lock(&lock);
        if (sr->sid == dbg_nr) {
            fprintf(stderr, "[DEBUG_READ_ALIGNMENT::%s_%d] NO_READ: %lu\n", __func__, dat->tid, sr->sid);
            for (j = 0, m = sr->n; j < m; ++j)
                fprintf(stderr, "[DEBUG_READ_ALIGNMENT::%s_%d] M_POS: %u %u %u\n", __func__, dat->tid, j, sr->m_pos[j]>>1, sr->m_pos[j]&1);
            for (j = 0, m = scm_v.n; j < m; ++j)
                fprintf(stderr, "[DEBUG_READ_ALIGNMENT::%s_%d] U_POS: %lu%c %lu %lu\n", __func__, dat->tid,
                        scm_v.a[j].uid>>1, "+-"[scm_v.a[j].uid&1], scm_v.a[j].u_pos, scm_v.a[j].s_pos);
        }
        pthread_mutex_unlock(&lock);
#endif

        if (frg_v.n == 0) continue;

        // sort fragments by s_beg - s_end
        qsort(frg_v.a, frg_v.n, sizeof(sr_frg_t), sr_frg_cmpfunc);

        // chaining
        m = frg_v.n;
        max_score = 0;
        for (j = 0; j < m; ++j) {
            frg = &frg_v.a[j];
            p = frg->s_end;
            // score = frg->score - ((frg->uid&1)? frg->u_beg : (utg[frg->uid>>1].n - frg->u_end - 1)) * clip_penalty; // utg clip
            u_clip = (frg->uid&1)? frg->u_beg : (utg[frg->uid>>1].n - frg->u_end - 1);
            if (u_clip > 0) // no clip allowed
                continue;
            score = frg->score;
            for (k = j + 1; k < m; ++k) {
                frg1 = &frg_v.a[k];
                u_clip = (frg1->uid&1)? (utg[frg1->uid>>1].n - frg1->u_end - 1) : frg1->u_beg;
                if (u_clip > 0)
                    continue; // no clip allowed
                arc = asmg_arc1(g->utg_asmg, frg->uid, frg1->uid);
                if (arc == 0)
                    continue;
                u_ovl = MIN(arc->lo, p + 1);
                p1 = frg1->s_beg;
                if (p1 > p + 1)
                    break;
                if (p1 + u_ovl != p + 1) // no overlap or gap allowed
                    continue;
                // score1 = score + frg1->score0 - ((frg1->uid&1)? (utg[frg1->uid>>1].n - frg1->u_end - 1) : frg1->u_beg) * clip_penalty;
                score1 = score + frg1->score0 - u_ovl * match_score;
                if (score1 < frg1->score || (score1 == frg1->score && frg1->prev.n == 0))
                    continue;
                if (score1 > frg1->score) {
                    // better chain
                    // clean prev
                    frg1->score = score1;
                    frg1->prev.n = 0;
                }
                kv_push(uint32_t, frg1->prev, j);
#ifdef DEBUG_READ_ALIGNMENT
                pthread_mutex_lock(&lock);
                if (sr->sid == dbg_nr) {
                    uint32_t j1;
                    for (j1 = 0; j1 < m; ++j1)
                        fprintf(stderr, "[DEBUG_READ_ALIGNMENT::%s_%d] FRG_V [%u %lu]: %lu%c [%lu] %lu %lu [%lu] %lu %lu [%lu] (%ld) <%lu:%u>\n",
                                __func__, dat->tid, j, k, frg_v.a[j1].uid>>1, "+-"[frg_v.a[j1].uid&1], utg[frg_v.a[j1].uid>>1].n,
                                frg_v.a[j1].u_beg, frg_v.a[j1].u_end, frg_v.a[j1].u_gap,
                                frg_v.a[j1].s_beg, frg_v.a[j1].s_end, frg_v.a[j1].s_cnt,
                                frg_v.a[j1].score, frg_v.a[j1].prev.n, frg_v.a[j1].prev.n > 0? frg_v.a[j1].prev.a[0] : 0);
                }
                pthread_mutex_unlock(&lock);
#endif
            }
            if (max_score < frg->score)
                max_score = frg->score;
        }
#ifdef DEBUG_READ_ALIGNMENT
        pthread_mutex_lock(&lock);
        if (sr->sid == dbg_nr) {
            for (j = 0; j < m; ++j)
                fprintf(stderr, "[DEBUG_READ_ALIGNMENT::%s_%d] FRG_V: %lu%c [%lu] %lu %lu [%lu] %lu %lu [%lu] (%ld) <%lu:%u>\n",
                        __func__, dat->tid, frg_v.a[j].uid>>1, "+-"[frg_v.a[j].uid&1], utg[frg_v.a[j].uid>>1].n,
                        frg_v.a[j].u_beg, frg_v.a[j].u_end, frg_v.a[j].u_gap,
                        frg_v.a[j].s_beg, frg_v.a[j].s_end, frg_v.a[j].s_cnt,
                        frg_v.a[j].score, frg_v.a[j].prev.n, frg_v.a[j].prev.n > 0? frg_v.a[j].prev.a[0] : 0);
        }
        pthread_mutex_unlock(&lock);
#endif

        // backtrace to get the alignment
        aln_v.n = 0;

        // only do alignment if score is no samller than the old alignment
        if (max_score >= (old_ra[i] >> 1)) {
            for (j = 0; j < m; ++j) {
                frg = &frg_v.a[j];
                if (frg->score < max_score)
                    continue;
                utg_v.n = 0;
                aln_frg_backtrace(j, 0, frg_v.a, &utg_v, &aln_v);
            }
        }

#ifdef DEBUG_READ_ALIGNMENT
        pthread_mutex_lock(&lock);
        if (sr->sid == dbg_nr) {
            fprintf(stderr, "[DEBUG_READ_ALIGNMENT::%s_%d] NO_ALN: %lu\n", __func__, dat->tid, aln_v.n);
            for (j = 0; j < aln_v.n; ++j) {
                v = &aln_v.a[j];
                for (k = 0; k < v->n; ++k) {
                    t = v->a[k];
                    fprintf(stderr, "[DEBUG_READ_ALIGNMENT::%s_%d] FRG_V %lu/%lu: %lu%c [%lu] %lu %lu [%lu] %lu %lu [%lu] (%ld)\n",
                            __func__, dat->tid, k, v->n,
                            frg_v.a[t].uid>>1, "+-"[frg_v.a[t].uid&1], utg[frg_v.a[t].uid>>1].n,
                            frg_v.a[t].u_beg, frg_v.a[t].u_end, frg_v.a[t].u_gap,
                            frg_v.a[t].s_beg, frg_v.a[t].s_end, frg_v.a[t].s_cnt,
                            frg_v.a[t].score);
                }
            }
        }
        pthread_mutex_unlock(&lock);
#endif

        // collect alignments
        n_a = 0;
        for (j = 0; j < aln_v.n; ++j) {
            v = &aln_v.a[j];
            s = 0;
            for (k = 0; k < v->n; ++k) {
                t = v->a[k];
                s += frg_v.a[t].s_cnt;
            }

            if ((double) s / sr->n < min_a_frac)
                continue;

            // add new alignment to ra_v
            kv_pushp(scg_ra_t, *ra_v, &ra);
            ra->sid = sr->sid;
            ra->n = v->n;
            MYMALLOC(ra->a, ra->n);
            for (k = 0; k < ra->n; ++k) {
                t = v->a[k];
                ra->a[k].uid = frg_v.a[t].uid;
                ra->a[k].u_beg = frg_v.a[t].u_beg;
                ra->a[k].u_end = frg_v.a[t].u_end;
                ra->a[k].s_beg = frg_v.a[t].s_beg;
                ra->a[k].s_end = frg_v.a[t].s_end;
            }
            ++n_a;
        }

        for (j = 0; j < n_a; ++j)
            ra_v->a[ra_v->n-j-1].s = 1.0 / n_a + max_score;

        if (n_a > 0)  ++n_m;
        if (n_a == 1) ++n_u;

        for (j = 0; j < aln_v.n; ++j)
            kv_destroy(aln_v.a[j]);

        for (j = 0; j < m; ++j)
            kv_destroy(frg_v.a[j].prev);
    }

    dat->m_stats[0] = n_m;
    dat->m_stats[1] = n_u;

    kv_destroy(scm_v);
    kv_destroy(frg_v);
    kv_destroy(utg_v);
    kv_destroy(aln_v);

    return NULL;
}

void scg_read_alignment(sr_v *sr, scg_ra_v *ra_v, scg_t *g, int n_threads, int for_unzip)
{
    int i;
    uint64_t j, b;

    if (n_threads == 0)
        n_threads = 1;

    if (pthread_mutex_init(&lock, NULL) != 0) {
        fprintf(stderr, "[E::%s] mutex init failed\n", __func__);
        return;
    }

    int64_t *old_ra;
    scg_ra_t *ra;
    MYCALLOC(old_ra, sr->n);
    if (for_unzip && ra_v->n > 0) {
        // for unzipping purpose
        // only map reads previously mapped to at least three unitigs
        // only map previouly mapped reads
        // only keep mapping results without mapping score decrease
        double fractpart, intpart;
        uint64_t sid;
        for (j = 0; j < ra_v->n; ++j) {
            ra = &ra_v->a[j];
            sid = ra->sid;
            if (ra->n > 2 && (old_ra[sid] & 1) == 0) {
                fractpart = modf(ra->s, &intpart);
                if (fractpart < DBL_EPSILON)
                    intpart -= 1;
                old_ra[sid] = (uint64_t) intpart << 1 | 1;
            }
        }
    } else {
        // map every reads
        for (j = 0; j < sr->n; ++j) old_ra[j] = 1;
    }

    b = sr->n / n_threads + 1;
    p_data_t *dat;
    MYCALLOC(dat, n_threads);
    for (i = 0; i < n_threads; ++i) {
        dat[i].tid = i;
        dat[i].sr = sr->a + i * b;
        dat[i].old_ra = old_ra + i * b;
        MYCALLOC(dat[i].ra_v, 1);
        dat[i].g = g;
        if (sr->n > (i + 1) * b) {
            dat[i].n = b;
        } else {
            dat[i].n = sr->n - i * b;
            break;
        }
    }

    pthread_t threads[n_threads];
    for (i = 1; i < n_threads; ++i)
        pthread_create(threads + i, NULL, scg_ra_analysis_thread, dat + i);

    scg_ra_analysis_thread(dat);

    for (i = 1; i < n_threads; ++i)
        pthread_join(threads[i], NULL);

    // clean old results in ra_v
    scg_ra_v_destroy(ra_v);
    // collect results
    ra_v->m = 0;
    for (i = 0; i < n_threads; ++i)
        ra_v->m += dat[i].ra_v->n;
    ra_v->n = 0;
    MYMALLOC(ra_v->a, ra_v->m);
    for (i = 0; i < n_threads; ++i) {
        // copy data
        memcpy(ra_v->a + ra_v->n, dat[i].ra_v->a, sizeof(scg_ra_t) * dat[i].ra_v->n);
        ra_v->n += dat[i].ra_v->n;
        free(dat[i].ra_v);
    }

    uint64_t n_r, n_m, n_u, n;
    n_r = n_m = n_u = 0;
    for (n = 0; n < sr->n; ++n)
        n_r += sr->a[n].n > 0;
    for (i = 0; i < n_threads; ++i) {
        n_m += dat[i].m_stats[0];
        n_u += dat[i].m_stats[1];
    }
    fprintf(stderr ,"[M::%s] %lu mappable reads, %lu mapped (%lu unique mapping)\n", __func__, n_r, n_m, n_u);

    pthread_mutex_destroy(&lock);
    free(old_ra);
    free(dat);
}

void scg_ra_print(scg_ra_t *ra, FILE *fo)
{
    uint32_t i;
    ra_frg_t *a;
    fprintf(fo, "RID %lu [N %u S %.3f]:", ra->sid, ra->n, ra->s);
    for (i = 0; i < ra->n; ++i) {
        a = &ra->a[i];
        fprintf(fo, " [u%lu%c %lu %lu %u %u]", a->uid>>1, "+-"[a->uid&1], a->u_beg, a->u_end, a->s_beg, a->s_end);
    }
    fprintf(fo, "\n");
}

void scg_rv_print(scg_ra_v *rv, FILE *fo)
{
    size_t i;
    for (i = 0; i < rv->n; ++i)
        scg_ra_print(&rv->a[i], fo);
}

