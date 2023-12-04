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

static void scg_ra_v_clean(scg_ra_v *ra_v)
{
    if (!ra_v) return;
    size_t i;
    for (i = 0; i < ra_v->n; ++i)
        free(ra_v->a[i].a);
    kv_destroy(*ra_v);
}

void scg_ra_v_destroy(scg_ra_v *ra_v)
{
    if (!ra_v) return;
    scg_ra_v_clean(ra_v);
    free(ra_v);
}

typedef struct { size_t n, m; uint32_t *a; } u32_v;
typedef struct { size_t n, m; u32_v *a; } v32_v;

typedef struct {
    uint64_t uid;
    uint64_t u_pos;
    uint64_t s_pos;
    uint64_t next;
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

static void set_fragment(sr_frg_t *frg, uint64_t u, 
        uint64_t s_beg, uint64_t s_end, uint64_t s_cnt,
        uint64_t u_beg, uint64_t u_end, uint64_t u_gap,
        int64_t score) {
    frg->uid = u;
    frg->s_beg = s_beg;
    frg->s_end = s_end;
    frg->s_cnt = s_cnt;
    frg->u_beg = u_beg;
    frg->u_end = u_end;
    frg->u_gap = u_gap;
    frg->score0 = score;
    frg->score = score;
    kv_init(frg->prev);
}

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
    uint64_t i, j, k, m, n, s, t, u, p, p1, u_beg, u_end, s_beg, s_end, s_cnt;
    int64_t score, score1, max_score, u_gap, u_clip, u_ovl, s_gap, *old_ra;
    uint32_t n_m, n_u, n_a;
    uint128_t x;
    sr_scm_t *scm;
    sr_frg_t *frg, *frg1;
    kvec_t(sr_scm_t) scm_v; // scm position triple (utg_id, utg_pos, sr_pos)
    kvec_t(sr_frg_t) frg_v;
    u32_v utg_v, pos_v, *v;
    v32_v aln_v;

    g = dat->g;
    ra_v = dat->ra_v;
    old_ra = dat->old_ra;
    utg = g->utg_asmg->vtx;

    kv_init(scm_v);
    kv_init(frg_v);
    kv_init(utg_v);
    kv_init(pos_v);
    kv_init(aln_v);
    n_m = n_u = 0;

#ifdef DEBUG_READ_ALIGNMENT
    uint64_t dbg_nr = 35; //16694; //31974; //7490; //23; //303579; //866; //910584;
#endif

    for (i = 0, n = dat->n; i < n; ++i) {
        if ((old_ra[i]&1) == 0)
            continue;

        sr = &dat->sr[i];

        if (sr->n == 0)
            continue;

        // find all syncmer positions
        scm_v.n = 0; // reset scm vector
        for (j = 0, m = sr->n; j < m; ++j) {
            s = sr->k_mer[j] >> 1;
            s_cnt = scm_utg_n(g, s);

            for (k = 0; k < s_cnt; ++k) {
                x = g->idx_u[s][k];
                u = scm_utg_uid(x);
                p = scm_utg_pos(x);
                t = scm_utg_rev(x)^(sr->m_pos[j]&1);
                // push scm
                kv_pushp(sr_scm_t, scm_v, &scm);
                scm->uid = u<<1 | t;
                scm->u_pos = t? utg[u].n - p - 1 : p;
                scm->s_pos = j;
                scm->next = 0xFFFFFFFFFFFFFFFEULL;
            }
        }

        if (scm_v.n == 0)
            continue;

        // sort scm triple by uid - s_pos - u_pos
        qsort(scm_v.a, scm_v.n, sizeof(sr_scm_t), sr_scm_cmpfunc);
        
        // collect all alignment fragments
        frg_v.n = 0; // reset frg vector
        for (j = 0, m = scm_v.n; j < m; ) {
            u = scm_v.a[j].uid;

            // find start position of next utg
            p = j;
            while (++p < m && scm_v.a[p].uid == u) {}

            // build position index
            pos_v.n = 0;
            kv_push(uint32_t, pos_v, j);
            p1 = scm_v.a[j].s_pos;
            for (k = j + 1; k < p; ++k) {
                if (scm_v.a[k].s_pos != p1) {
                    kv_push(uint32_t, pos_v, k);
                    p1 = scm_v.a[k].s_pos;
                }
            }
            kv_push(uint32_t, pos_v, p);

            // for each scm find next mapping positions
            for (k = 0; k < pos_v.n - 2; ++k) {
                s = pos_v.a[k + 1];
                t = pos_v.a[k + 2];
                uint64_t s1 = pos_v.a[k], t1 = s;
                while (s1 < s) {
                    // find closest larger
                    while (t1 < t && scm_v.a[t1].u_pos <= scm_v.a[s1].u_pos)
                        ++t1;
                    if (t1 < t && scm_v.a[t1].u_pos > scm_v.a[s1].u_pos)
                        scm_v.a[s1].next = t1 << 1;
                    
                    ++s1;
                }
            }
            
            // backtrack to get alignment fragments
            for (k = j; k < p; ++k) {
                s = k;
                if (scm_v.a[s].next & 1) // not a starting point
                    continue;
                u_beg = scm_v.a[s].u_pos;
                s_beg = scm_v.a[s].s_pos;
                s_cnt = 1;
                // calculate gaps and score
                u_gap = s_gap = 0;
                while (1) {
                    t = scm_v.a[s].next >> 1;
                    if (t == 0x7FFFFFFFFFFFFFFFULL)
                        break;
                    u_gap += labs(scm_v.a[t].u_pos - scm_v.a[s].u_pos) - 1;
                    s_gap += labs(scm_v.a[t].s_pos - scm_v.a[s].s_pos) - 1;
                    scm_v.a[s].next |= 1;
                    ++s_cnt;
                    s = t;
                }
                if (s_cnt == 1) // singletons will be processed later
                    continue;
                scm_v.a[s].next |= 1;
                u_end = scm_v.a[s].u_pos;
                s_end = scm_v.a[s].s_pos;
                u_gap = MAX(u_gap, s_gap);
                if (u_gap < 0) u_gap = 0;
                score = (int64_t) s_cnt * match_score - u_gap * gap_penalty;

                if (score >= 0) {
                    // push fragment
                    kv_pushp(sr_frg_t, frg_v, &frg);
                    set_fragment(frg, u, s_beg, s_end, s_cnt, u_beg, u_end, u_gap, score);
                }
            }

            // add singleton scm
            for (k = j; k < p; ++k) {
                if (scm_v.a[k].next == 0xFFFFFFFFFFFFFFFEULL) {
                    scm = &scm_v.a[k];
                    kv_pushp(sr_frg_t, frg_v, &frg);
                    set_fragment(frg, u, scm->s_pos, scm->s_pos, 1, scm->u_pos, scm->u_pos, 0, 1);
                    continue;
                }
            }

            // moving to next utg
            j = p;
        }

        /***
        // collect all alignment fragments
        frg_v.n = 0; // reset frg vector
        for (j = 0, m = scm_v.n; j < m; ) {
            s = j; // fragment start position
            u = scm_v.a[j].uid;
            p = scm_v.a[j].u_pos;

            // check alignment along a utg
            int d, shift = 0, incr = 0;
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
                    shift = 1;
                    break;
                }
                p = p1;
                incr = d;
            }

            // calculate gaps and score
            u_gap = s_gap = 0;
            for (k = s + 1; k < j; ++k) {
                u_gap += labs(scm_v.a[k].u_pos - scm_v.a[k-1].u_pos) - 1;
                s_gap += labs(scm_v.a[k].s_pos - scm_v.a[k-1].s_pos) - 1;
            }
            u_gap = MAX(u_gap, s_gap);
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

            if (j > 0) j -= shift;
        }
        **/

#ifdef DEBUG_READ_ALIGNMENT
        pthread_mutex_lock(&lock);
        if (sr->sid == dbg_nr) {
            fprintf(stderr, "[DEBUG_READ_ALIGNMENT::%s_%d] NO_READ: %lu\n", __func__, dat->tid, sr->sid);
            for (j = 0, m = sr->n; j < m; ++j)
                fprintf(stderr, "[DEBUG_READ_ALIGNMENT::%s_%d] M_POS: %lu %u s%lu%c\n", __func__, dat->tid, 
                        j, sr->m_pos[j]>>1, sr->k_mer[j]>>1, "+-"[sr->m_pos[j]&1]);
            for (j = 0, m = scm_v.n; j < m; ++j)
                fprintf(stderr, "[DEBUG_READ_ALIGNMENT::%s_%d] U_POS [%lu]: %lu%c %lu %lu [%lu]\n", 
                        __func__, dat->tid, j,
                        scm_v.a[j].uid>>1, "+-"[scm_v.a[j].uid&1], scm_v.a[j].u_pos, scm_v.a[j].s_pos,
                        scm_v.a[j].next >> 1);
        }
        
        if (sr->sid == dbg_nr) {
            for (j = 0, m = frg_v.n; j < m; ++j)
                fprintf(stderr, "[DEBUG_READ_ALIGNMENT::%s_%d] FRG_V: %lu%c [%lu] %lu %lu [%lu] %lu %lu [%lu] (%ld)\n",
                        __func__, dat->tid, frg_v.a[j].uid>>1, "+-"[frg_v.a[j].uid&1], utg[frg_v.a[j].uid>>1].n,
                        frg_v.a[j].u_beg, frg_v.a[j].u_end, frg_v.a[j].u_gap,
                        frg_v.a[j].s_beg, frg_v.a[j].s_end, frg_v.a[j].s_cnt,
                        frg_v.a[j].score);
        }
        pthread_mutex_unlock(&lock);
#endif

        if (frg_v.n == 0) continue;

        // sort fragments by s_beg - s_end
        qsort(frg_v.a, frg_v.n, sizeof(sr_frg_t), sr_frg_cmpfunc);

        // chaining
        m = frg_v.n;
        for (j = 0; j < m; ++j) {
            frg = &frg_v.a[j];
            p = frg->s_end;
            // score = frg->score - ((frg->uid&1)? frg->u_beg : (utg[frg->uid>>1].n - frg->u_end - 1)) * clip_penalty; // utg clip
            // u_clip = (frg->uid&1)? frg->u_beg : (utg[frg->uid>>1].n - frg->u_end - 1);
            u_clip = utg[frg->uid>>1].n - frg->u_end - 1;
            if (u_clip > 0) // no clip allowed
                continue;
            score = frg->score;
            for (k = j + 1; k < m; ++k) {
                frg1 = &frg_v.a[k];
                // u_clip = (frg1->uid&1)? (utg[frg1->uid>>1].n - frg1->u_end - 1) : frg1->u_beg;
                u_clip = frg1->u_beg;
                if (u_clip > 0)
                    continue; // no clip allowed
                arc = asmg_arc1(g->utg_asmg, frg->uid, frg1->uid);
                if (arc == 0)
                    continue;
                u_ovl = MIN(arc->ln, p + 1);
                p1 = frg1->s_beg;
                if (p1 > p + 1)
                    break;
                if (p1 + u_ovl != p + 1) // no overlap or gap allowed
                    continue;
                // score1 = score + frg1->score0 - ((frg1->uid&1)? (utg[frg1->uid>>1].n - frg1->u_end - 1) : frg1->u_beg) * clip_penalty;
                score1 = score + frg1->score0 - u_ovl * match_score;
                if (score1 <= score ||
                        score1 < frg1->score || 
                        (score1 == frg1->score && frg1->prev.n == 0))
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
                        fprintf(stderr, "[DEBUG_READ_ALIGNMENT::%s_%d] FRG_V [%lu %lu]: %lu%c [%lu] %lu %lu [%lu] %lu %lu [%lu] (%ld) <%lu:%u>\n",
                                __func__, dat->tid, j, k, frg_v.a[j1].uid>>1, "+-"[frg_v.a[j1].uid&1], utg[frg_v.a[j1].uid>>1].n,
                                frg_v.a[j1].u_beg, frg_v.a[j1].u_end, frg_v.a[j1].u_gap,
                                frg_v.a[j1].s_beg, frg_v.a[j1].s_end, frg_v.a[j1].s_cnt,
                                frg_v.a[j1].score, frg_v.a[j1].prev.n, frg_v.a[j1].prev.n > 0? frg_v.a[j1].prev.a[0] : 0);
                }
                pthread_mutex_unlock(&lock);
#endif
            }
        }
        
        max_score = 0;
        for (j = 0; j < m; ++j) {
            frg = &frg_v.a[j];
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
    kv_destroy(pos_v);
    kv_destroy(aln_v);

    return NULL;
}

void scg_read_alignment(sr_db_t *sr_db, scg_ra_v *ra_v, scg_t *g, int n_threads, int for_unzip)
{
    if (sr_db->n == 0 || !asmg_vtx_n1(g->utg_asmg)) return;

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
    MYCALLOC(old_ra, sr_db->n);
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
        for (j = 0; j < sr_db->n; ++j) old_ra[j] = 1;
    }

    b = sr_db->n / n_threads + 1;
    p_data_t *dat;
    MYCALLOC(dat, n_threads);
    for (i = 0; i < n_threads; ++i) {
        dat[i].tid = i;
        dat[i].sr = sr_db->a + i * b;
        dat[i].old_ra = old_ra + i * b;
        MYCALLOC(dat[i].ra_v, 1);
        dat[i].g = g;
        if (sr_db->n >= (i + 1) * b) {
            dat[i].n = b;
        } else {
            dat[i].n = sr_db->n - i * b;
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
    scg_ra_v_clean(ra_v);
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
        free(dat[i].ra_v->a);
        free(dat[i].ra_v);
    }

    uint64_t n_r, n_m, n_u, n;
    n_r = n_m = n_u = 0;
    for (n = 0; n < sr_db->n; ++n)
        n_r += sr_db->a[n].n > 0;
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

