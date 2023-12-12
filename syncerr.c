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
 * 25/07/23 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>

#include "kvec.h"
#include "kstring.h"
#include "kthread.h"

#include "syncmer.h"
#include "syncasm.h"
#include "levdist.h"
#include "graph.h"
#include "misc.h"

#undef DEBUG_FIND_SYNCERR
#undef DEBUG_SYNCMER_CORRECTION
#undef DEBUG_SYNCMER_EC_DETAIL
#undef DEBUG_TRACE_ED_ALIGNMENT

typedef kvec_t(uint64_t) kvec64_t;
typedef kvec_t(uint32_t) kvec32_t;

#define EC_FAILURE 0
#define EC_SUCCESS 1
#define EC_AMBISNQ 2 // ambiguous syncmer path
#define EC_AMBISEQ 3 // ambigouous sequence (ambiguous syncmer path for sure)

#ifdef DEBUG_SYNCMER_CORRECTION
static const char *const EC_STATUS[] = {"FAILURE", "SUCCESS", "AMBISNQ", "AMBISEQ"};
#endif

typedef struct {
    int status, n_path, edist, s_edist;
    kstring_t c_seq, opt_seq;
    kvec64_t c_path, opt_path;
} dfs_info_t;

typedef struct {
    wf_config_t *conf;
    kstring_t *seq, *c_seq;
    kvec64_t *c_kmer;
    kvec32_t *c_mpos;
    dfs_info_t *dfs;
    long stats[11]; // tail_err ec_status[4] middle_err ec_status[4] overlap_err
} ec_cached_t;

typedef struct {
    sr_db_t *sr_db;
    scg_t *g;
    double max_edist;
    ec_cached_t *cache;
    FILE *fo;
} ec_shared_t;

pthread_mutex_t mutex;

static void dfs_info_destroy(dfs_info_t *dfs)
{
    if (!dfs) return;
    free(dfs->c_seq.s);
    free(dfs->opt_seq.s);
    free(dfs->c_path.a);
    free(dfs->opt_path.a);
    free(dfs);
}

static void dfs_info_reset(dfs_info_t *dfs)
{
    if (!dfs) return;
    dfs->status = EC_FAILURE;
    dfs->n_path = 0;
    dfs->edist = INT_MAX;
    dfs->s_edist = INT_MAX;
    dfs->c_seq.l = 0;
    dfs->opt_seq.l = 0;
    dfs->c_path.n = 0;
    dfs->opt_path.n = 0;
}

static inline int kputsn_rev(const char *p, int l, kstring_t *s)
{
    if (s->l + l + 1 >= s->m) {
        s->m = s->l + l + 2;
        kroundup32(s->m);
        MYREALLOC(s->s, s->m);
    }
    int i, j;
    for (i = l - 1, j = s->l; i >= 0; --i, ++j) s->s[j] = seq_nt4_comp_table[(int) p[i]];
    s->l += l;
    s->s[s->l] = 0;
    return l;
}

static inline int kputsn1(const char *p, int l, kstring_t *s, int r)
{
    return r? kputsn_rev(p, l, s) : kputsn(p, l, s);
}

static int kvcmp(kvec64_t *a, kvec64_t *b)
{
    if (a->n != b->n)
        return 1;
    size_t i, n;
    for (i = 0, n = a->n; i < n; ++i)
        if (a->a[i] != b->a[i])
            return 1;
    return 0;
}

#define MAX_DFS_PATH 10000

static void dfs_search(asmg_t *asmg, syncmer_t *scms, dfs_info_t *dfs_info, uint64_t sink, wf_config_t *conf)
{
    // bounded by MAX_DFS_PATH
    if (dfs_info->n_path >= MAX_DFS_PATH)
        return;

    uint64_t i, source, arc_n, w, l0, n0, d0;
    int l_seq, ls, s0, t_end0, q_end0, score;
    asmg_arc_t *arc, *a;
    kstring_t *c_seq;
    kvec64_t *c_path;
    char *k_seq;
    wf_diag1_t *wf_diag1;

    c_seq = &dfs_info->c_seq;
    l0 = c_seq->l;
    c_path = &dfs_info->c_path;
    n0 = c_path->n;
    source = c_path->a[n0-1];
    arc = asmg_arc_a(asmg, source);
    arc_n = asmg_arc_n(asmg, source);
    t_end0 = conf->t_end;
    q_end0 = conf->q_end;
    s0 = conf->score;
    d0 = conf->wf_diag->n;
    // make a copy of DP matrix
    MYMALLOC(wf_diag1, d0);
    memcpy(wf_diag1, conf->wf_diag->a, sizeof(wf_diag1_t) * d0);

    // process each arc
    for (i = 0; i < arc_n; ++i) {
        a = &arc[i];
        if (a->del)
            continue;
        w = a->w;
        ls = a->ls;
        l_seq = asmg->vtx[w>>1].len;
        k_seq = asmg->vtx[w>>1].seq;

        // push sink id
        kv_push(uint64_t, *c_path, w);

        // push sequence
        if (w&1)
            kputsn_rev(k_seq, l_seq - ls, c_seq);
        else
            kputsn(&k_seq[ls], l_seq - ls, c_seq);

        conf->qs = c_seq->s;
        conf->ql = c_seq->l;
        // conf->t_end and conf->q_end will be zero if not aligned
        wf_ed_core(conf);

#ifdef DEBUG_SYNCMER_EC_DETAIL
        fprintf(stderr, "[DEBUG_SYNCMER_EC_DETAIL::%s] u%lu%c u%lu%c L_SEQ=%d LS=%d EXT=%d\n", 
                __func__, source>>1, "+-"[source&1], w>>1, "+-"[w&1], l_seq, ls, l_seq - ls);
        fprintf(stderr, "[DEBUG_SYNCMER_EC_DETAIL::%s] ED_ALIGNMENT: ED=%d BW=%d tL=%d t_EN=%d qL=%d q_EN=%d\n",
                __func__, conf->score, conf->bw, conf->tl, conf->t_end, conf->ql, conf->q_end);
        fprintf(stderr, "[DEBUG_SYNCMER_EC_DETAIL::%s] ERROR CORRECTION RESULT - ORIGINAL  SEQ: %.*s\n",
                __func__, conf->tl, conf->ts);
        fprintf(stderr, "[DEBUG_SYNCMER_EC_DETAIL::%s] ERROR CORRECTION RESULT - CORRECTED SEQ: %.*s\n",
                __func__, (int) conf->ql, conf->qs);
#endif

        // calculate real score considering target sequence clipping
        score = conf->score + conf->tl - conf->t_end;
        if (score <= conf->bw && (sink == UINT64_MAX || sink == w)) {
#ifdef DEBUG_SYNCMER_EC_DETAIL
            fprintf(stderr, "[DEBUG_SYNCMER_EC_DETAIL::%s] ERROR CORRECTION RESULT - OPTIMUM ALIGNMENT: score=%d\n",
                    __func__, score);
#endif
            // this is an optimum alignment
            // also update optimum path info when score == edist 
            // to get the complete synmcer list for trailing errors
#ifdef DEBUG_SYNCMER_CORRECTION
            // reset edist
            // this is necessary only ifdef DEBUG_SYNCMER_CORRECTION
            // as edist could be updated before this
            if (!dfs_info->status)
                dfs_info->edist = INT_MAX;
#endif
            dfs_info->status = EC_SUCCESS; // update EC_STATUS
            if (score <= dfs_info->edist) {
                if (conf->t_end > t_end0) // otherwise only an extension
                    dfs_info->s_edist = dfs_info->edist;
                dfs_info->edist = score;
                // for trailing errors the last syncmer may need to be removed
                // because the syncmer sequence is only partially mapped
                if (sink == UINT64_MAX && conf->q_end < conf->ql)
                    --c_path->n;
                // if the secondary alignment is as good as the optimum
                if (dfs_info->edist == dfs_info->s_edist) {
                    // need to compare to the query sequence to decide
                    // if the alignment sequence is ambiguous
                    if (conf->q_end != dfs_info->opt_seq.l || 
                            strncmp(c_seq->s, dfs_info->opt_seq.s, conf->q_end))
                        dfs_info->status = EC_AMBISEQ;
                    // need to compare to the optimum path to decide
                    // if syncmer path is ambiguous
                    if (dfs_info->status == EC_SUCCESS && 
                            kvcmp(c_path, &dfs_info->opt_path))
                        dfs_info->status = EC_AMBISNQ;
                }
                dfs_info->opt_seq.l = 0;
                kputsn(c_seq->s, conf->q_end, &dfs_info->opt_seq);
                kv_copy(uint64_t, dfs_info->opt_path, *c_path);
            } else if (score < dfs_info->s_edist) {
                dfs_info->s_edist = score;
            }
        }
#ifdef DEBUG_SYNCMER_CORRECTION
        else if (!dfs_info->status) { // in debug mode this will be copied if EC_FAILURE
            if (score <= dfs_info->edist) {
                dfs_info->edist = score;
                dfs_info->opt_seq.l = 0;
                // the alignment may fail at this point
                // copy the entire sequence instead of the mapped piece
                kputsn(c_seq->s, c_seq->l, &dfs_info->opt_seq);
                kv_copy(uint64_t, dfs_info->opt_path, *c_path);
            }
        }
#endif
        
        if (conf->score <= conf->bw &&                     // stop dfs if edit distance exceeds threshold
                conf->ql - l_seq <= conf->tl + conf->bw && // stop dfs if target and query sequence overlap exceeds threshold
                ((sink != UINT64_MAX && sink != w) ||      // stop dfs if target end reached and sink reached for middle errors
                 conf->t_end < conf->tl))                  // stop dfs if target end reached for trailing errors
            // further dfs
            dfs_search(asmg, scms, dfs_info, sink, conf);
        else
            // increase path counter when a dfs search stopped
            dfs_info->n_path++;

        // reset conf and dfs_info
        c_path->n = n0;
        c_seq->l = l0;
        conf->t_end = t_end0;
        conf->q_end = q_end0;
        conf->score = s0;
        conf->wf_diag->n = d0;
        memcpy(conf->wf_diag->a, wf_diag1, sizeof(wf_diag1_t) * d0);
    }
    
    free(wf_diag1);
}

int error_correction_by_graph_path_search(asmg_t *asmg, syncmer_t *scms, uint64_t source, uint64_t sink, wf_config_t *conf, dfs_info_t *dfs_info) {
    
    if (conf->tl < 0) return 0;
    
    dfs_info_reset(dfs_info);
    kv_push(uint64_t, dfs_info->c_path, source);
    dfs_search(asmg, scms, dfs_info, sink, conf);

#ifdef DEBUG_SYNCMER_CORRECTION
    if (dfs_info->n_path >= MAX_DFS_PATH)
        fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::WARN::%s] exceeds DFS path limit %d - error correction might not be accurate\n",
                __func__, MAX_DFS_PATH);

    // print debug information even when error correction failed
    fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] ERROR CORRECTION RESULT - %s; N_PATH: %d; EDIST: %d [%.6f]; S_EDIST: %d [%.6f]\n",
            __func__, EC_STATUS[dfs_info->status], dfs_info->n_path, dfs_info->edist, 
            (double) dfs_info->edist/conf->tl, dfs_info->s_edist, (double) dfs_info->s_edist/conf->tl);
    fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] ERROR CORRECTION RESULT - OPT PATH:", __func__);
    size_t i;
    for (i = 0; i < dfs_info->opt_path.n; ++i) fprintf(stderr, " u%lu%c", dfs_info->opt_path.a[i]>>1, "+-"[dfs_info->opt_path.a[i]&1]);
    fputc('\n', stderr);
#ifdef DEBUG_TRACE_ED_ALIGNMENT
    if (dfs_info->opt_seq.l > 0) {
        int score = 0, t_end = 0, q_end = 0, n_cigar = 0;
        // keep in mind that when error correction failed the bandwidth parameter does not apply
        uint32_t *cigar = wf_ed(conf->tl, conf->ts, dfs_info->opt_seq.l, dfs_info->opt_seq.s, 1, dfs_info->status? conf->bw : -1, &score, &t_end, &q_end, &n_cigar);
        fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] ERROR CORRECTION RESULT - ", __func__);
        wf_print_cigar(cigar, n_cigar, stderr);
        fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] ED_ALIGNMENT: ED=%d BW=%d tL=%d t_EN=%d qL=%d q_EN=%d\n",
                __func__, score, dfs_info->status? conf->bw : -1, conf->tl, t_end, (int) dfs_info->opt_seq.l, q_end);
        wf_print_alignment(conf->ts, conf->tl, dfs_info->opt_seq.s, dfs_info->opt_seq.l, cigar, n_cigar, 0, stderr);
        free(cigar);
    }
#else
    fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] ERROR CORRECTION RESULT - ORIGINAL  SEQ: %.*s\n",
            __func__, conf->tl, conf->ts);
    fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] ERROR CORRECTION RESULT - CORRECTED SEQ: %.*s\n",
            __func__, (int) dfs_info->opt_seq.l, dfs_info->opt_seq.s);
#endif
#endif

    return dfs_info->status;
}

#define MIN_ERR_SEQ_LEN 10
#define MIN_ERR_BASE 6

#define MASK_ONE 0xFFFFFFFFFFFFFFFEULL

static void read_error_correction_analysis_thread(void *_data, long i, int tid) // kt_for() callback
{
    ec_shared_t *shared = (ec_shared_t *) _data;
    sr_t *sr;
    asmg_t *asmg;
    syncmer_t *scms;
    uint32_t *m_pos, beg_pos, end_pos;
    uint64_t *k_mer, beg_utg, end_utg;
    int32_t j, l, r, n, n_scm, kmer_size, beg, end;
    int updated;
    double max_edist;
    long *stats;
    FILE *fo;
    wf_config_t *conf;
    dfs_info_t *dfs;
    kstring_t *seq, *c_seq;
    kvec64_t *c_kmer;
    kvec32_t *c_mpos;

    asmg = shared->g->utg_asmg;
    scms = shared->g->scm_db->a;
    kmer_size = shared->sr_db->k;
    sr = &shared->sr_db->a[i];
    k_mer = sr->k_mer;
    m_pos = sr->m_pos;
    n_scm = sr->n;
    max_edist = shared->max_edist;
    fo = shared->fo;
    stats = shared->cache[tid].stats;
    conf = shared->cache[tid].conf;
    dfs = shared->cache[tid].dfs;
    seq = shared->cache[tid].seq;
    c_seq = shared->cache[tid].c_seq;
    c_kmer = shared->cache[tid].c_kmer;
    c_mpos = shared->cache[tid].c_mpos;

#ifdef DEBUG_SYNCMER_CORRECTION
    // read information for debugging
    fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] %lu %s %u %u:", __func__, 
            sr->sid, sr->sname, sr->hoco_l, sr->n);
    for (j = 0; j < n_scm; ++j) {
        if (j > 0) { // print arc coverage
            fprintf(stderr, " [%u]", asmg_arc(asmg, 
                        (k_mer[j-1]&MASK_ONE) | (m_pos[j-1]&1),
                        (k_mer[j]&MASK_ONE) | (m_pos[j]&1))->cov);
        }
        fprintf(stderr, " (%d u%lu%c %u)", j, k_mer[j]>>1, "+-"[m_pos[j]&1], scms[k_mer[j]>>1].cov);
    }
    fputc('\n', stderr);
    // position of syncmers on read
    fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] %lu %s %u %u:", __func__,
            sr->sid, sr->sname, sr->hoco_l, sr->n);
    for (j = 0; j < n_scm; ++j) fprintf(stderr, " %u", m_pos[j]>>1);
    fputc('\n', stderr);
#endif

    // find error syncmer blocks and do correction
    int err_c1;
    // reset preallocated sequence and path for corrections
    c_seq->l = 0;
    c_kmer->n = 0;
    c_mpos->n = 0;
    updated = 1;
#ifdef DEBUG_SYNCMER_CORRECTION
    int err_block = 0, edist = 0;
#endif
    beg = -1; // beg is a bad syncmer
    while (1) {
        beg_pos = beg < 1? 0 : ((m_pos[beg - 1] >> 1) + kmer_size);
        beg_pos += MIN_ERR_SEQ_LEN;

        // find the rigth boundary - end is a good syncmer
        // the distance between beg and end syncmer is no less than MIN_ERR_SEQ_LEN
        // with this restriction a good syncmer could be included for error correction
        for (end = beg + 1; end < n_scm; ++end)
            if (!scms[k_mer[end]>>1].del && !(k_mer[end]&1) && (m_pos[end] >> 1) >= beg_pos) // good syncmer
                break;

#ifdef DEBUG_SYNCMER_CORRECTION
        fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] RID %lu ERROR BLOCK %d: %d %d\n", __func__, sr->sid, ++err_block, beg, end);
#endif

        // do error correction for reads with at least one syncmer anchor
        if (beg >= 0 || end < n_scm) {
            // do graph search for error correction
            // with restrictions in distance and source and sink utg nodes
            // here assembly utgs are syncmer singletons, i.e., utg_id is scm_id
            // always guarantee the source utg exist
            if (beg < 0) {
                // no left boundary
                // start from end instead
                // also do reverse complementary
                beg = end; // now beg is a good syncmer
                beg_utg = (k_mer[beg]&MASK_ONE) | (!(m_pos[beg]&1));
                beg_pos = 0;
                end_utg = UINT64_MAX;
                l = m_pos[beg] >> 1;
                r = 1;
            } else {
                // beg >= 0
                // in fact beg > 0
                --beg; // now beg is a good syncmer
                beg_utg = (k_mer[beg]&MASK_ONE) | (m_pos[beg]&1);
                // beg is a correct syncmer
                // start from the next base after this sycnmer
                beg_pos = (m_pos[beg]>>1) + kmer_size;
                if (end >= n_scm) {
                    // no right boundary
                    end_utg = UINT64_MAX;
                    l = (int) sr->hoco_l - beg_pos;
                } else {
                    // bounded by correct syncmer
                    end_utg = (k_mer[end]&MASK_ONE) | (m_pos[end]&1);
                    l = (int) (m_pos[end]>>1) - beg_pos;
                }
                r = 0;
            }
            
            // the error sequence to correct
            assert(l >= 0);
            ks_resize(seq, l);
            get_kmer_dna_seq(sr->hoco_s, beg_pos, l, r, seq->s);
            seq->l = l;
            if (l >= MIN_ERR_SEQ_LEN) {
                // do error correction
                // edit distance configuration
                conf->ts = seq->s;
                conf->tl = seq->l;
                conf->score = 0;
                conf->qs = 0;
                conf->ql = 0;
                conf->is_ext = 1;
                conf->bw = ceil(l * max_edist);
                // TODO is this good?
                if (conf->bw < MIN_ERR_BASE)
                    conf->bw = MIN_ERR_BASE;
                if (!conf->wf_diag)
                    MYCALLOC(conf->wf_diag, 1);
                if (conf->wf_diag->m < l * 4) {
                    conf->wf_diag->m = l * 4;
                    MYREALLOC(conf->wf_diag->a, l * 4);
                }
                conf->wf_diag->n = 1;
                conf->wf_diag->a[0].d =  0;
                conf->wf_diag->a[0].k = -1;

                err_c1 = error_correction_by_graph_path_search(asmg, scms, beg_utg, end_utg, conf, dfs);

#ifdef DEBUG_SYNCMER_CORRECTION
                fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] RID %lu ERROR BLOCK %d [%s] (PATH = %d): %u %d (u%lu%c u%lu%c) %.*s\n",
                        __func__, sr->sid, err_block, EC_STATUS[err_c1], dfs->n_path, 
                        beg_pos, l, beg_utg>>1, "+-"[beg_utg&1], end_utg>>1, "+-"[end_utg&1], l, seq->s);
#endif
                if (err_c1)
                    assert(beg_utg == dfs->opt_path.a[0] && (end_utg == UINT64_MAX || end_utg == dfs->opt_path.a[dfs->opt_path.n-1]));

                if (end_utg == UINT64_MAX) {
                    ++stats[0];
                    ++stats[1 + err_c1];
                } else {
                    ++stats[5];
                    ++stats[6 + err_c1];
                }
            } else {
                err_c1 = EC_FAILURE;
                ++stats[10];
#ifdef DEBUG_SYNCMER_CORRECTION
                fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::WARN::%s] RID %lu ERROR BLOCK %d FLANKING SYNCMER OVERLAPPED: %u %d (%lu%c %lu%c)\n",
                        __func__, sr->sid, err_block, beg_pos, l, beg_utg>>1, "+-"[beg_utg&1], end_utg>>1, "+-"[end_utg&1]);
                print_aligned_syncmers_on_seq(sr, kmer_size, 0, UINT32_MAX, stderr);
#endif
            }

            // update c_path
            if (err_c1 == EC_SUCCESS) {
                n = dfs->opt_path.n; // n > 0
                // append the corrected syncmer list
                if (r) {
                    // end_utg == UINT64_MAX
                    for (j = n - 1; j > 0; --j) {
                        kv_push(uint64_t, *c_kmer, (dfs->opt_path.a[j]&MASK_ONE) | 1);
                        kv_push(uint32_t, *c_mpos, UINT32_MAX ^ (dfs->opt_path.a[j]&1));
                    }
                } else {
                    for (j = 1; j < n - 1; ++j) {
                        kv_push(uint64_t, *c_kmer, (dfs->opt_path.a[j]&MASK_ONE) | 1);
                        kv_push(uint32_t, *c_mpos, MASK_ONE | (dfs->opt_path.a[j]&1));
                    }
                    if (end_utg == UINT64_MAX && n > 1) {
                        // j == n-1
                        kv_push(uint64_t, *c_kmer, (dfs->opt_path.a[j]&MASK_ONE) | 1);
                        kv_push(uint32_t, *c_mpos, MASK_ONE | (dfs->opt_path.a[j]&1));
                    }
                }
            } else {
                // append the original sycnmer list
                if (r) {
                    kv_pushn(uint64_t, *c_kmer, k_mer, beg);
                    kv_pushn(uint32_t, *c_mpos, m_pos, beg);
                } else if (beg + 1 < n_scm) {
                    kv_pushn(uint64_t, *c_kmer, &k_mer[beg + 1], end - beg - 1);
                    kv_pushn(uint32_t, *c_mpos, &m_pos[beg + 1], end - beg - 1);
                }
            }
            // update c_seq if necessary
            if (fo) {
                if (err_c1 == EC_SUCCESS || err_c1 == EC_AMBISNQ) {
                    // append the corrected sequence
                    if (r)
                        kputsn_rev(dfs->opt_seq.s, dfs->opt_seq.l, c_seq);
                    else
                        kputsn(dfs->opt_seq.s, dfs->opt_seq.l, c_seq);
                } else {
                    // append the original sequence
                    if (r)
                        kputsn_rev(seq->s, seq->l, c_seq);
                    else
                        kputsn(seq->s, seq->l, c_seq);
                }
            }

#ifdef DEBUG_SYNCMER_CORRECTION
            if (err_c1 == EC_SUCCESS || err_c1 == EC_AMBISNQ)
                edist += dfs->edist;
#endif
        } else {
            // no good syncmers on the read
#ifdef DEBUG_SYNCMER_CORRECTION
            pthread_mutex_lock(&mutex);
            fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] %s RID %lu no correct syncmers L=%u N=%u:", __func__, sr->sname, sr->sid, sr->hoco_l, sr->n);
            for (j = 0; j < n_scm; ++j) fprintf(stderr, " %u", scms[k_mer[j]>>1].cov);
            fputc('\n', stderr);
            pthread_mutex_unlock(&mutex);
#endif
            updated = 0;
        }

        // end is always a good syncmer
        // find next bad syncmer as beg
        for (beg = end + 1; beg < n_scm; ++beg)
            if (scms[k_mer[beg]>>1].del || (k_mer[end]&1))
                // bad syncmer
                break;
        
        if (beg > n_scm) break;

        // append kmer and mpos between [end, beg-1]
        kv_pushn(uint64_t, *c_kmer, &k_mer[end], beg - end);
        kv_pushn(uint32_t, *c_mpos, &m_pos[end], beg - end);

        // append sequence between [end, beg-1]
        if (fo) {
            beg_pos = m_pos[end] >> 1;
            end_pos = m_pos[beg-1] >> 1;
            l = end_pos + kmer_size - beg_pos;
            ks_resize(c_seq, c_seq->l + l);
            get_kmer_dna_seq(sr->hoco_s, beg_pos, l, 0, &c_seq->s[c_seq->l]);
            c_seq->l += l;
        }
    }
    
    if (updated) {
        // here is to update syncmer seq for the read
        // TODO need to consider data consistency
        size_t s, n_c = c_kmer->n;
        MYREALLOC(sr->k_mer, n_c);
        memcpy(sr->k_mer, c_kmer->a, sizeof(uint64_t) * n_c);
        MYREALLOC(sr->m_pos, n_c);
        memcpy(sr->m_pos, c_mpos->a, sizeof(uint32_t) * n_c);
        MYREALLOC(sr->s_mer, n_c);
        for (s = 0; s < n_c; ++s)
            sr->s_mer[s] = scms[sr->k_mer[s]>>1].s;
        sr->n = n_c;
    }

    if (fo) {
        if (!updated) {
            l = sr->hoco_l;
            ks_resize(c_seq, l);
            get_kmer_dna_seq(sr->hoco_s, 0, l, 0, c_seq->s);
            c_seq->l = l;
        }
        pthread_mutex_lock(&mutex);
        fprintf(fo, ">%s\n%.*s\n", sr->sname, (int) c_seq->l, c_seq->s);
        pthread_mutex_unlock(&mutex);
    }

#ifdef DEBUG_SYNCMER_CORRECTION
    // read information for debugging
    fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] ERROR CORRECTION RESULT - %lu %s %u %u:", __func__,
            sr->sid, sr->sname, sr->hoco_l, sr->n);
    n_scm = sr->n;
    k_mer = sr->k_mer;
    m_pos = sr->m_pos;
    for (j = 0; j < n_scm; ++j) {
        if (j > 0) { // print arc coverage
            fprintf(stderr, " [%u]", asmg_arc(asmg,
                        (k_mer[j-1]&MASK_ONE) | (m_pos[j-1]&1),
                        (k_mer[j]&MASK_ONE) | (m_pos[j]&1))->cov);
        }
        fprintf(stderr, " (%d u%lu%c %u)", j, k_mer[j]>>1, "+-"[m_pos[j]&1], scms[k_mer[j]>>1].cov);
    }
    fputc('\n', stderr);
    // position of syncmers on read
    fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] ERROR CORRECTION RESULT - %lu %s %u %u:", __func__,
            sr->sid, sr->sname, sr->hoco_l, sr->n);
    for (j = 0; j < n_scm; ++j) fprintf(stderr, " %u", m_pos[j]>>1);
    fputc('\n', stderr);

    if (fo) {
        int score = 0, t_end = 0, q_end = 0, n_cigar = 0;
        l = sr->hoco_l;
        ks_resize(seq, l);
        get_kmer_dna_seq(sr->hoco_s, 0, l, 0, seq->s);
        seq->l = l;
        uint32_t *cigar = wf_ed(seq->l, seq->s, c_seq->l, c_seq->s, 1, -1, &score, &t_end, &q_end, &n_cigar);
        fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] ERROR CORRECTION RESULT - ", __func__);
        wf_print_cigar(cigar, n_cigar, stderr);
        fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] ED_ALIGNMENT: ED=%d BW=%d tL=%d t_EN=%d qL=%d q_EN=%d\n",
                __func__, score, -1, (int) seq->l, t_end, (int) c_seq->l, q_end);
        wf_print_alignment(seq->s, seq->l, c_seq->s, c_seq->l, cigar, n_cigar, 0, stderr);
        free(cigar);
        if (score != edist)
            fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::WARN::%s] RID %lu INCONSISTENT ERROR CORRECTION RESULT - EC=%d, EDIST=%d\n", 
                    __func__, sr->sid, score, edist);
    }
#endif

    return;
}

// find and mark candidate syncmers for error correction
// error syncmer condidates will be marked by syncmer->del
// syncmers with a coverage below err_mer_c will be directly marked
// syncmers with a coverage above max_err_c will never be marked
// syncmers with a coverage fall into the range (err_mer_c, max_err_c) will be checked
// and will be markeed if all the incoming and outgoing arcs are unreliable
// arcs with a coverage below err_arc_c will be considered unreliable
// arcs with a coverage below err_arc_f * MIN(V, W) will also be considered unreliable
// where V and W is the syncmer coverage of the source and target vertex respectively
int64_t find_error_syncmers(scg_t *g, uint32_t err_mer_c, uint32_t max_err_c, uint32_t err_arc_c, double max_arc_f, int del_err)
{
    asmg_t *asmg;
    syncmer_t *scm;
    asmg_arc_t *a;
    size_t i, j, k, n_scm, n_arc;
    uint32_t nv, nw;
    int b[2];

    asmg = g->utg_asmg;
    n_scm = g->scm_db->n;
    scm = g->scm_db->a;

    for (i = 0; i < n_scm; ++i) {
        if (scm[i].del || scm[i].cov >= max_err_c)
            continue;
        if (scm[i].cov < err_mer_c) {
            scm[i].del = 1;
            continue;
        }
        nv = scm[i].cov;
        b[0] = b[1] = -1;
        for (k = 0; k < 2; ++k) {
            a = asmg_arc_a(asmg, i << 1 | k);
            n_arc = asmg_arc_n1(asmg, i << 1 | k);
            if (n_arc) b[k] = 0;
            else continue;
            n_arc = asmg_arc_n(asmg, i << 1 | k);
            for (j = 0; j < n_arc; ++j) {
                if (a[j].del) continue;
                nw = scm[a[j].w >> 1].cov;
                if (a[j].cov >= err_arc_c && a[j].cov >= MIN(nv, nw) * max_arc_f) {
                    // good arc found
                    b[k] = 1;
                    break;
                }
            }
        }
        if (!b[0] || !b[1]) scm[i].del = 1;
    }

    int64_t n_err = 0;
    uint32_t max_c = 0;
    for (i = 0; i < n_scm; ++i) {
        if (scm[i].del) {
            if (scm[i].cov > max_c)
                max_c = scm[i].cov;
            ++n_err;
#ifdef DEBUG_FIND_SYNCERR
            fprintf(stderr, "[DEBUG_FIND_SYNCERR::%s] error syncmer: %lu %u\n", __func__, i, scm[i].cov);
            if (scm[i].cov > err_mer_c) {
                for (k = 0; k < 2; ++k) {
                    a = asmg_arc_a(asmg, i << 1 | k);
                    n_arc = asmg_arc_n(asmg, i << 1 | k);
                    for (j = 0; j < n_arc; ++j) {
                        if (a[j].del || a[j].comp) continue;
                        fprintf(stderr, "[DEBUG_FIND_SYNCERR::%s] error syncmer arc list: "
                                "%lu %lu%c -> %lu%c v_cov = %u w_cov = %u arc_cov = %u\n", 
                                __func__, j, 
                                a[j].v>>1, "+-"[a[j].v&1], 
                                a[j].w>>1, "+-"[a[j].w&1], 
                                scm[a[j].v>>1].cov, scm[a[j].w>>1].cov, a[j].cov);
                    }
                }
            }
#endif
        }
    }
    
    if (del_err) {
        for (i = 0; i < n_scm; ++i)
            if (scm[i].del)
                asmg_vtx_del(asmg, i, 1);
    }

    fprintf(stderr, "[M::%s] error syncmer candidates: num = %ld, max_c = %u\n", __func__, n_err, max_c);

    return n_err;
}

/***
static int arc_cov_cmpfunc(const void *a, const void *b)
{
    uint64_t x, y;
    x = ((asmg_arc_t *) a)->cov;
    y = ((asmg_arc_t *) b)->cov;
    return (x < y) - (x > y);
}
**/

static void update_syncmer_db(sr_db_t *sr_db, syncmer_db_t *scm_db)
{
    uint64_t i, j, k, n, m;
    syncmer_t *scms;
    uint64_t *k_mer, sid;
    uint32_t *m_pos, *c_cov;
    // clean scm_db
    free(scm_db->c); scm_db->c = 0;
    free(scm_db->h); scm_db->h = 0;
    scms = scm_db->a;
    for (i = 0, n = scm_db->n; i < n; ++i)
        scms[i].cov = 0;
    // collect syncmer coverage
    for (i = 0, n = sr_db->n; i < n; ++i) {
        k_mer = sr_db->a[i].k_mer;
        for (j = 0, m = sr_db->a[i].n; j < m; ++j)
            ++scms[k_mer[j]>>1].cov;
    }
    // reallocte k_mer position array
    for (i = 0, n = scm_db->n; i < n; ++i) {
        free(scms[i].m_pos);
        MYMALLOC(scms[i].m_pos, scms[i].cov);
        // set cov to zero which will be recalculated
        scms[i].cov = 0;
    }
    // collect syncmer positions on reads
    MYCALLOC(c_cov, scm_db->n);
    for (i = 0, n = sr_db->n; i < n; ++i) {
        sid = sr_db->a[i].sid;
        k_mer = sr_db->a[i].k_mer;
        m_pos = sr_db->a[i].m_pos;
        for (j = 0, m = sr_db->a[i].n; j < m; ++j) {
            k = k_mer[j] >> 1;
            scms[k].m_pos[scms[k].cov++] = (sid << 32) | (j << 1) | (m_pos[j] & 1);
            if(!(m_pos[j] & 1)) ++c_cov[k];
        }
    }
    // mark syncmers with no coverage as deleted
    // this is necessay as there due to the corner case 
    // for a good syncmer, all copies were error-corrected to others
    // meanwhile other syncmers were error-corrected to become it
    // FIXME the simple solution here is to delete this syncmer
    for (i = 0, n = scm_db->n; i < n; ++i)
        scms[i].del = !c_cov[i];
    free(c_cov);
}

// read error correction in hoco space by aligning to the syncmer graph
// the syncmer consensus sequences need to be presented in input syncmer graph
// max_edist = 0.02 max edit distance to correct a read
void read_error_correction(sr_db_t *sr_db, scg_t *g, double max_edist, uint32_t err_mer_c, uint32_t max_err_c, 
        uint32_t err_arc_c, double max_arc_f, int n_threads, FILE *fo, int verbose)
{
    if (n_threads <= 0)
        n_threads = 1;

#ifdef DEBUG_SYNCMER_CORRECTION
    if (n_threads > 1) {
        n_threads = 1;
        fprintf(stderr, "[DEBUG_SYNCMER_CORRECTION::%s] set thread number to one for debugging mode\n", __func__);
    }
#endif

    double realtime1 = realtime();
    double cputime1 = cputime();

    // mark errors in syncmer database
    find_error_syncmers(g, err_mer_c, max_err_c, err_arc_c, max_arc_f, 1);
    
    // sort arcs by coverage for each source vertex
    // DFS then will visit the high coverage paths first
    // this will hopefully increase the chance of error correction for complex regions capped by MAX_DFS_PATH
    // TODO does this really help?
    /***
    uint64_t i, n, na;
    asmg_arc_t *a;
    asmg_t *asmg;
    asmg = g->utg_asmg;
    for (i = 0, n = 2 * asmg->n_vtx; i < n; ++i) {
        na = asmg_arc_n(asmg, i);
        if (na) {
            a = asmg_arc_a(asmg, i);
            qsort(a, na, sizeof(asmg_arc_t), arc_cov_cmpfunc);
        }
    }
    **/

    int t;
    // shared data
    ec_shared_t shared;

    ec_cached_t *cached;
    MYCALLOC(cached, n_threads);
    for (t = 0; t < n_threads; ++t) {
        MYCALLOC(cached[t].conf, 1);
        MYCALLOC(cached[t].seq, 1);
        MYCALLOC(cached[t].c_seq, 1);
        MYCALLOC(cached[t].c_kmer, 1);
        MYCALLOC(cached[t].c_mpos, 1);
        MYCALLOC(cached[t].dfs, 1);
        MYBZERO(cached[t].stats, 11);
    }

    shared.sr_db = sr_db;
    shared.g = g;
    shared.max_edist = max_edist;
    shared.cache = cached;
    shared.fo = fo;
    
    if (pthread_mutex_init(&mutex, NULL) != 0) {
        fprintf(stderr, "[E::%s] pthread mutex init failed\n", __func__);
        goto do_clean;
    }
    kt_for(n_threads, read_error_correction_analysis_thread, &shared, sr_db->n);
    pthread_mutex_destroy(&mutex);
    
    /***
    uint64_t rids[] = {0, 9, 15, 21, 31}; //{9007, 58477, 59602, 89548, 129398, 192, 15778, 24013, 25430, 44358, 102676, 103790, 129066, 140012};
    for (int r = 0; r < sizeof(rids) / sizeof(*rids); ++r)
        read_error_correction_analysis_thread(&shared, rids[r], 0);
    **/

    int j;
    long *stats;
    stats = cached[0].stats;
    for (t = 1; t < n_threads; ++t)
        for (j = 0; j < 11; ++j)
            stats[j] += cached[t].stats[j];
    
    // update syncmer database
    update_syncmer_db(sr_db, g->scm_db);

    // print error correction summary results
    fprintf(stderr, "[M::%s] Error Correction Summary Results\n", __func__);
    fprintf(stderr, "[M::%s] total number of error blocks : %ld\n", __func__, stats[0] + stats[5] + stats[10]);
    fprintf(stderr, "[M::%s]                - uncorrected : %ld\n", __func__, stats[1] + stats[6]);
    fprintf(stderr, "[M::%s]                  - corrected : %ld\n", __func__, stats[2] + stats[7]);
    fprintf(stderr, "[M::%s]             - ambiguous seqs : %ld\n", __func__, stats[3] + stats[8]);
    fprintf(stderr, "[M::%s]             - ambiguous path : %ld\n", __func__, stats[4] + stats[9]);
    // more verbose information
    if (verbose) {
        fprintf(stderr, "[M::%s] error blocks in the tail end : %ld\n", __func__, stats[0]);
        fprintf(stderr, "[M::%s]                - uncorrected : %ld\n", __func__, stats[1]);
        fprintf(stderr, "[M::%s]                  - corrected : %ld\n", __func__, stats[2]);
        fprintf(stderr, "[M::%s]             - ambiguous seqs : %ld\n", __func__, stats[3]);
        fprintf(stderr, "[M::%s]             - ambiguous path : %ld\n", __func__, stats[4]);
        fprintf(stderr, "[M::%s]   error blocks in the middle : %ld\n", __func__, stats[5]);
        fprintf(stderr, "[M::%s]                - uncorrected : %ld\n", __func__, stats[6]);
        fprintf(stderr, "[M::%s]                  - corrected : %ld\n", __func__, stats[7]);
        fprintf(stderr, "[M::%s]             - ambiguous seqs : %ld\n", __func__, stats[8]);
        fprintf(stderr, "[M::%s]             - ambiguous path : %ld\n", __func__, stats[9]);
        fprintf(stderr, "[M::%s]      error blocks overlapped : %ld\n", __func__, stats[10]);
        fprintf(stderr, "[M::%s]   error correction  CPU time : %.3f sec\n", __func__, cputime() - cputime1);
        fprintf(stderr, "[M::%s]   error correction real time : %.3f sec\n", __func__, realtime() - realtime1);
    }

do_clean:
    for (t = 0; t < n_threads; ++t) {
        // do free
        wf_config_destroy(cached[t].conf, 0);
        dfs_info_destroy(cached[t].dfs);
        free(cached[t].seq->s);
        free(cached[t].seq);
        free(cached[t].c_seq->s);
        free(cached[t].c_seq);
        free(cached[t].c_kmer->a);
        free(cached[t].c_kmer);
        free(cached[t].c_mpos->a);
        free(cached[t].c_mpos);
    }
    free(cached);
}



