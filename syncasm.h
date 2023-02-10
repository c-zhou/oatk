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

#ifndef __SYNCASM_H__
#define __SYNCASM_H__

#include "syncmer.h"
#include "graph.h"

#define scg_utg_t asmg_vtx_t

#define scm_utg_pos(u) ((uint64_t)((u)&0xFFFFFFFFFULL))
#define scm_utg_uid(u) ((uint64_t)(((u)>>36)&0x3FFFFFFFFFFULL))
#define scm_utg_scm(u) ((uint64_t)((u)>>78))
#define scm_utg_rev(u) ((uint64_t)(((u)>>78)&1))
#define scm_utg_sid(u) ((uint64_t)((u)>>79))
#define scm_utg_n(g, s) (((g)->idx_u[(s)+1])-((g)->idx_u[(s)]))

typedef struct {
    /*** syncmer graph ***/
    uint64_t m_scm, n_scm;
    syncmer_t *scm;
    void *h_scm; // scm to id map
    /*** unitig graph ***/
    asmg_t *utg_asmg;
    /*** auxilliary information ***/
    // utg scm index scm_id[49]|scm_rev[1]|utg_id[42]|utg_pos[36]
    uint128_t *scm_u;
    // utg scm starting positions
    // of size n_scm + 1
    // idx_u[n_scm] - idx_u[0] is the size of scm_u
    uint128_t **idx_u;
} scg_t;

// read alignment to syncmer graph
typedef struct {
    uint64_t uid; // utg id << 1 | strand
    uint64_t u_beg;
    uint64_t u_end; // inclusive
    uint32_t s_beg;
    uint32_t s_end; // inclusive
} ra_frg_t;

typedef struct {
    uint64_t sid; // seq id
    uint32_t n; // number fragment
    ra_frg_t *a; // alignment fragment
    double s; // alignment score s.(1/n)
} scg_ra_t;

typedef struct {size_t n, m; scg_ra_t *a;} scg_ra_v;

#ifdef __cplusplus
extern "C" {
#endif

void scg_destroy(scg_t *g);
scg_t *make_syncmer_graph(sr_v *sr, uint32_t min_k_cov, uint32_t min_a_cov);
void scg_arc_coverage(scg_t *scg, sr_v *sr);
void process_mergeable_unitigs(scg_t *g);
int scg_is_empty(scg_t *scg);
void scg_stat(scg_t *scg, FILE *fo, uint64_t *stats);
void scg_consensus(sr_v *sr, scg_t *scg, int w, int save_seq, FILE *fo);
void scg_print_unitig_syncmer_list(scg_t *g, FILE *fo);
int64_t scg_syncmer_consensus(sr_v *sr, syncmer_t *scm, int rev, int64_t beg, int w, kstring_t *c_seq);
int64_t scg_unitig_consensus(sr_v *sr, scg_utg_t *utg, syncmer_t *scm, int w, kstring_t *c_seq);
int scg_multiplex(scg_t *g, scg_ra_v *ra_v, double r_thresh);
void scg_demultiplex(scg_t *g);
void scg_ra_v_destroy(scg_ra_v *ra_v);
void scg_read_alignment(sr_v *sr, scg_ra_v *ra_v, scg_t *g, int n_threads, int for_unzip);
void scg_ra_print(scg_ra_t *ra, FILE *fo);
void scg_rv_print(scg_ra_v *rv, FILE *fo);
#ifdef __cplusplus
}
#endif

#endif

