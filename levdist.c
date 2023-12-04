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
 * 24/07/23 - Chenxi Zhou: Created                                               *
 *                                                                               *
 * The core part was adapted from Heng Li's implemetation of the                 *
 * Landau-Vishkin/Myers86 algorithm to calculate edit distance                   *
 *                                                                               *
 * https://github.com/lh3/lv89                                                   *
 *                                                                               *
 * Modified to adapt to graph extension alignment                                *
 *                                                                               *
 *********************************************************************************/

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "kstring.h"
#include "levdist.h"
#include "misc.h"

void wf_tb_destroy(wf_tb_t *wf_tb)
{
    if (wf_tb) {
        size_t i;
        for (i = 0; i < wf_tb->n; ++i)
            free(wf_tb->a[i].a);
        free(wf_tb->a);
    }
    free(wf_tb);
}

void wf_config_destroy(wf_config_t *conf, int clean_seq)
{
    if (clean_seq) {
        free(conf->ts);
        free(conf->qs);
    }
    if (conf->wf_diag) {
        free(conf->wf_diag->a);
        free(conf->wf_diag);
    }
    wf_tb_destroy(conf->wf_tb);
    free(conf->cigar);
    free(conf);
}

// Extend a diagonal along exact matches. This is a bottleneck and could be made faster with padding.
static inline int32_t wf_extend(int32_t tl, const char *ts, int32_t ql, const char *qs, const wf_diag1_t *p)
{
    int32_t k = p->k;
    const char *ts_, *qs_;
    uint64_t cmp = 0;
    int32_t max_k = (ql - p->d < tl? ql - p->d : tl) - 1;
    ts_ = ts + 1;
    qs_ = qs + p->d + 1;
    while (k + 7 < max_k) {
        uint64_t x = *(uint64_t*)(ts_ + k); // warning: unaligned memory access
        uint64_t y = *(uint64_t*)(qs_ + k);
        cmp = x ^ y;
        if (cmp == 0) k += 8;
        else break;
    }
    if (cmp)
        k += __builtin_ctzl(cmp) >> 3; // on x86, this is done via the BSR instruction: https://www.felixcloutier.com/x86/bsr
    else if (k + 7 >= max_k)
        while (k < max_k && *(ts_ + k) == *(qs_ + k)) // use this for generic CPUs. It is slightly faster than the unoptimized version
            ++k;
    return k;
}

// filter diagonals with a fixed bandwidth
static inline void wf_prune_bw(int32_t tl, int32_t ql, int32_t is_ext, int32_t bw, int32_t *st, int32_t *en, const wf_diag1_t *b)
{
    int32_t min_d, max_d, s = *st, e = *en;
    if (is_ext) {
        min_d = -bw, max_d = bw;
    } else {
        min_d = ql < tl? ql - tl - bw : tl - ql - bw;
        max_d = tl > ql? tl - ql + bw : ql - tl + bw;
    }
    min_d = min_d > -tl? min_d : -tl;
    max_d = max_d >  ql? max_d :  ql;
    while (b[s].d < min_d) ++s;
    while (b[e - 1].d > max_d) --e;
    *st = s, *en = e;
}

static void wf_tb_add(wf_tb_t *tb, int32_t m)
{
    wf_tb1_t *p;
    if (tb->n == tb->m) {
        tb->m = tb->m + (tb->m>>1) + 4;
        MYREALLOC(tb->a, tb->m);
    }
    p = &tb->a[tb->n++];
    p->n = m;
    MYCALLOC(p->a, (m + 31) / 32);
}

static void wf_cigar_push1(wf_cigar_t *c, int32_t op, int32_t len)
{
    if (c->n && op == (c->cigar[c->n-1]&0xf)) {
        c->cigar[c->n-1] += len<<4;
    } else {
        if (c->n == c->m) {
            c->m = c->m + (c->m>>1) + 4;
            MYREALLOC(c->cigar, c->m);
        }
        c->cigar[c->n++] = len<<4 | op;
    }
}

static void wf_cigar_push(wf_cigar_t *c, int32_t n_cigar, const uint32_t *cigar)
{
    if (n_cigar == 0) return;
    wf_cigar_push1(c, cigar[0]&0xf, cigar[0]>>4);
    if (c->n + n_cigar - 1 > c->m) {
        c->m = c->n + n_cigar - 1;
        kroundup32(c->m);
        MYREALLOC(c->cigar, c->m);
    }
    memcpy(&c->cigar[c->n], &cigar[1], sizeof(uint32_t) * (n_cigar - 1));
    c->n += n_cigar - 1;
}

/*
 * Basic LV89
 */
static int wf_step_basic(wf_tb_t *tb, int32_t is_ext, int32_t bw, int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t n, wf_diag1_t *a, int32_t *t_end, int32_t *q_end)
{
    int32_t j, st = 0, en = n + 2;
    wf_diag1_t *b = a + n + 2; // temporary array

    // wfa_extend
    *t_end = *q_end = -1;
    for (j = 0; j < n; ++j) {
        wf_diag1_t *p = &a[j];
        int32_t k = p->k;
        if (k >= tl || k + p->d >= ql) continue;
        k = wf_extend(tl, ts, ql, qs, p);
        if (k + p->d == ql - 1 || k == tl - 1) {
            if (is_ext || (k + p->d == ql - 1 && k == tl - 1)) {
                *t_end = k, *q_end = k + p->d;
                return -1;
            }
        }
        p->k = k;
    }

    // wfa_next
    b[0].d = a[0].d - 1;
    b[0].p = 1;
    b[0].k = a[0].k + 1;
    b[1].d = a[0].d;
    b[1].p =  n == 1 || a[0].k > a[1].k? 0 : 1;
    b[1].k = (n == 1 || a[0].k > a[1].k? a[0].k : a[1].k) + 1;
    if (tb) { // with traceback
        for (j = 1; j < n - 1; ++j) {
            int32_t k = a[j-1].k, p = -1;
            p = k > a[j].k + 1? p : 0;
            k = k > a[j].k + 1? k : a[j].k + 1;
            p = k > a[j+1].k + 1? p : 1;
            k = k > a[j+1].k + 1? k : a[j+1].k + 1;
            b[j+1].d = a[j].d, b[j+1].k = k, b[j+1].p = p;
        }
    } else { // simpler without traceback
        for (j = 1; j < n - 1; ++j) {
            int32_t k = a[j-1].k;
            k = k > a[j].k + 1? k : a[j].k + 1;
            k = k > a[j+1].k + 1? k : a[j+1].k + 1;
            b[j+1].d = a[j].d, b[j+1].k = k;
        }
    }
    if (n >= 2) {
        b[n].d = a[n-1].d;
        b[n].p = a[n-2].k > a[n-1].k + 1? -1 : 0;
        b[n].k = a[n-2].k > a[n-1].k + 1? a[n-2].k : a[n-1].k + 1;
    }
    b[n+1].d = a[n-1].d + 1;
    b[n+1].p = -1;
    b[n+1].k = a[n-1].k;

    if (bw < 0 || n < bw + bw + 1) {
        if (b[0].d < -tl) ++st;
        if (b[n+1].d > ql) --en;
    } else wf_prune_bw(tl, ql, is_ext, bw, &st, &en, b);
    if (tb) { // keep traceback information
        wf_tb1_t *q;
        wf_tb_add(tb, n + 2);
        q = &tb->a[tb->n - 1];
        q->d0 = b[st].d;
        for (j = 0; j < en - st; ++j)
            q->a[j>>5] |= (uint64_t)(b[j+st].p + 1) << (j&0x1f)*2;
    }
    memcpy(a, &b[st], (en - st) * sizeof(*a));
    return en - st;
}

// traceback
static uint32_t *wf_traceback(int32_t t_end, const char *ts, int32_t q_end, const char *qs, wf_tb_t *tb, int32_t *n_cigar)
{
    wf_cigar_t cigar = {0,0,0};
    int32_t i = q_end, k = t_end, s = tb->n - 1;
    for (;;) {
        int32_t k0 = k, j, pre;
        while (i >= 0 && k >= 0 && qs[i] == ts[k])
            --i, --k;
        if (k0 - k > 0)    
            wf_cigar_push1(&cigar, 7, k0 - k);
        if (i < 0 || k < 0) break;
        assert(s >= 0);
        j = i - k - tb->a[s].d0;
        assert(j < tb->a[s].n);
        pre = (int32_t)(tb->a[s].a[j>>5] >> (j&0x1f)*2 & 0x3) - 1;
        if (pre == 0) {
            wf_cigar_push1(&cigar, 8, 1);
            --i, --k;
        } else if (pre < 0) {
            wf_cigar_push1(&cigar, 1, 1);
            --i;
        } else {
            wf_cigar_push1(&cigar, 2, 1);
            --k;
        }
        --s;
    }
    if (i >= 0) wf_cigar_push1(&cigar, 1, i + 1);
    else if (k >= 0) wf_cigar_push1(&cigar, 2, k + 1);
    for (i = 0; i < cigar.n>>1; ++i) {
        uint32_t t = cigar.cigar[i];
        cigar.cigar[i] = cigar.cigar[cigar.n - i - 1];
        cigar.cigar[cigar.n - i - 1] = t;
    }
    *n_cigar = cigar.n;
    return cigar.cigar;
}

void wf_ed_core(wf_config_t *conf)
{
    char *ts, *qs;
    int32_t s, n, tl, ql, t_end, q_end, is_ext, bw;
    wf_diag1_t *a;
    wf_tb_t *tb;

    ts = conf->ts;
    qs = conf->qs;
    tl = conf->tl;
    ql = conf->ql;
    t_end = conf->t_end;
    q_end = conf->q_end;

    s = conf->score;
    n = conf->wf_diag->n;
    if (conf->wf_diag->m < 2 * (tl + ql + 2)) {
        conf->wf_diag->m = 2 * (tl + ql + 2);
        MYREALLOC(conf->wf_diag->a, conf->wf_diag->m);
    }
    a = conf->wf_diag->a;
    n = conf->wf_diag->n;
    tb = conf->wf_tb;
    
    is_ext = conf->is_ext;
    bw = conf->bw;

    int na = n;
    while (1) {
        na = wf_step_basic(tb, is_ext, bw, tl, ts, ql, qs, n, a, &t_end, &q_end);
        if (na < 0) break;
        ++s;
        n = na;
        if (bw >= 0 && s > bw) break;
    }
    
    if (tb) {
        free(conf->cigar);
        conf->cigar = wf_traceback(t_end, ts, q_end, qs, tb, &conf->n_cigar);
    }
    conf->t_end = t_end + 1;
    conf->q_end = q_end + 1;
    conf->score = s;
    conf->wf_diag->a = a;
    conf->wf_diag->n = n;
}

uint32_t *wf_ed(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t is_ext, int32_t bw, int32_t *score, int32_t *t_endl, int32_t *q_endl, int32_t *n_cigar)
{
    int32_t s = 0, n = 1, t_end = -1, q_end = -1, i;
    wf_diag1_t *a;
    uint32_t *cigar = 0;
    wf_tb_t tb = {0, 0, 0};
    assert(tl > 0 && ql > 0);
    MYMALLOC(a, 2 * (tl + ql + 2)); // without CIGAR, this would be all the memory needed
    a[0].d = 0, a[0].k = -1;
    while (1) {
        n = wf_step_basic(n_cigar? &tb : 0, is_ext, bw, tl, ts, ql, qs, n, a, &t_end, &q_end);
        if (n < 0) break;
        ++s;
    }
    free(a);
    if (n_cigar) { // generate CIGAR
        cigar = wf_traceback(t_end, ts, q_end, qs, &tb, n_cigar);
        for (i = 0; i < tb.n; ++i) free(tb.a[i].a);
        free(tb.a);
    }
    *score = s, *t_endl = t_end + 1, *q_endl = q_end + 1;
    return cigar;
}

void wf_print_cigar(uint32_t *cigar, int32_t n_cigar, FILE *fo)
{
    if (!cigar || !n_cigar)
        return;

    int32_t i;
    fprintf(fo, "CIGAR_STR [%d]: ", n_cigar);
    for (i = 0; i < n_cigar; ++i) {
        fprintf(fo, "%d%c", cigar[i]>>4, "MIDNSHP=XB"[cigar[i]&0xf]);
    }
    fputc('\n', fo);
}

static inline void kputcn_(int c, int l, kstring_t *s)
{
    if (s->l + l > s->m) {
        s->m = s->l + l;
        kroundup32(s->m);
        MYREALLOC(s->s, s->m);
    }
    while(l-- > 0) s->s[s->l++] = c;
}

void wf_print_alignment(char *ts, int tl, char *qs, int ql, uint32_t *cigar, int32_t n_cigar, int lwd, FILE *fo)
{
    if (!cigar || !n_cigar)
        return;
    
    int32_t i, ti, qi;
    uint32_t cl, l;
    uint8_t c;
    kstring_t t_seq = {0, 0, 0}, q_seq = {0, 0, 0}, a_seq = {0, 0, 0};
    ti = qi = l = 0;
    for (i = 0; i < n_cigar; ++i) {
        c = "MIDNSHP=XB"[cigar[i]&0xf];
        cl = cigar[i]>>4;
        switch (c) {
            case 'I':
                kputcn_('-', cl, &t_seq);
                kputsn_(&qs[qi], cl, &q_seq);
                kputcn_(' ', cl, &a_seq);
                qi += cl;
                break;
            case 'D':
                kputsn_(&ts[ti], cl, &t_seq);
                kputcn_('-', cl, &q_seq);
                kputcn_(' ', cl, &a_seq);
                ti += cl;
                break;
            case 'M':
            case '=':
            case 'X':
                kputsn_(&ts[ti], cl, &t_seq);
                if (c == 'M')
                    for (l = 0; l < cl; ++l)
                        kputc_(ts[ti+l] == qs[qi+l]? '|' : '*', &a_seq);
                else
                    kputcn_(c == '='? '|' : '*', cl, &a_seq);
                kputsn_(&qs[qi], cl, &q_seq);
                ti += cl;
                qi += cl;
                break;
            default:
                fprintf(stderr, "[E::%s] unsupported cigar element '%c'\n", __func__, c);
                goto do_clean;
        }
    }

    if (lwd <= 0) {
        fprintf(fo, "target: %.*s\n", (int) t_seq.l, t_seq.s);
        fprintf(fo, "        %.*s\n", (int) a_seq.l, a_seq.s);
        fprintf(fo, " query: %.*s\n", (int) q_seq.l, q_seq.s);
    } else {
        uint32_t j, t_beg, t_end, q_beg, q_end;
        int b = 0;
        l = t_seq.l;
        while (l > 0) { l /= 10; ++b; }
        kstring_t pad = {0, 0, 0};
        kputcn_(' ', b + b + 14, &pad);
        l = 0;
        t_beg = q_beg = 0;
        while (l < t_seq.l) {
            cl = l;
            l += lwd;
            if (l > t_seq.l)
                l = t_seq.l;
            t_end = t_beg + l - cl;
            q_end = q_beg + l - cl;
            for (j = cl; j < l; ++j) {
                t_end -= (t_seq.s[j] == '-');
                q_end -= (q_seq.s[j] == '-');
            }
            fprintf(fo, "target [%*.u - %-*.u]: %.*s\n", b, t_beg + 1, b, t_end, (int) (l - cl), &t_seq.s[cl]);
            fprintf(fo, "%.*s%.*s\n", (int) pad.l, pad.s, (int) (l - cl), &a_seq.s[cl]);
            fprintf(fo, " query [%*.u - %-*.u]: %.*s\n", b, q_beg + 1, b, q_end, (int) (l - cl), &q_seq.s[cl]);
            fputc('\n', fo);
            t_beg = t_end;
            q_beg = q_end;
        }
        free(pad.s);
    }

do_clean:
    free(t_seq.s); free(q_seq.s); free(a_seq.s);
}

#if defined LEVDIST_TEST_NAIVE // compile with gcc -DLEVDIST_TEST_NAIVE levdist.c
int main(int argc, char *argv[])
{
    char tsDefault[] = "AATGCTCTCATGACATATGAGATAGATACATAGAGACAGATATAGATACACACAGAGATATATGACGTCTGTATGCTCTCTCTCATAGATATACTCTGTAGACTGTCATATACATGCAGAAAAA";
    char qsDefault[] = "CGCTCTCATGACANATGAGATAGATACATAGAGNCAGATATAGATACACACAGTTT";
    int is_ext = 1;

    char *ts, *qs;
    if (argc > 2) {
        ts = argv[1];
        qs = argv[2];
        if (argc > 3) is_ext = atoi(argv[3]);
    } else {
        ts = tsDefault;
        qs = qsDefault;
    }

    int tl = strlen(ts), ql = strlen(qs);

    wf_config_t *conf;
    MYCALLOC(conf, 1);
    conf->ts = ts;
    conf->qs = qs;
    conf->tl = tl;
    conf->ql = ql;
    conf->is_ext = is_ext;
    conf->bw = -1;
    MYCALLOC(conf->wf_diag, 1);
    MYCALLOC(conf->wf_tb, 1);

    conf->wf_diag->n = 1;
    conf->wf_diag->m = 2 * (tl + ql + 2);
    MYMALLOC(conf->wf_diag->a, conf->wf_diag->m);
    conf->wf_diag->a[0].d =  0;
    conf->wf_diag->a[0].k = -1;
    
    conf->ql = ql;
    conf->tl = tl;
    wf_ed_core(conf);
    fprintf(stdout, "[M::%s] TARGET %s\n", __func__, ts);
    fprintf(stdout, "[M::%s] QUERY  %s\n", __func__, qs);
    if (conf->bw >= 0 && conf->score > conf->bw) {
        fprintf(stdout, "[M::%s] no alignment with edit distance threshold BW=%d\n", __func__, conf->bw);
    } else {
        fprintf(stdout, "[M::%s] ED=%d tL=%d t_EN=%d qL=%d q_EN=%d\n", __func__, conf->score, conf->tl, conf->t_end, conf->ql, conf->q_end);
        wf_print_cigar(conf->cigar, conf->n_cigar, stdout);
        wf_print_alignment(ts, conf->tl, qs, conf->ql, conf->cigar, conf->n_cigar, 100, stdout);
    }

    wf_config_destroy(conf, 0);

    return 0;
}
#elif defined LEVDIST_TEST_STEP // compile with gcc -DLEVDIST_TEST_STEP levdist.c kopen.c -lz
#include <time.h>
#include <zlib.h>
#include "ketopt.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)
int main(int argc, char *argv[])
{
    argc--; argv++;
    if (argc < 2) {
        fprintf(stderr, "[E::%s] at least need two input files\n", __func__);
        return 1;
    }

    gzFile fp1, fp2;
    kseq_t *ks1, *ks2;

    fp1 = gzopen(argv[0], "r");
    fp2 = gzopen(argv[1], "r");
    assert(fp1 && fp2);
    ks1 = kseq_init(fp1);
    ks2 = kseq_init(fp2);
    
    int s = 0, n = 0;
    if (argc > 2)
        s = atoi(argv[2]);
    if (argc > 3)
        n = atoi(argv[3]);

    if (!s) s = 200; // max step size
    if (!n) n = 20;  // number trials

    srand(time(NULL));

    while (kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0) {
        wf_config_t *conf;
        MYCALLOC(conf, 1);
        conf->ts = ks1->seq.s;
        conf->qs = ks2->seq.s;
        conf->tl = ks1->seq.l;
        conf->ql = ks2->seq.l;
        conf->is_ext = 1;
        conf->bw = -1;
        MYCALLOC(conf->wf_diag, 1);
        MYCALLOC(conf->wf_tb, 1);
        // conf->wf_tb = 0;

        conf->score = 0;
        conf->wf_diag->n = 1;
        conf->wf_diag->m = 2 * (conf->tl + conf->ql + 2);
        MYMALLOC(conf->wf_diag->a, conf->wf_diag->m);
        conf->wf_diag->a[0].d =  0;
        conf->wf_diag->a[0].k = -1;
        
        wf_ed_core(conf);
        fprintf(stdout, "[M::%s] ED=%d tL=%d t_EN=%d qL=%d q_EN=%d\n", __func__, conf->score, conf->tl, conf->t_end, conf->ql, conf->q_end);
        if (conf->wf_tb) {
            wf_print_cigar(conf->cigar, conf->n_cigar, stdout);
            wf_print_alignment(conf->ts, conf->tl, conf->qs, conf->ql, conf->cigar, conf->n_cigar, 100, stdout);        
        }

        int i, q_endl, step, score = conf->score, t_end = conf->t_end, q_end = conf->q_end;
        // no alignment backtrack
        wf_tb_destroy(conf->wf_tb);
        conf->wf_tb = 0;
        for (i = 0; i < n; ++i) {
            // reset DP matrix
            conf->score = 0;
            conf->wf_diag->n = 1;
            conf->wf_diag->a[0].d =  0;
            conf->wf_diag->a[0].k = -1;
            /***
            conf->ql = ks2->seq.l;
            wf_ed_core(conf);
            **/
            // randomly do alignment in stepwise mode
            q_endl = 0;
            step = 0;
            while (q_endl < ks2->seq.l) {
                q_endl += rand() % s;
                if (q_endl > ks2->seq.l)
                    q_endl = ks2->seq.l;
                conf->ql = q_endl;
                wf_ed_core(conf);
                ++step;
            }
            // end of trial
            if (conf->score != score || conf->t_end != t_end || conf->q_end != q_end) {
                fprintf(stdout, "[M::%s] inconsistent alignment in trial %d: ED=%d tL=%d t_EN=%d qL=%d q_EN=%d\n", 
                        __func__, i, conf->score, conf->tl, conf->t_end, conf->ql, conf->q_end);
            } else {
                fprintf(stdout, "[M::%s] successful trial %d in %d steps\n", __func__, i, step);
            }
        }

        wf_config_destroy(conf, 0);
    }
    kseq_destroy(ks1);
    kseq_destroy(ks2);
    gzclose(fp1);
    gzclose(fp2);

    return 0;
}
#endif

