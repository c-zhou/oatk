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
 * 11/08/22 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdio.h>
#include <assert.h>
#include <pthread.h>
#include <math.h>

#include "misc.h"
#include "khashl.h"
#include "kvec.h"
#include "syncmer.h"
#include "MurmurHash3.h"

#undef DEBUG_KMER_EXTRACTION
#undef DEBUG_S_KMER_GROUP

unsigned char seq_nt4_table[256] = {
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

char char_nt4_table[4] = {'A', 'C', 'G', 'T'};

// complementary seq of four base pairs packed in a uint8_t
// (A,C,G,T) -> (0,1,2,3)
// AAAA=0b00000000 -> comp_seq8_table[0]=0b11111111
/***
unsigned char *make_comp_seq8_table()
{
    unsigned char comp_seq8_table[256];
    int comp_seq_table[4] = {3, 2, 1, 0};
    uint8_t i, j, c;
    for (i = 0; ; ++i) {
        c = 0;
        for (j = 0; j < 4; ++j)
            c = c << 2 | comp_seq_table[(i >> j * 2) & 3];
        comp_seq8_table[i] = c;
        if (i == 255) break;
    }
}
**/
static unsigned char comp_seq8_table[256] = {
    255, 191, 127, 63,  239, 175, 111, 47,  223, 159, 95, 31,  207, 143, 79, 15,
    251, 187, 123, 59,  235, 171, 107, 43,  219, 155, 91, 27,  203, 139, 75, 11,
    247, 183, 119, 55,  231, 167, 103, 39,  215, 151, 87, 23,  199, 135, 71, 7,
    243, 179, 115, 51,  227, 163, 99,  35,  211, 147, 83, 19,  195, 131, 67, 3,
    254, 190, 126, 62,  238, 174, 110, 46,  222, 158, 94, 30,  206, 142, 78, 14,
    250, 186, 122, 58,  234, 170, 106, 42,  218, 154, 90, 26,  202, 138, 74, 10,
    246, 182, 118, 54,  230, 166, 102, 38,  214, 150, 86, 22,  198, 134, 70, 6,
    242, 178, 114, 50,  226, 162, 98,  34,  210, 146, 82, 18,  194, 130, 66, 2,
    253, 189, 125, 61,  237, 173, 109, 45,  221, 157, 93, 29,  205, 141, 77, 13,
    249, 185, 121, 57,  233, 169, 105, 41,  217, 153, 89, 25,  201, 137, 73, 9,
    245, 181, 117, 53,  229, 165, 101, 37,  213, 149, 85, 21,  197, 133, 69, 5,
    241, 177, 113, 49,  225, 161, 97,  33,  209, 145, 81, 17,  193, 129, 65, 1,
    252, 188, 124, 60,  236, 172, 108, 44,  220, 156, 92, 28,  204, 140, 76, 12,
    248, 184, 120, 56,  232, 168, 104, 40,  216, 152, 88, 24,  200, 136, 72, 8,
    244, 180, 116, 52,  228, 164, 100, 36,  212, 148, 84, 20,  196, 132, 68, 4,
    240, 176, 112, 48,  224, 160, 96,  32,  208, 144, 80, 16,  192, 128, 64, 0
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

static uint32_t murmur3_seed = 1234;
static uint8_t lmask[4] = {255, 192, 240, 252};

#ifdef DEBUG_KMER_EXTRACTION
static uint128_t hash128(uint64_t sid, uint8_t *s, uint32_t p, int w)
#else
static uint128_t hash128(uint8_t *s, uint32_t p, int w)
#endif
{
    uint128_t h128;
    uint8_t *key;
    int i, j, rev, res, p0, p1, b, c;

    rev = p & 1, p >>= 1;
    p0 = p, p1 = p + w - 1;

    // assert(p0 >= 0 && p1 < l);

    // res = rev? 6 - p1 % 4 * 2 : p0 % 4 * 2; // shift bits
    res = rev? ((p1&3)^3)<<1 : (p0&3)<<1;
    b = p1 / 4 - p0 / 4 + 1; // number bytes holding syncmer
    MYMALLOC(key, b);
    memcpy(key, s + p0 / 4, b);
    // key[b-1] needs bit shift
    // when the last byte is partially filled
    // if (l < (p1 / 4 + 1) * 4) key[b-1] <<= 8 - l % 4 * 2;
    // do reverse complementary
    if (rev) {
        for (i = 0, j = b - 1; i < j; ++i, --j) {
            c = key[i];
            key[i] = comp_seq8_table[key[j]];
            key[j] = comp_seq8_table[c];
        }
        if (i == j) key[i] = comp_seq8_table[key[i]];
    }
    // do bit shift to align syncmer bytes
    for (i = 0; i < b-1; ++i) {
        key[i] <<= res;
        key[i] |= key[i+1] >> (8-res);
    }
    key[b-1] <<= res;
    // mask the lower bits
    key[b-1] &= lmask[w&3];

    // murmur3 hashing
    MurmurHash3_x64_128(key, (w - 1) / 4 + 1, murmur3_seed, &h128);

#ifdef DEBUG_KMER_EXTRACTION
    fprintf(stderr, "[DEBUG_KMER_EXTRACTION::%s] sid:%lu p0:%d p1:%d rev:%d ", __func__, sid, p0, p1, rev);
    for (i = 0; i < p1 - p0 + 1; ++i)
        // fputc(char_nt4_table[(key[i / 4] >> (3 - i % 4) * 2) & 3], stderr);
        fputc(char_nt4_table[(key[i/4]>>(((i&3)^3)<<1)) & 3], stderr);

    fputc('\n', stderr);
#endif

    free(key);
    
    return h128;
}

typedef struct {
    int n_reads;
    uint64_t *sid;
    char **name;
    char **seq;
    int *len;
    int k, w;
    sr_v sr;
} p_data_t;

static inline int q_next(int i, int q)
{
    ++i;
    return i == q? 0 : i;
}

static void *sr_read_analysis_thread(void *args)
{
    p_data_t *dat = (p_data_t *) args;
    
    int len, k, w; 
    char *seq;
    k = dat->k, w = dat->w;
    assert(k > 0 && k < 32 && w > k);

    int r;
    for (r = 0; r < dat->n_reads; ++r) {
        sr_t sr;
        sr.sid = dat->sid[r];
        sr.sname = dat->name[r];

        seq = dat->seq[r];
        len = dat->len[r];

        kvec_t(uint8_t) hoco_s, ho_rl;
        kvec_t(uint32_t) ho_l_rl, n_nucl, m_pos;
        kvec_t(uint64_t) s_mer;
        kvec_t(uint128_t) k_mer_h;
    
        kv_init(hoco_s);
        kv_init(ho_rl);
        kv_init(ho_l_rl);
        kv_init(n_nucl);
        kv_init(m_pos);
        kv_init(s_mer);
        kv_init(k_mer_h);

        int q = w - k + 1; // buf q size
        uint64_t m, s, z, mz, shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0, 0}, buf_m[q], buf_s[q];
        int i, j, l, c, neq, rl, hoco_l, buf_pos, mz_pos;

        MYBONE(buf_m, q);
        MYBONE(buf_s, q);
        mz = UINT64_MAX; // minimizer
        l = hoco_l = buf_pos = mz_pos = 0;
        for (i = 0; i < len; ++i) {
            c = seq_nt4_table[(uint8_t) seq[i]];
            m = s = UINT64_MAX;
            if ((hoco_l++ & 3) == 0) kv_push(uint8_t, hoco_s, 0);
            if (c < 4) { // not an ambiguous base
                // hoco_s.a[hoco_s.n - 1] = hoco_s.a[hoco_s.n - 1] << 2 | c;
                if (c) hoco_s.a[hoco_s.n - 1] |= c << ((((hoco_l-1)&3)^3)<<1); // 6 - ((hoco_l - 1) % 4) << 1;
                // hpc
                rl = 1;
                if (i + 1 < len && seq_nt4_table[(uint8_t)seq[i + 1]] == c) {
                    for (rl = 2; i + rl < len; ++rl)
                        if (seq_nt4_table[(uint8_t)seq[i + rl]] != c)
                            break;
                    i += rl - 1; // put $i at the end of the current homopolymer run
                }
            
                if (rl > 255)
                    kv_push(uint32_t, ho_l_rl, rl - 1);
                rl = MIN(rl, 256);
                kv_push(uint8_t, ho_rl, rl - 1);
            
                ++l;
                kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
                kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
                if (kmer[0] != kmer[1]) { // skip "symmetric k-mers" as we don't know it strand
                    z = kmer[0] < kmer[1]? 0 : 1; // strand
                    if (l >= k) {
                        m = hash64(kmer[z], mask);
                        s = kmer[z] << 1 | z;
                    }
                }
            } else {
                // ambiguous bases are converted to base 'A'
                // ambiguous bases are not homopolymer compressed
                // hoco_s.a[hoco_s.n - 1] <<= 2;
                kv_push(uint32_t, n_nucl, i);
                l = 0;
            }

            if (buf_pos == mz_pos && mz != UINT64_MAX && l > w) {
                // open syncmer
                z = buf_s[buf_pos] & 1;
                kv_push(uint64_t, s_mer, buf_s[buf_pos]);
                kv_push(uint32_t, m_pos, (hoco_l - w - 1) << 1 | z);
#ifdef DEBUG_KMER_EXTRACTION
                kv_push(uint128_t, k_mer_h, hash128(sr.sid, hoco_s.a, m_pos.a[m_pos.n-1], w));
#else
                kv_push(uint128_t, k_mer_h, hash128(hoco_s.a, m_pos.a[m_pos.n-1], w));
#endif
                // remove syncmers at the same position on a read
                // this is possible as a syncmer could start and end with the same smer
                if (m_pos.n >= 2 && m_pos.a[m_pos.n-1] >> 1 == m_pos.a[m_pos.n-2] >> 1) s_mer.n -= 2, m_pos.n -= 2, k_mer_h.n -= 2;
            }
            
            buf_m[buf_pos] = m;
            buf_s[buf_pos] = s;
            if (m <= mz && m != UINT64_MAX) {
                if (l >= w) {
                    // close syncmer
                    z = s & 1;
                    kv_push(uint64_t, s_mer, s^1);
                    kv_push(uint32_t, m_pos, (hoco_l - w) << 1 | z);
#ifdef DEBUG_KMER_EXTRACTION
                    kv_push(uint128_t, k_mer_h, hash128(sr.sid, hoco_s.a, m_pos.a[m_pos.n-1], w));
#else
                    kv_push(uint128_t, k_mer_h, hash128(hoco_s.a, m_pos.a[m_pos.n-1], w));
#endif
                }
                if (m < mz) mz = m, mz_pos = buf_pos;
            }
            if (m >= mz && buf_pos == mz_pos) {
                // m has been processed in the last 'if' statement if m == mz
                neq = m != mz;
                // update minimizer
                // find the first for identical minimizers
                for (j = buf_pos + 1, mz = UINT64_MAX; j < q; ++j)
                    if (mz > buf_m[j]) mz = buf_m[j], mz_pos = j;
                for (j = 0; j <= buf_pos; ++j)
                    if (mz > buf_m[j]) mz = buf_m[j], mz_pos = j;
                if (neq && ((mz_pos == q_next(buf_pos, q) && mz == m) || mz_pos == buf_pos) && mz != UINT64_MAX && l >= w) {
                    // newly added S-mer is a minimizer
                    // close syncmer
                    z = s & 1;
                    kv_push(uint64_t, s_mer, s^1);
                    kv_push(uint32_t, m_pos, (hoco_l - w) << 1 | z);
#ifdef DEBUG_KMER_EXTRACTION
                    kv_push(uint128_t, k_mer_h, hash128(sr.sid, hoco_s.a, m_pos.a[m_pos.n-1], w));
#else
                    kv_push(uint128_t, k_mer_h, hash128(hoco_s.a, m_pos.a[m_pos.n-1], w));
#endif
                }
            }
            
            buf_pos = q_next(buf_pos, q);
        }

        // for the last open mer
        if (buf_pos == mz_pos && mz != UINT64_MAX && l >= w) { // not (l > w) as l no self increment yet
            // open syncmer
            z = buf_s[buf_pos] & 1;
            kv_push(uint64_t, s_mer, buf_s[buf_pos]);
            kv_push(uint32_t, m_pos, (hoco_l - w) << 1 | z); // not (hoco_l - w - 1) as hoco_l no self increment yet
#ifdef DEBUG_KMER_EXTRACTION
            kv_push(uint128_t, k_mer_h, hash128(sr.sid, hoco_s.a, m_pos.a[m_pos.n-1], w));
#else
            kv_push(uint128_t, k_mer_h, hash128(hoco_s.a, m_pos.a[m_pos.n-1], w));
#endif
            if (m_pos.n >= 2 && m_pos.a[m_pos.n-1] >> 1 == m_pos.a[m_pos.n-2] >> 1) s_mer.n -= 2, m_pos.n -= 2, k_mer_h.n -= 2;
        }
    
        if (hoco_s.n) MYREALLOC(hoco_s.a, hoco_s.n);
        if (ho_rl.n) MYREALLOC(ho_rl.a, ho_rl.n);
        if (ho_l_rl.n) MYREALLOC(ho_l_rl.a, ho_l_rl.n);
        if (n_nucl.n) MYREALLOC(n_nucl.a, n_nucl.n);
        if (m_pos.n) MYREALLOC(m_pos.a, m_pos.n);
        if (s_mer.n) MYREALLOC(s_mer.a, s_mer.n);
        if (k_mer_h.n) MYREALLOC(k_mer_h.a, k_mer_h.n);

        sr.hoco_l = hoco_l;
        sr.hoco_s = hoco_s.a;
        sr.ho_rl = ho_rl.a;
        sr.ho_l_rl = ho_l_rl.a;
        sr.n_nucl = n_nucl.a;
        sr.n = m_pos.n;
        sr.m_pos = m_pos.a;
        sr.s_mer = s_mer.a;
        sr.k_mer_h = k_mer_h.a;
        kv_push(sr_t, dat->sr, sr);

        free(seq);
    }
    
    return NULL;
}

static void do_analysis(p_data_t *dat, pthread_t *threads, int n_threads, sr_v *sr)
{
    int t;
    uint64_t i;
    for (t = 1; t < n_threads; ++t)
        pthread_create(threads + t, NULL, sr_read_analysis_thread, dat + t);

    sr_read_analysis_thread(dat);

    for (t = 1; t < n_threads; ++t)
        pthread_join(threads[t], NULL);

    for (t = 0; t < n_threads; ++t) {
        for (i = 0; i < dat[t].sr.n; ++i)
            kv_push(sr_t, *sr, dat[t].sr.a[i]);
        dat[t].n_reads = 0;
        dat[t].sr.n = 0;
    }

    return;
}

static void sr_read_single_thread(sstream_t *s_stream, sr_v *sr, int k, int w)
{
    p_data_t *dat;
    MYCALLOC(dat, 1);
    dat->n_reads = 1;
    MYMALLOC(dat->sid, 1);
    MYMALLOC(dat->name, 1);
    MYMALLOC(dat->seq, 1);
    MYMALLOC(dat->len, 1);
    dat->k = k;
    dat->w = w;
    kv_init(dat->sr);

    int l;
    uint64_t i;
    i = 0;
    while ((l = sstream_read(s_stream)) >= 0) {
        dat->sid[0] = i++;
        dat->name[0] = strdup(s_stream->s->ks->name.s);
        dat->seq[0] = strdup(s_stream->s->ks->seq.s);
        dat->len[0] = l;
        sr_read_analysis_thread(dat);
    }

    if (sr->n) free(sr->a);
    sr->a = dat->sr.a;
    sr->n = dat->sr.n;
    sr->m = dat->sr.m;

    free(dat->sid);
    free(dat->name);
    free(dat->seq);
    free(dat->len);
    free(dat);

    return;
}

void sr_read(sstream_t *s_stream, sr_v *sr, int k, int w, int n_threads)
{
    if (n_threads == 1) {
        sr_read_single_thread(s_stream, sr, k, w);
        return;
    }

    p_data_t *dat;
    int t;
    uint64_t batch_n;
    
    MYCALLOC(dat, n_threads);
    batch_n = 10000;
    for (t = 0; t < n_threads; ++t) {
        dat[t].n_reads = 0;
        MYMALLOC(dat[t].sid, batch_n);
        MYMALLOC(dat[t].name, batch_n);
        MYMALLOC(dat[t].seq, batch_n);
        MYMALLOC(dat[t].len, batch_n);
        dat[t].k = k;
        dat[t].w = w;
        kv_init(dat[t].sr);
        kv_resize(sr_t, dat[t].sr, batch_n);
    }

    pthread_t threads[n_threads];
    int l;
    uint64_t i, j, n;
    i = j = n = 0;
    while ((l = sstream_read(s_stream)) >= 0) {
        t = j / batch_n;
        n = dat[t].n_reads;
        dat[t].sid[n] = i;
        dat[t].name[n] = strdup(s_stream->s->ks->name.s);
        dat[t].seq[n] = strdup(s_stream->s->ks->seq.s);
        dat[t].len[n] = l;
        dat[t].n_reads++;

        ++i, ++j;
        if (j == batch_n * n_threads) {
            do_analysis(dat, threads, n_threads, sr);
            j = 0;
        }
    }
    if (j > 0) do_analysis(dat, threads, n_threads, sr);

    for (t = 0; t < n_threads; ++t) {
        free(dat[t].sid);
        free(dat[t].name);
        free(dat[t].seq);
        free(dat[t].len);
        kv_destroy(dat[t].sr);
    }
    free(dat);

    return;
}

KHASHL_MAP_INIT(KH_LOCAL, kh_ctab_t, kh_ctab, khint_t, int, kh_hash_uint32, kh_eq_generic)

static int syncmer_s_cmpfunc(const void *a, const void *b)
{
    uint64_t x, y;
    uint128_t s, t;
    x = ((syncmer_t *) a)->s;
    y = ((syncmer_t *) b)->s;
    s = ((syncmer_t *) a)->h;
    t = ((syncmer_t *) b)->h;
    return x == y? (s == t? 0 : (s > t? 1 : -1)) : (x > y? 1 : -1);
}

static int syncmer_h_cmpfunc(const void *a, const void *b)
{
    uint128_t x, y;
    uint64_t s, t;
    x = ((syncmer_t *) a)->h;
    y = ((syncmer_t *) b)->h;
    s = ((syncmer_t *) a)->s;
    t = ((syncmer_t *) b)->s;
    return x == y? (s == t? 0 : (s > t? 1 : -1)) : (x > y? 1 : -1);
}

static int int64_cmpfunc(const void *a, const void *b)
{
    int64_t x, y;
    x = *(int64_t *) a;
    y = *(int64_t *) b;
    return (x > y) - (x < y);
}

#ifdef DEBUG_S_KMER_GROUP
static int uint64_r_cmpfunc(const void *a, const void *b)
{
    uint64_t x, y;
    x = *(uint64_t *) a;
    y = *(uint64_t *) b;
    return (x < y) - (x > y);
}
#endif

static void kh_ctab_put1(kh_ctab_t *ctab, int c)
{
    int absent;
    khint_t k;
    k = kh_ctab_put(ctab, c, &absent);
    if (absent)
        kh_val(ctab, k) = 1;
    else
        ++kh_val(ctab, k);
}

static double kh_ctab_stat(kh_ctab_t *ctab, int *uniq, int *singleton)
{
    khint_t k;
    double sum;
    int s, c, n;
    sum = .0;
    n = c = 0;
    if (kh_size(ctab) == 0)
        goto do_assign;

    for (k = (khint_t) 0; k != kh_end(ctab); ++k) {
        if (kh_exist(ctab, k)) {
            s = kh_key(ctab, k);
            c = kh_val(ctab, k);
            sum += (double) s * c;
            n += c;
        }
    }
    k = kh_ctab_get(ctab, 1);
    if (kh_exist(ctab, k)) c = kh_val(ctab, k);

do_assign:
    if (uniq) *uniq = n;
    if (singleton) *singleton = c;
    
    return sum / n;
}

static int64_t *kh_ctab_cnt(kh_ctab_t *ctab, int max_n)
{
    int64_t *cnt;
    int s, c;
    khint_t k;

    MYCALLOC(cnt, max_n + 1);
    for (k = (khint_t) 0; k != kh_end(ctab); ++k) {
        if (kh_exist(ctab, k)) {
            s = kh_key(ctab, k);
            c = kh_val(ctab, k);
            if (s < max_n)
                cnt[s] = c;
            else
                cnt[max_n] += c;
        }
    }

    return cnt;
}

static void hist_plot(int64_t *hist, uint32_t n, const char *h, FILE *fo)
{
    if (n < 5) return;

    uint32_t i, j, b, p_cnt, cnts[n];
    int32_t c, c_digits;
    double tot_cnt, cnt;
    
    tot_cnt = .0;
    cnts[0] = cnts[1] = cnts[2] = 0;
    for (i = 3; i < n; ++i) {
        cnts[i] = (uint32_t) hist[i];
        tot_cnt += cnts[i];
    }
    tot_cnt *= .99;
    cnt = .0;
    b = 0;
    for (i = 0; i < n; ++i) {
        cnt += cnts[i];
        if (cnt >= tot_cnt) {
            b = i + 1;
            break;
        }
    }
    p_cnt = 0;
    for (i = 0; i < b; ++i)
        if (p_cnt < cnts[i])
            p_cnt = cnts[i];

    c_digits = 0;
    for (i = 0; i < b; ++i) {
        c = hist[i] >> 32;
        int d = c > 0? 0 : 1;
        do {
            c /= 10;
            ++d;
        } while (c != 0);
        if (d > c_digits) c_digits = d;
    }
    if (b < n) ++c_digits;

    double per_dot = MAX(1, p_cnt / 100);
    uint32_t d;
    for (i = 0; i < b; ++i) {
        c = hist[i] >> 32;
        cnt = (uint32_t) hist[i];
        fprintf(fo, "[M::%s] [%s] %*d: ", __func__, h, c_digits, c);
        d = (uint32_t) (cnt / per_dot);
        for (j = 0; j < MIN(d, 100); ++j) fputc('*', fo);
        d = cnt/per_dot > 100? (uint32_t) log10(cnt/per_dot/100) : 0;
        for (j = 0; j < MAX(d, 0); ++j) fputc('+', fo);
        fprintf(fo, " %d\n", (int) cnt);
    }
    if (b < n) {
        c = hist[b-1] >> 32;
        cnt = .0;
        for (i = b; i < n; ++i)
            cnt += (uint32_t) hist[i];
        fprintf(fo, "[M::%s] [%s] >%*d: ", __func__, h, c_digits-1, c);
        d = (uint32_t) (cnt / per_dot);
        for (j = 0; j < MIN(d, 100); ++j) fputc('*', fo);
        d = cnt/per_dot > 100? (uint32_t) log10(cnt/per_dot/100) : 0;
        for (j = 0; j < MAX(d, 0); ++j) fputc('+', fo);
        fprintf(fo, " %d\n", (int) cnt);
    }
}

static void kh_ctab_print(kh_ctab_t *ctab, const char *h, FILE *fo, int more)
{
    khint_t k;
    uint32_t i, n = kh_size(ctab);
    int64_t hist[n];
    int s, c;
    i = 0;
    for (k = (khint_t) 0; k != kh_end(ctab); ++k) {
        if (kh_exist(ctab, k)) {
            s = kh_key(ctab, k);
            c = kh_val(ctab, k);
            hist[i++] = (int64_t) s << 32 | c;
        }
    }

    qsort(hist, n, sizeof(int64_t), int64_cmpfunc);
    hist_plot(hist, n, h, fo);

    if (more > 0) {
        for (i = 0; i < n; ++i)
            fprintf(fo, "[M::%s] [%s CNTS] %ld %d\n", __func__, h, hist[i] >> 32, (int) hist[i]);
    }
}

#define MAX_DEPTH 1000
#define LOWEST_CUT 5

static void ha_hist_line(int c, int x, int exceed, int64_t cnt)
{
    int j;
    if (c >= 0)
        fprintf(stderr, "[M::%s] %5d: ", __func__, c);
    else 
        fprintf(stderr, "[M::%s] %5s: ", __func__, "rest");
    for (j = 0; j < x; ++j) fputc('*', stderr);
    if (exceed) fputc('>', stderr);
    fprintf(stderr, " %lld\n", (long long)cnt);
}

static int ha_analyze_count(int n_cnt, int start_cnt, const int64_t *cnt, int *peak_het, int verbose)
{
    // TODO code review
    const int hist_max = 100;
    int i, start, low_i, max_i, max2_i, max3_i;
    int64_t max, max2, max3, min;

    // determine the start point
    assert(n_cnt > start_cnt);
    *peak_het = -1;
    start = cnt[1] > 0? 1 : 2;

    // find the low point from the left
    low_i = start > start_cnt? start : start_cnt;
    for (i = low_i + 1; i < n_cnt; ++i)
        if (cnt[i] > cnt[i-1])
            break;
    low_i = i - 1;
    if (verbose > 0) fprintf(stderr, "[M::%s] lowest: count[%d] = %ld\n", __func__, low_i, (long)cnt[low_i]);
    if (low_i == n_cnt - 1) return -1; // low coverage

    // find the highest peak
    max_i = low_i + 1, max = cnt[max_i];
    for (i = low_i + 1; i < n_cnt; ++i)
        if (cnt[i] > max)
            max = cnt[i], max_i = i;
    if (verbose > 0) fprintf(stderr, "[M::%s] highest: count[%d] = %ld\n", __func__, max_i, (long)cnt[max_i]);

    // print histogram
    for (i = start; i < n_cnt; ++i) {        
        int x, exceed = 0;
        x = (int)((double)hist_max * cnt[i] / cnt[max_i] + .499);
        if (x > hist_max) exceed = 1, x = hist_max; // may happen if cnt[2] is higher
        if (i > max_i && x == 0) break;
        if (verbose > 0) ha_hist_line(i, x, exceed, cnt[i]);
    }
    {
        int x, exceed = 0;
        int64_t rest = 0;
        for (; i < n_cnt; ++i) rest += cnt[i];
        x = (int)((double)hist_max * rest / cnt[max_i] + .499);
        if (x > hist_max) exceed = 1, x = hist_max;
        if(verbose > 0) ha_hist_line(-1, x, exceed, rest);
    }

    // look for smaller peak on the low end
    max2 = -1; max2_i = -1;
    for (i = max_i - 1; i > low_i; --i)
        if (cnt[i] >= cnt[i-1] && cnt[i] >= cnt[i+1]) 
            if (cnt[i] > max2)
                max2 = cnt[i], max2_i = i;
    if (max2_i > low_i && max2_i < max_i) {
        for (i = max2_i + 1, min = max; i < max_i; ++i)
            if (cnt[i] < min)
                min = cnt[i];
        if (max2 < max * 0.05 || min > max2 * 0.95)
            max2 = -1, max2_i = -1;
    }
    if (max2 > 0)
        if (verbose > 0) fprintf(stderr, "[M::%s] left: count[%d] = %ld\n", __func__, max2_i, (long)cnt[max2_i]);
    else 
        if (verbose > 0) fprintf(stderr, "[M::%s] left: none\n", __func__);

    // look for smaller peak on the high end
    max3 = -1; max3_i = -1;
    for (i = max_i + 1; i < n_cnt - 1; ++i)
        if (cnt[i] >= cnt[i-1] && cnt[i] >= cnt[i+1])
            if (cnt[i] > max3)
                max3 = cnt[i], max3_i = i;
    if (max3_i > max_i) {
        for (i = max_i + 1, min = max; i < max3_i; ++i)
            if (cnt[i] < min)
                min = cnt[i];
        if (max3 < max * 0.05 || min > max3 * 0.95 || max3_i > max_i * 2.5)
            max3 = -1, max3_i = -1;
    }
    if (max3 > 0)
        if (verbose > 0) fprintf(stderr, "[M::%s] right: count[%d] = %ld\n", __func__, max3_i, (long)cnt[max3_i]);
    else
        if (verbose > 0) fprintf(stderr, "[M::%s] right: none\n", __func__);
    if (max3_i > 0) {
        *peak_het = max_i;
        return max3_i;
    } else {
        if (max2_i > 0) *peak_het = max2_i;
        return max_i;
    }
}

void sr_stat(sr_v *sr, sr_stat_t *stats, int w, FILE *fo, int verbose)
{
    size_t i, j, n, m;
    int c, p0, p1;
    uint64_t s;
    uint128_t h;
    kvec_t(syncmer_t) syncmers;
    kh_ctab_t *dist_ctab, *smer_ctab, *kmer_ctab, *smer_k_ctab;

    dist_ctab = kh_ctab_init();
    smer_ctab = kh_ctab_init();
    kmer_ctab = kh_ctab_init();
    smer_k_ctab = kh_ctab_init();
    kv_init(syncmers);
    n = sr->n;
    m = 0;
    for (i = 0; i < n; ++i) {
        sr_t s = sr->a[i];
        m += s.n;
        p0 = p1 = -1;
        for (j = 0; j < s.n; ++j) {
            syncmer_t a = {s.k_mer_h[j], s.s_mer[j], 0, 0, 0};
            kv_push(syncmer_t, syncmers, a);
            p0 = p1;
            p1 = s.m_pos[j] >> 1;
            if (p0 >= 0) kh_ctab_put1(dist_ctab, p1 - p0 - w);
        }
    }

    double dist;
    dist = kh_ctab_stat(dist_ctab, NULL, NULL);
    
    assert(m == syncmers.n);
    if (syncmers.n == 0) {
        fprintf(fo, "[M::%s] empty syncmer collection\n", __func__);
        return;
    }
    
    syncmer_t *a;

    qsort(syncmers.a, syncmers.n, sizeof(syncmer_t), syncmer_s_cmpfunc);
    a = syncmers.a;
    s = a[0].s, c = 1;
    for (i = 1; i < syncmers.n; ++i) {
        if (a[i].s != s) {
            kh_ctab_put1(smer_ctab, c);
            s = a[i].s, c = 1;
        } else {
            ++c;
        }
    }
    kh_ctab_put1(smer_ctab, c);

#ifdef DEBUG_S_KMER_GROUP
    kvec_t(uint64_t) k_cnt;
    kv_init(k_cnt);
    a = syncmers.a;
    s = a[0].s, h = a[0].h, c = 1;
    for (i = 1; i < syncmers.n; ++i) {
        if (a[i].h != h) {
            kv_push(uint64_t, k_cnt, c);
            if (a[i].s != s) {
                qsort(k_cnt.a, k_cnt.n, sizeof(uint64_t), uint64_r_cmpfunc);
                fprintf(stderr, "[DEBUG_S_KMER_GROUP::%s] %lu", __func__, k_cnt.a[0]);
                for (j = 1; j < k_cnt.n; ++j) fprintf(stderr, " %lu", k_cnt.a[j]);
                fprintf(stderr, "\n");
                kv_size(k_cnt) = 0;
                s = a[i].s;
            }
            h = a[i].h, c = 1;
        } else {
            ++c;
        }
    }
    kv_push(uint64_t, k_cnt, c);
    qsort(k_cnt.a, k_cnt.n, sizeof(uint64_t), uint64_r_cmpfunc);
    fprintf(stderr, "[DEBUG_S_KMER_GROUP::%s] %lu", __func__, k_cnt.a[0]);
    for (j = 1; j < k_cnt.n; ++j) fprintf(stderr, " %lu", k_cnt.a[j]);
    fprintf(stderr, "\n");
    kv_destroy(k_cnt);
#endif

    int smer1, smeru, s_peak_hom, s_peak_het;
    int64_t *s_cnts;
    double smera;
    
    smera = kh_ctab_stat(smer_ctab, &smeru, &smer1);
    s_cnts = kh_ctab_cnt(smer_ctab, MAX_DEPTH);
    s_peak_hom = s_peak_het = 0;
    s_peak_hom = ha_analyze_count(MAX_DEPTH+1, LOWEST_CUT, s_cnts, &s_peak_het, verbose-1);

    qsort(syncmers.a, syncmers.n, sizeof(syncmer_t), syncmer_h_cmpfunc);
    a = syncmers.a;
    h = a[0].h, c = 1;
    for (i = 1; i < syncmers.n; ++i) {
        if (a[i].h != h) {
            kh_ctab_put1(kmer_ctab, c);
            h = a[i].h, c = 1;
        } else {
            ++c;
        }
    }
    kh_ctab_put1(kmer_ctab, c);
    
    int kmer1, kmeru, k_peak_hom, k_peak_het;
    int64_t *k_cnts;
    double kmera;
    
    kmera = kh_ctab_stat(kmer_ctab, &kmeru, &kmer1);
    k_cnts = kh_ctab_cnt(kmer_ctab, MAX_DEPTH);
    k_peak_hom = k_peak_het = 0;
    k_peak_hom = ha_analyze_count(MAX_DEPTH+1, LOWEST_CUT, k_cnts, &k_peak_het, verbose-1);

    fprintf(fo, "[M::%s] number syncmers collected: %lu\n", __func__, m);
    fprintf(fo, "[M::%s] number syncmers per read: %.3f\n", __func__, (double) m / n);
    fprintf(fo, "[M::%s] average kmer space: %.3f\n", __func__, dist);
    fprintf(fo, "[M::%s] number uniqe smer: %d; singletons: %d (%.3f%%)\n", __func__, smeru, smer1, (double) smer1 * 100 / smeru);
    fprintf(fo, "[M::%s] average smer count: %.3f\n", __func__, smera);
    fprintf(fo, "[M::%s] smer peak_hom: %d; peak_het: %d\n", __func__, s_peak_hom, s_peak_het);
    fprintf(fo, "[M::%s] number uniqe kmer: %d; singletons: %d (%.3f%%)\n", __func__, kmeru, kmer1, (double) kmer1 * 100 / kmeru);
    fprintf(fo, "[M::%s] average kmer count: %.3f\n", __func__, kmera);
    fprintf(fo, "[M::%s] kmer peak_hom: %d; peak_het: %d\n", __func__, k_peak_hom, k_peak_het);
    
    stats->syncmer_n = m;
    stats->syncmer_per_read = (double) m / n;
    stats->syncmer_avg_dist = dist;
    stats->smer_unique = smeru;
    stats->smer_singleton = smer1;
    stats->smer_avg_cnt = smera;
    stats->smer_peak_hom = s_peak_hom;
    stats->smer_peak_het = s_peak_het;
    stats->kmer_unique = kmeru;
    stats->kmer_singleton = kmer1;
    stats->kmer_avg_cnt = kmera;
    stats->kmer_peak_hom = k_peak_hom;
    stats->kmer_peak_het = k_peak_het;

    if (verbose > 1) {
        kh_ctab_print(dist_ctab, "DIST", fo, verbose - 1);
        kh_ctab_print(smer_ctab, "SMER", fo, verbose - 1);
        kh_ctab_print(kmer_ctab, "KMER", fo, verbose - 1);
    }

    kh_ctab_destroy(dist_ctab);
    kh_ctab_destroy(smer_ctab);
    kh_ctab_destroy(kmer_ctab);
    kh_ctab_destroy(smer_k_ctab);
    free(s_cnts);
    free(k_cnts);
    kv_destroy(syncmers);

    return;
}

int sr_validate(sr_v *sr)
{
    if (sr->n > MAX_RD_NUM) {
        fprintf(stderr, "[E::%s] read number exceeds the limit %llu\n", __func__, MAX_RD_NUM);
        return 1;
    }
    size_t i;
    for (i = 0; i < sr->n; ++i) {
        sr_t *s = &sr->a[i];
        if (s->n > MAX_RD_SCM) {
            fprintf(stderr, "[E::%s] syncmer number (%u) on read exceeds the limit %llu: %s\n", __func__, s->n, MAX_RD_SCM, s->sname);
            return 2;
        }
    }
    return 0;
}

void sr_destroy(sr_t *sr)
{
    if (sr->sname) free(sr->sname);
    if (sr->hoco_s) free(sr->hoco_s);
    if (sr->ho_rl) free(sr->ho_rl);
    if (sr->ho_l_rl) free(sr->ho_l_rl);
    if (sr->n_nucl) free(sr->n_nucl);
    if (sr->s_mer) free(sr->s_mer);
    if (sr->k_mer_h) free(sr->k_mer_h);
    if (sr->m_pos) free(sr->m_pos);
}

void sr_v_destroy(sr_v *sr_v)
{
    size_t i;
    for (i = 0; i < sr_v->n; ++i)
        sr_destroy(&sr_v->a[i]);
    kv_destroy(*sr_v);
}

static void fputs_smer(uint64_t s, int k, FILE *fo)
{
    int i;
    for (i = 0; i < k; ++i)
        fputc(char_nt4_table[(s>>(k-i-1)*2) & 3ULL], fo);
}

static void fputs_kmer(uint8_t *s, uint32_t p, int w, FILE *fo)
{
    int i, j, rev, res, p0, p1, b, c;
    uint8_t *key;

    rev = p & 1, p >>= 1;
    p0 = p, p1 = p + w - 1;

    // res = rev? 6 - p1 % 4 * 2 : p0 % 4 * 2; // shift bits
    res = rev? ((p1&3)^3)<<1 : (p0&3)<<1;
    b = p1 / 4 - p0 / 4 + 1; // number bytes holding syncmer
    MYMALLOC(key, b);
    memcpy(key, s + p0 / 4, b);
    
    // key[b-1] needs bit shift
    // when the last byte is partially filled
    // if (l < (p1 / 4 + 1) * 4) key[b-1] <<= 8 - l % 4 * 2;
    // do reverse complementary
    if (rev) {
        for (i = 0, j = b - 1; i < j; ++i, --j) {
            c = key[i];
            key[i] = comp_seq8_table[key[j]];
            key[j] = comp_seq8_table[c];
        }
        if (i == j) key[i] = comp_seq8_table[key[i]];
    }
    // do bit shift to align syncmer bytes
    for (i = 0; i < b-1; ++i) {
        key[i] <<= res;
        key[i] |= key[i+1] >> (8-res);
    }
    key[b-1] <<= res;
    // mask the lower bits
    key[b-1] &= lmask[w&3];

    for (i = 0; i < p1 - p0 + 1; ++i)
        // fputc(char_nt4_table[(key[i / 4] >> (3 - i % 4) * 2) & 3], fo);
        fputc(char_nt4_table[(key[i/4]>>(((i&3)^3)<<1)) & 3], fo);
    free(key);

    return;
}

void print_syncmer_on_seq(sr_t *sr, uint32_t n, int k, int w, FILE *fo)
{
    if (n >= sr->n) return;
    fprintf(fo, ">%lu_%d_%u_%lu_%u\t", sr->sid, n, sr->m_pos[n] >> 1, sr->s_mer[n] & 1, sr->m_pos[n] & 1);
    fprintf(fo, "RD:Z:%lu\t", sr->sid);
    fprintf(fo, "MM:Z:");
    fputs_smer(sr->s_mer[n] >> 1, k, fo);
    fputc('\t', fo);
    fprintf(fo, "KH:Z:%064lu%064lu\n", (uint64_t) (sr->k_mer_h[n] >> 64), (uint64_t) sr->k_mer_h[n]);
    fputs_kmer(sr->hoco_s, sr->m_pos[n], w, fo);
    fputc('\n', fo);
}

void print_all_syncmers_on_seq(sr_t *sr, int k, int w, FILE *fo)
{
    uint32_t i;
    for (i = 0; i < sr->n; ++i) print_syncmer_on_seq(sr, i, k, w, fo);
}

void print_seq(sr_t *sr, FILE *fo)
{
    uint32_t i;
    for (i = 0; i < sr->hoco_l; ++i) fputc(char_nt4_table[(sr->hoco_s[i/4]>>(((i&3)^3)<<1)) & 3], fo);
    fputc('\n', fo);
}

