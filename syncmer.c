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
#include <math.h>
#include <pthread.h>

#include "misc.h"
#include "khashl.h"
#include "kvec.h"
#include "syncmer.h"

#undef DEBUG_KMER_EXTRACTION
#undef DEBUG_S_KMER_GROUP
#undef DEBUG_LINK_COVERAGE
#undef DEBUG_CHECK_HASH_COLLISION

#define DO_HOCO_COMPRESSION

const unsigned char seq_nt4_table[256] = {
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

const char seq_nt4_comp_table[128] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
    32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
    48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', 123, 124, 125, 126, 127
};

const char char_nt4_table[4] = {'A', 'C', 'G', 'T'};

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

static uint8_t lmask[4] = {255, 192, 240, 252};
static uint64_t murmur3_seed = 1234;

uint64_t MurmurHash64A(const void *key, uint32_t len, uint64_t seed)
{
    const uint64_t m = 0xc6a4a7935bd1e995LLU;
    const int r = 47;

    uint64_t h = seed ^ (len * m);

    const uint64_t *data = (const uint64_t *)key;
    const uint64_t *end = (len >> 3) + data;

    while(data != end) {
        uint64_t k = *data++;

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    const unsigned char *data2 = (const unsigned char *) data;

    switch(len & 7) {
        case 7: h ^= (uint64_t) (data2[6]) << 48;
        case 6: h ^= (uint64_t) (data2[5]) << 40;
        case 5: h ^= (uint64_t) (data2[4]) << 32;
        case 4: h ^= (uint64_t) (data2[3]) << 24;
        case 3: h ^= (uint64_t) (data2[2]) << 16;
        case 2: h ^= (uint64_t) (data2[1]) << 8;
        case 1: h ^= (uint64_t) (data2[0]);
                h *= m;
    };

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

#ifdef DEBUG_KMER_EXTRACTION
static uint64_t kmer_hash64(uint64_t sid, uint8_t *s, uint32_t p, int w)
#else
static uint64_t kmer_hash64(uint8_t *s, uint32_t p, int w)
#endif
{
    uint64_t h64;
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

    h64 = MurmurHash64A(key, (w - 1) / 4 + 1, murmur3_seed);

#ifdef DEBUG_KMER_EXTRACTION
    fprintf(stderr, "[DEBUG_KMER_EXTRACTION::%s] sid:%lu p0:%d p1:%d rev:%d ", __func__, sid, p0, p1, rev);
    for (i = 0; i < p1 - p0 + 1; ++i)
        // fputc(char_nt4_table[(key[i / 4] >> (3 - i % 4) * 2) & 3], stderr);
        fputc(char_nt4_table[(key[i/4]>>(((i&3)^3)<<1)) & 3], stderr);

    fputc('\n', stderr);
#endif

    free(key);
    
    return h64;
}

typedef struct {
    int n_reads;
    uint64_t *sid;
    char **name;
    char **seq;
    int *len;
    sr_db_t *sr_db;
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
    k = dat->sr_db->s, w = dat->sr_db->k;
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
        kvec_t(uint64_t) k_mer;
    
        kv_init(hoco_s);
        kv_init(ho_rl);
        kv_init(ho_l_rl);
        kv_init(n_nucl);
        kv_init(m_pos);
        kv_init(s_mer);
        kv_init(k_mer);

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
                rl = 1;
#ifdef DO_HOCO_COMPRESSION    
                // hpc
                if (i + 1 < len && seq_nt4_table[(uint8_t)seq[i + 1]] == c) {
                    for (rl = 2; i + rl < len; ++rl)
                        if (seq_nt4_table[(uint8_t)seq[i + rl]] != c)
                            break;
                    i += rl - 1; // put $i at the end of the current homopolymer run
                }
#endif
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
                kv_push(uint64_t, k_mer, kmer_hash64(sr.sid, hoco_s.a, m_pos.a[m_pos.n-1], w));
#else
                kv_push(uint64_t, k_mer, kmer_hash64(hoco_s.a, m_pos.a[m_pos.n-1], w));
#endif
                // remove syncmers at the same position on a read
                // this is possible as a syncmer could start and end with the same smer
                if (m_pos.n >= 2 && m_pos.a[m_pos.n-1] >> 1 == m_pos.a[m_pos.n-2] >> 1) s_mer.n -= 2, m_pos.n -= 2, k_mer.n -= 2;
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
                    kv_push(uint64_t, k_mer, kmer_hash64(sr.sid, hoco_s.a, m_pos.a[m_pos.n-1], w));
#else
                    kv_push(uint64_t, k_mer, kmer_hash64(hoco_s.a, m_pos.a[m_pos.n-1], w));
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
                    kv_push(uint64_t, k_mer, kmer_hash64(sr.sid, hoco_s.a, m_pos.a[m_pos.n-1], w));
#else
                    kv_push(uint64_t, k_mer, kmer_hash64(hoco_s.a, m_pos.a[m_pos.n-1], w));
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
            kv_push(uint64_t, k_mer, kmer_hash64(sr.sid, hoco_s.a, m_pos.a[m_pos.n-1], w));
#else
            kv_push(uint64_t, k_mer, kmer_hash64(hoco_s.a, m_pos.a[m_pos.n-1], w));
#endif
            if (m_pos.n >= 2 && m_pos.a[m_pos.n-1] >> 1 == m_pos.a[m_pos.n-2] >> 1) s_mer.n -= 2, m_pos.n -= 2, k_mer.n -= 2;
        }
    
        if (hoco_s.n) MYREALLOC(hoco_s.a, hoco_s.n);
        if (ho_rl.n) MYREALLOC(ho_rl.a, ho_rl.n);
        if (ho_l_rl.n) MYREALLOC(ho_l_rl.a, ho_l_rl.n);
        if (n_nucl.n) MYREALLOC(n_nucl.a, n_nucl.n);
        if (m_pos.n) MYREALLOC(m_pos.a, m_pos.n);
        if (s_mer.n) MYREALLOC(s_mer.a, s_mer.n);
        if (k_mer.n) MYREALLOC(k_mer.a, k_mer.n);

        sr.hoco_l = hoco_l;
        sr.hoco_s = hoco_s.a;
        sr.ho_rl = ho_rl.a;
        sr.ho_l_rl = ho_l_rl.a;
        sr.n_nucl = n_nucl.a;
        sr.n = m_pos.n;
        sr.m_pos = m_pos.a;
        sr.s_mer = s_mer.a;
        sr.k_mer = k_mer.a;
        kv_push(sr_t, *dat->sr_db, sr);

        free(seq);
    }
    
    return NULL;
}

static void do_analysis(p_data_t *dat, pthread_t *threads, int n_threads, sr_db_t *sr_db)
{
    int t;
    for (t = 1; t < n_threads; ++t)
        pthread_create(threads + t, NULL, sr_read_analysis_thread, &dat[t]);

    sr_read_analysis_thread(dat);

    for (t = 1; t < n_threads; ++t)
        pthread_join(threads[t], NULL);

    for (t = 0; t < n_threads; ++t) {
        kv_pushn(sr_t, *sr_db, dat[t].sr_db->a, dat[t].sr_db->n);
        dat[t].n_reads = 0;
        dat[t].sr_db->n = 0;
    }

    return;
}

static void sr_read_single_thread(sstream_t *s_stream, sr_db_t *sr_db, size_t mD)
{
    sr_db_clean(sr_db);
    sr_db_init(sr_db, sr_db->k, sr_db->s);

    if (mD == 0)
        mD = SIZE_MAX;

    p_data_t *dat;
    MYCALLOC(dat, 1);
    dat->n_reads = 1;
    MYMALLOC(dat->sid, 1);
    MYMALLOC(dat->name, 1);
    MYMALLOC(dat->seq, 1);
    MYMALLOC(dat->len, 1);
    dat->sr_db = sr_db;

    int l;
    uint64_t i;
    size_t D;
    i = 0;
    D = 0;
    while ((l = sstream_read(s_stream)) >= 0) {
        dat->sid[0] = i++;
        dat->name[0] = strdup(s_stream->s->ks->name.s);
        dat->seq[0] = strdup(s_stream->s->ks->seq.s);
        dat->len[0] = l;
        sr_read_analysis_thread(dat);
        D += l;
        if (D >= mD) {
            fprintf(stderr, "[M::%s] data limit (%lu) reached. Discard the remaining sequences...\n", __func__, mD);
            break;
        }
    }

    free(dat->sid);
    free(dat->name);
    free(dat->seq);
    free(dat->len);
    free(dat);

    return;
}

void sr_read(sstream_t *s_stream, sr_db_t *sr_db, size_t mD, int n_threads)
{
    if (n_threads == 1) {
        sr_read_single_thread(s_stream, sr_db, mD);
        return;
    }

    sr_db_clean(sr_db);
    sr_db_init(sr_db, sr_db->k, sr_db->s);
    
    if (mD == 0)
        mD = SIZE_MAX;

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
        MYMALLOC(dat[t].sr_db, 1);
        sr_db_init(dat[t].sr_db, sr_db->k, sr_db->s);
        kv_resize(sr_t, *dat[t].sr_db, batch_n);
    }

    pthread_t threads[n_threads];
    int l;
    uint64_t i, j, n;
    size_t D;
    i = j = n = 0;
    D = 0;
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
            do_analysis(dat, threads, n_threads, sr_db);
            j = 0;
        }

        D += l;
        if (D >= mD) {
            fprintf(stderr, "[M::%s] data limit (%lu) reached. Discard the remaining sequences...\n", __func__, mD);
            break;
        }
    }
    if (j > 0) do_analysis(dat, threads, n_threads, sr_db);

    for (t = 0; t < n_threads; ++t) {
        free(dat[t].sid);
        free(dat[t].name);
        free(dat[t].seq);
        free(dat[t].len);
        free(dat[t].sr_db->a);
        free(dat[t].sr_db);
    }
    free(dat);

    return;
}

KHASHL_MAP_INIT(KH_LOCAL, kh_ctab_t, kh_ctab, khint_t, int, kh_hash_uint32, kh_eq_generic)

static int syncmer_s_cmpfunc(const void *a, const void *b)
{
    uint64_t x, y;
    uint64_t s, t;
    x = ((syncmer_t *) a)->s;
    y = ((syncmer_t *) b)->s;
    s = ((syncmer_t *) a)->h;
    t = ((syncmer_t *) b)->h;
    return x == y? ((s > t) - (s < t)) : ((x > y) - (x < y));
}

static int syncmer_h_cmpfunc(const void *a, const void *b)
{
    uint64_t x, y;
    uint64_t s, t;
    x = ((syncmer_t *) a)->h;
    y = ((syncmer_t *) b)->h;
    s = ((syncmer_t *) a)->s;
    t = ((syncmer_t *) b)->s;
    return x == y? ((s > t) - (s < t)) : ((x > y) - (x < y));
}

static int int64_cmpfunc(const void *a, const void *b)
{
    int64_t x, y;
    x = *(int64_t *) a;
    y = *(int64_t *) b;
    return (x > y) - (x < y);
}

static int uint128_cmpfunc(const void *a, const void *b)
{
    uint128_t x, y;
    x = *(uint128_t *) a;
    y = *(uint128_t *) b;
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
    if (k < kh_end(ctab))
        c = kh_val(ctab, k);

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
    if (verbose > 0) {
        if (max2 > 0)
            fprintf(stderr, "[M::%s] left: count[%d] = %ld\n", __func__, max2_i, (long)cnt[max2_i]);
        else
            fprintf(stderr, "[M::%s] left: none\n", __func__);
    }
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
    if (verbose > 0) {
        if (max3 > 0)
            fprintf(stderr, "[M::%s] right: count[%d] = %ld\n", __func__, max3_i, (long)cnt[max3_i]);
        else
            fprintf(stderr, "[M::%s] right: none\n", __func__);
    }
    if (max3_i > 0) {
        *peak_het = max_i;
        return max3_i;
    } else {
        if (max2_i > 0) *peak_het = max2_i;
        return max_i;
    }
}

void sr_db_stat(sr_db_t *sr_db, FILE *fo, int verbose)
{
    size_t i, j, n, m;
    int c, w, p0, p1;
    uint64_t s;
    uint64_t h;
    kvec_t(syncmer_t) syncmers;
    kh_ctab_t *dist_ctab, *smer_ctab, *kmer_ctab, *smer_k_ctab;

    sr_stat_t *stats;
    stats = sr_db->stats;
    if (!stats) {
        MYCALLOC(stats, 1);
        sr_db->stats = stats;
    }

    dist_ctab = kh_ctab_init();
    smer_ctab = kh_ctab_init();
    kmer_ctab = kh_ctab_init();
    smer_k_ctab = kh_ctab_init();
    kv_init(syncmers);
    w = sr_db->k;
    n = sr_db->n;
    m = 0;
    for (i = 0; i < n; ++i) {
        sr_t s = sr_db->a[i];
        m += s.n;
        p0 = p1 = MAX_RD_LEN;
        for (j = 0; j < s.n; ++j) {
            syncmer_t a = {s.k_mer[j] >> 1, s.s_mer[j], 0, 0, 0};
            kv_push(syncmer_t, syncmers, a);
            p0 = p1;
            p1 = s.m_pos[j] >> 1;
            if (p0 != (int) MAX_RD_LEN && p1 != (int) MAX_RD_LEN)
                kh_ctab_put1(dist_ctab, p1 - p0 - w);
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

int sr_db_validate(sr_db_t *sr_db)
{
    if (sr_db->n > MAX_RD_NUM) {
        fprintf(stderr, "[E::%s] read number exceeds the limit %llu\n", __func__, MAX_RD_NUM);
        return 1;
    }
    size_t i;
    for (i = 0; i < sr_db->n; ++i) {
        sr_t *s = &sr_db->a[i];
        if (s->n > MAX_RD_SCM) {
            fprintf(stderr, "[E::%s] syncmer number (%u) on read exceeds the limit %llu: %s\n", __func__, s->n, MAX_RD_SCM, s->sname);
            return 2;
        }
    }
    return 0;
}

void sr_destroy(sr_t *sr)
{
    if (!sr) return;
    if (sr->sname) free(sr->sname);
    if (sr->hoco_s) free(sr->hoco_s);
    if (sr->ho_rl) free(sr->ho_rl);
    if (sr->ho_l_rl) free(sr->ho_l_rl);
    if (sr->n_nucl) free(sr->n_nucl);
    if (sr->s_mer) free(sr->s_mer);
    if (sr->k_mer) free(sr->k_mer);
    if (sr->m_pos) free(sr->m_pos);
}

void sr_db_init(sr_db_t *sr_db, int k, int s)
{
    if (!sr_db) return;
    kv_init(*sr_db);
    sr_db->k = k;
    sr_db->s = s;
    sr_db->stats = 0;
}

void sr_db_clean(sr_db_t *sr_db)
{
    if (!sr_db) return;
    size_t i;
    for (i = 0; i < sr_db->n; ++i)
        sr_destroy(&sr_db->a[i]);
    kv_destroy(*sr_db);
    free(sr_db->stats);
}

void sr_db_destroy(sr_db_t *sr_db)
{
    if (!sr_db) return;
    sr_db_clean(sr_db);
    free(sr_db);
}

void syncmer_db_init(syncmer_db_t *scm_db)
{
    if (!scm_db) return;
    kv_init(*scm_db);
    scm_db->c = 0;
    scm_db->h = 0;
}

void syncmer_db_clean(syncmer_db_t *scm_db)
{
    if (!scm_db) return;
    size_t i;
    for (i = 0; i < scm_db->n; ++i)
        free(scm_db->a[i].m_pos);
    kv_destroy(*scm_db);
    free(scm_db->c);
    free(scm_db->h);
}

void syncmer_db_destroy(syncmer_db_t *scm_db)
{
    if (!scm_db) return;
    syncmer_db_clean(scm_db);
    free(scm_db);
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
    if ((sr->m_pos[n]>>1) == MAX_RD_LEN) return; // this is a corrected mer
    fprintf(fo, ">%lu_%d_%u_%lu_%u\t", sr->sid, n, sr->m_pos[n] >> 1, sr->s_mer[n] & 1, sr->m_pos[n] & 1);
    fprintf(fo, "RD:Z:%lu\t", sr->sid);
    fprintf(fo, "MM:Z:");
    fputs_smer(sr->s_mer[n] >> 1, k, fo);
    fputc('\t', fo);
    fprintf(fo, "KH:Z:%lu\n", sr->k_mer[n]);
    fputs_kmer(sr->hoco_s, sr->m_pos[n], w, fo);
    fputc('\n', fo);
}

void print_all_syncmers_on_seq(sr_t *sr, int k, int w, FILE *fo)
{
    uint32_t i;
    for (i = 0; i < sr->n; ++i) print_syncmer_on_seq(sr, i, k, w, fo);
}

void print_aligned_syncmers_on_seq(sr_t *sr, int w, uint32_t beg, uint32_t end, FILE *fo)
{
    uint32_t i, j, l, n, p;
    kstring_t s = {0, 0, 0};
    l = sr->hoco_l;
    n = sr->n;
    end = MIN(end, n);
    get_hoco_seq(sr, &s);
    fprintf(fo, "%.*s\n", (int) s.l, s.s);
    for (i = beg; i < end; ++i) {
        p = sr->m_pos[i]>>1;
        if (p == MAX_RD_LEN) continue;
        for (j = 0; j < p; ++j) fputc('*', fo);
        fprintf(fo, "%.*s", w, &s.s[p]);
        for (j = p + w; j < l; ++j) fputc('*', fo);
        fputc('\n', fo);
    }
    free(s.s);
}

void print_hoco_seq(sr_t *sr, FILE *fo)
{
    uint32_t i;
    for (i = 0; i < sr->hoco_l; ++i) fputc(char_nt4_table[(sr->hoco_s[i/4]>>(((i&3)^3)<<1)) & 3], fo);
    fputc('\n', fo);
}

void get_hoco_seq(sr_t *sr, kstring_t *s)
{
    uint32_t i;
    s->l = 0;
    for (i = 0; i < sr->hoco_l; ++i) kputc_(char_nt4_table[(sr->hoco_s[i/4]>>(((i&3)^3)<<1)) & 3], s);
}

void get_kmer_seq(uint8_t *hoco_s, uint32_t pos, int l, uint32_t rev, uint8_t *kmer_s)
{
    int i, j;
    uint8_t t;
    uint32_t p;
    for (i = 0; i < l; ++i) {
        p = pos + i;
        kmer_s[i] = (hoco_s[p/4]>>(((p&3)^3)<<1))&3;
    }
    if (rev) {
        for (i = 0, j = l - 1; i < j; ++i, --j) {
            t = kmer_s[i];
            kmer_s[i] = (kmer_s[j]^3)&3;
            kmer_s[j] = (t^3)&3;
        }
        if (i == j) kmer_s[i] = (kmer_s[i]^3)&3;
    }
}

void get_kmer_dna_seq(uint8_t *hoco_s, uint32_t pos, int l, uint32_t rev, char *dna_seq)
{
    int i, j;
    uint32_t p;
    for (i = 0; i < l; ++i) {
        p = pos + i;
        dna_seq[i] = char_nt4_table[(hoco_s[p/4]>>(((p&3)^3)<<1))&3];
    }
    if (rev) {
        char t;
        for (i = 0, j = l - 1; i < j; ++i, --j) {
            t = dna_seq[i];
            dna_seq[i] = seq_nt4_comp_table[(int) dna_seq[j]];
            dna_seq[j] = seq_nt4_comp_table[(int) t];
        }
        if (i == j) dna_seq[i] = seq_nt4_comp_table[(int) dna_seq[j]];
    }
}

static inline int equal_array(void *a, void *b, int n)
{
    uint64_t *x = (uint64_t *) a;
    uint64_t *y = (uint64_t *) b;
    int i;
    for (i = 0; i < n; ++i)
        if (x[i] != y[i])
            return 0;
    return 1;
}

// process a cluster of kmers with the same hash value
// check hash collisions
// add kmer to database
static void process_kmer_cluster(uint128_t *scm, uint32_t n, syncmer_db_t *scm_db, sr_db_t *sr_db)
{
    int n_clus, *clus;
    MYMALLOC(clus, n);

    if (n == 1) {
        // no hash collision for sure
        n_clus = 1;
        clus[0] = 0;
    } else {
        int i, j, k, b, c, B, C, p0, p1, rev, res;
        uint8_t *kmer, *kmer_list;
        uint32_t n_s;
        uint64_t sid;
        size_t s;

        k = sr_db->k;
        B = (((k - 1) >> 3) + 1) << 3; // maxmum number bytes to hold syncmer 64bit aligned
        C = B >> 3; // number of comparsions of 64bit interger
        n_clus = 0;
        kmer_list = 0;
        MYMALLOC(kmer, B);

        for (s = 0; s < n; ++s) {
            sid = (uint64_t) scm[s] >> 32;
            n_s = (uint32_t) ((uint64_t) scm[s]) >> 1;
            rev = scm[s] & 1;
            p0 = sr_db->a[sid].m_pos[n_s] >> 1;
            p1 = p0 + k - 1;
            res = rev? ((p1&3)^3)<<1 : (p0&3)<<1;
            b = p1 / 4 - p0 / 4 + 1; // actual number bytes holding syncmer
            memcpy(kmer, sr_db->a[sid].hoco_s + p0 / 4, b);
            // do reverse complementary
            if (rev) {
                for (i = 0, j = b - 1; i < j; ++i, --j) {
                    c = kmer[i];
                    kmer[i] = comp_seq8_table[kmer[j]];
                    kmer[j] = comp_seq8_table[c];
                }
                if (i == j) kmer[i] = comp_seq8_table[kmer[i]];
            }
            // do bit shift to align syncmer bytes
            for (i = 0; i < b - 1; ++i) {
                kmer[i] <<= res;
                kmer[i] |= kmer[i+1] >> (8-res);
            }
            kmer[b-1] <<= res;
            // mask the lower bits
            kmer[b-1] &= lmask[k&3];
            // mask the remaining bytes
            MYBZERO(&kmer[b], B - b);
            // compare to the existing kmers to find hash collision
            for (i = 0; i < n_clus; ++i)
                if (equal_array(kmer, &kmer_list[i * B], C))
                    break;
            clus[s] = i;
            if (i >= n_clus) {
                // new kmer
                // add to kmer list
                ++n_clus;
                MYREALLOC(kmer_list, n_clus * B);
                memcpy(&kmer_list[(n_clus - 1) * B], kmer, B);
            }
        }
#ifdef DEBUG_CHECK_HASH_COLLISION
        uint64_t h64 = (uint64_t) (scm[0] >> 64);
        for (i = 0; i < n_clus; ++i)
            assert(h64 == MurmurHash64A(&kmer_list[i * B], (k - 1) / 4 + 1, murmur3_seed));
#endif
        free(kmer);
        free(kmer_list);
    }

    // add each cluster to syncmer database
    uint32_t s, *cnts;
    MYCALLOC(cnts, n_clus);
    for (s = 0; s < n; ++s)
        ++cnts[clus[s]];

    syncmer_t *syncmer;
    int i;
    for (i = 0; i < n_clus; ++i) {
        kv_pushp(syncmer_t, *scm_db, &syncmer);
        syncmer->h = (uint64_t) (scm[0] >> 64);
        syncmer->s = UINT64_MAX; // UINT64_MAX cannot be a smer as the first bit of smer is always zero
        syncmer->cov = 0;
        syncmer->del = 0;
        MYMALLOC(syncmer->m_pos, cnts[i]);
    }

    uint64_t smer;
    for (s = 0; s < n; ++s) {
        // update syncmer database
        syncmer = &scm_db->a[scm_db->n - n_clus + clus[s]];
        syncmer->m_pos[syncmer->cov++] = (uint64_t) scm[s];
        smer = sr_db->a[(uint64_t) scm[s] >> 32].s_mer[(uint32_t) ((uint64_t) scm[s]) >> 1];
        if (syncmer->s == UINT64_MAX) {
            syncmer->s = smer;
        } else if (syncmer->s != smer) {
            fprintf(stderr, "[E::%s] identical kmers have different smers\n", __func__);
            fprintf(stderr, "[E::%s] kmer hash  : %lu\n", __func__, syncmer->h);
            fprintf(stderr, "[E::%s] smer code 0: %lu; read id: %lu\n", __func__, syncmer->s, syncmer->m_pos[0] >> 32);
            fprintf(stderr, "[E::%s] smer code 1: %lu; read id: %lu\n", __func__, smer, (uint64_t) scm[s] >> 32);
            exit(EXIT_FAILURE);
        }
        // update read database
        sr_db->a[(uint64_t) scm[s] >> 32].k_mer[(uint32_t) ((uint64_t) scm[s]) >> 1] = (scm_db->n - n_clus + clus[s]) << 1;
    }

#ifdef DEBUG_CHECK_HASH_COLLISION
    if (n_clus > 1) {
        fprintf(stderr, "[DEBUG_CHECK_HASH_COLLISION::%s] hash collision: hash = %lu; clus = %d\n", 
                __func__, (uint64_t) (scm[0] >> 64), n_clus);
        for (i = 0; i < n_clus; ++i)
            fprintf(stderr, "[DEBUG_CHECK_HASH_COLLISION::%s] clus %d: read id %lu\n", 
                    __func__, i, scm_db->a[scm_db->n - n_clus + i].m_pos[0] >> 32);
    }
#endif    

    free(cnts);
    free(clus);
}

// make syncmer database from reads
// change the read kmer hash to kmer id
syncmer_db_t *collect_syncmer_from_reads(sr_db_t *sr_db)
{
    size_t i, j;
    uint64_t n1, n2;
    sr_t *s;
    kvec_t(uint128_t) scm;
    kv_init(scm);
    n1 = 0;
    for (i = 0; i < sr_db->n; ++i) {
        s = &sr_db->a[i];
        assert(s->sid == i);
        n1 += s->n;
        for (j = 0; j < s->n; ++j) {
            uint128_t s1 = (uint128_t) s->k_mer[j] << 64 | ((s->sid << 32) | (j << 1) | (s->m_pos[j] & 1));
            kv_push(uint128_t, scm, s1);
        }
    }
    if (scm.n == 0) {
        kv_destroy(scm);
        return 0;
    }

    qsort(scm.a, scm.n, sizeof(uint128_t), uint128_cmpfunc);

    // pack syncmers by kmer hash
    // check hash collisions
    syncmer_db_t *scm_db;
    size_t last;
    uint64_t h64;
    MYMALLOC(scm_db, 1);
    syncmer_db_init(scm_db);
    h64 = (uint64_t) (scm.a[0] >> 64);
    for (i = 1, last = 0; i < scm.n; ++i) {
        if ((uint64_t) (scm.a[i] >> 64) != h64) {
            // process kmer cluster
            process_kmer_cluster(&scm.a[last], i - last, scm_db, sr_db);

            last = i;
            h64 = (uint64_t) (scm.a[i] >> 64);
        }
    }
    process_kmer_cluster(&scm.a[last], i - last, scm_db, sr_db);
    kv_destroy(scm);

    MYREALLOC(scm_db->a, scm_db->n);
    scm_db->m = scm_db->n;
    MYMALLOC(scm_db->c, scm_db->n);
    for (i = 0; i < scm_db->n; ++i) scm_db->c[i] = 1;

    n2 = 0;
    for (i = 0; i < scm_db->n; ++i) n2 += scm_db->a[i].cov;
    assert(n1 == n2);

    return scm_db;
}

static kh_inline khint_t kh_hash_uint128(uint128_t key)
{
    khint_t k1 = kh_hash_uint64((khint64_t) key);
    khint_t k2 = kh_hash_uint64((khint64_t) (key >> 64));
    return kh_hash_uint64((khint64_t) ((uint64_t) k1 << 32 | k2));
}

KHASHL_MAP_INIT(KH_LOCAL, kh_scm_t, kh_scm, uint128_t, uint64_t, kh_hash_uint128, kh_eq_generic)
KHASHL_MAP_INIT(KH_LOCAL, kh_kcn_t, kh_kcn, uint64_t, int, kh_hash_uint64, kh_eq_generic)

static inline void kh_scm_add(kh_scm_t *h, uint128_t v_p, uint64_t c)
{
    int absent;
    khint_t k;
    k = kh_scm_put(h, v_p, &absent);
    if (absent)
        kh_val(h, k) = 0;
    kh_val(h, k) += c;
    return;
}

static inline void kh_kcn_put1(kh_kcn_t *k_cn, uint64_t key)
{
    int absent;
    khint_t k;
    k = kh_kcn_put(k_cn, key, &absent);
    if (absent)
        kh_val(k_cn, k) = 1;
    else
        ++kh_val(k_cn, k);
}

static inline int kh_kcn_get1(kh_kcn_t *k_cn, uint64_t key)
{
    khint_t k;
    k = kh_kcn_get(k_cn, key);
    if (k < kh_end(k_cn))
        return kh_val(k_cn, k);
    else
        return 0;
}

typedef struct {uint32_t c, l; double f;} pt1_t;

static int pt1_f_cmpfunc(const void *a, const void *b)
{
    if (((pt1_t *) a)->f == ((pt1_t *) b)->f)
        return (((pt1_t *) a)->c > ((pt1_t *) b)->c) - (((pt1_t *) a)->c < ((pt1_t *) b)->c);
    return (((pt1_t *) a)->f > ((pt1_t *) b)->f) - (((pt1_t *) a)->f < ((pt1_t *) b)->f);
}

/***
static int pt1_c_cmpfunc(const void *a, const void *b)
{
    if (((pt1_t *) a)->c == ((pt1_t *) b)->c)
        return (((pt1_t *) a)->l > ((pt1_t *) b)->l) - (((pt1_t *) a)->l < ((pt1_t *) b)->l);
    return (((pt1_t *) a)->c > ((pt1_t *) b)->c) - (((pt1_t *) a)->c < ((pt1_t *) b)->c);
}

static int pt1_l_cmpfunc(const void *a, const void *b)
{
    if (((pt1_t *) a)->l == ((pt1_t *) b)->l)
        return (((pt1_t *) a)->c > ((pt1_t *) b)->c) - (((pt1_t *) a)->c < ((pt1_t *) b)->c);
    return (((pt1_t *) a)->l > ((pt1_t *) b)->l) - (((pt1_t *) a)->l < ((pt1_t *) b)->l);
}
**/

// this subroutine is to analyse the relationship between 
// read coverage and link counts given the distance (measured by syncmer distance)
// tries to construct a series of linear regressions 
// N_LINK = beta * N_COV_d
// for disntace d = 2...D (upper bound D subject to data amount)
int syncmer_link_coverage_analysis(sr_db_t *sr_db, syncmer_db_t *scm_db, uint32_t min_k_cov, uint32_t min_n_seq, 
        uint32_t min_pt, double min_f, double **_beta, double **_bse, double **_r2, int verbose)
{
    uint64_t i, j;
    syncmer_t *scm;
    
    min_pt = MAX(min_pt, 30);
    min_f  = MAX(min_f,  .0);

    scm = scm_db->a;

    // collect read length - syncmer count
    kh_ctab_t *rl_ctab;
    int64_t *rl_cnts, *rd_cnts;
    uint32_t max_n;
    rl_ctab = kh_ctab_init();
    max_n = 0;
    for (i = 0; i < sr_db->n; ++i) {
        max_n = MAX(max_n, sr_db->a[i].n);
        kh_ctab_put1(rl_ctab, sr_db->a[i].n);
    }
    if (max_n == 0) return 0;

    rl_cnts = kh_ctab_cnt(rl_ctab, max_n);
    kh_ctab_destroy(rl_ctab);
    
    // rl_cnts cumsum from end to beg
    // rl_cnts size is max_n+1
    for (i = 0; i < max_n; ++i)
        rl_cnts[max_n-i-1] += rl_cnts[max_n-i];

#ifdef DEBUG_LINK_COVERAGE
    for (i = 0; i <= max_n; ++i)
        fprintf(stderr, "[DEBUG_LINK_COVERAGE::%s] rl_cnts: %lu %ld\n", __func__, i, rl_cnts[i]);
#endif

    kh_scm_t *a_cov; // arc cov
    kh_kcn_t *k_cn; // scm copy number
    sr_t *s;
    uint32_t a, c, v_v, n1, *pt_n;
    uint64_t v0, v1, *scm_id;
    uint128_t v_p;
    double *beta, *bse, *r2, c0, c1;
    kvec_t(pt1_t) pt1;
    khint_t k;

    a_cov = kh_scm_init();
    k_cn = kh_kcn_init();
    kv_init(pt1);
    MYMALLOC(scm_id, max_n);
    MYCALLOC(rd_cnts, max_n+1);
    MYCALLOC(beta, max_n);
    MYCALLOC(bse, max_n);
    MYCALLOC(r2, max_n);
    MYCALLOC(pt_n, max_n);
    n1 = 0;
    for (i = 2; i < max_n; ++i) {
        if (rl_cnts[i] < min_n_seq) break;
        // do coverage analysis
        // the gap will be i-2
        kh_scm_m_clear(a_cov);
        for (j = 0; j < sr_db->n; ++j) {
            s = &sr_db->a[j];
            if (s->n < i) continue;
            // collect scm_id along reads
            for (a = 0; a < s->n; ++a)
                scm_id[a] = s->k_mer[a] >> 1;
            // collect coverage
            for (a = i-1; a < s->n; ++a) {
                if (scm[scm_id[a+1-i]].cov < min_k_cov || 
                        scm[scm_id[a]].cov < min_k_cov)
                    continue;
                v0 = scm_id[a+1-i] << 1 | (s->m_pos[a+1-i] & 1);
                v1 = scm_id[a] << 1 | (s->m_pos[a] & 1);
                v0 <= v1? kh_scm_add(a_cov, (uint128_t) v0 << 64 | v1, 1) :
                    kh_scm_add(a_cov, ((uint128_t) v1^1) << 64 | (v0^1), 1);
                // rd_cnts collects the pair number of a given gap in the data
                ++rd_cnts[i];
            }
        }

        if (i == 2) {
            // kmer copy number estimation
            for (k = (khint_t) 0; k < kh_end(a_cov); ++k) {
                if (!kh_exist(a_cov, k)) continue;
                v_p = kh_key(a_cov, k);
                kh_kcn_put1(k_cn, (uint64_t) (v_p >> 65));
                kh_kcn_put1(k_cn, ((uint64_t) v_p) >> 1);
            }

#ifdef DEBUG_LINK_COVERAGE
            for (k = (khint_t) 0; k < kh_end(k_cn); ++k) {
                if (!kh_exist(k_cn, k)) continue;
                fprintf(stderr, "[DEBUG_LINK_COVERAGE::%s] k_cn: %lu %d\n", __func__, i-2, kh_val(k_cn, k));
            }
#endif
        }

#ifdef DEBUG_LINK_COVERAGE
        if ((i-2) % 3 == 0) {
            // print arc cov stats
            for (k = (khint_t) 0; k < kh_end(a_cov); ++k) {
                if (!kh_exist(a_cov, k)) continue;
                v_p = kh_key(a_cov, k);
                v_v = kh_val(a_cov, k);
                v0 = (uint64_t) (v_p >> 65);
                v1 = ((uint64_t) v_p) >> 1;
                c0 = MAX(2, kh_kcn_get1(k_cn, v0)) / 2.0;
                c1 = MAX(2, kh_kcn_get1(k_cn, v1)) / 2.0;
                c = (uint32_t) (MIN(scm[v0].cov/c0, scm[v1].cov/c1));
                v_v = MIN(v_v, c);
                fprintf(stderr, "[DEBUG_LINK_COVERAGE::%s] lc_cnts: %lu %u %u\n", __func__, i-2, c, v_v);
            }
        }
#endif

        // collect covs
        pt1.n = 0;
        for (k = (khint_t) 0; k < kh_end(a_cov); ++k) {
            if (!kh_exist(a_cov, k)) continue;
            v_p = kh_key(a_cov, k);
            v_v = kh_val(a_cov, k);
            v0 = (uint64_t) (v_p >> 65);
            v1 = ((uint64_t) v_p) >> 1;
            c0 = MAX(2, kh_kcn_get1(k_cn, v0)) / 2.0;
            c1 = MAX(2, kh_kcn_get1(k_cn, v1)) / 2.0;
            c = (uint32_t) (MIN(scm[v0].cov/c0, scm[v1].cov/c1));
            v_v = MIN(v_v, c);
            pt1_t pt = {c, v_v, (double) v_v / c};
            kv_push(pt1_t, pt1, pt);
        }
        
        uint32_t beg, end;
        /***
        beg = floor(pt1.n * .01);
        end =  ceil(pt1.n * .99);
        // sort pt by c
        qsort(pt1.a, pt1.n, sizeof(pt1_t), pt1_c_cmpfunc);
        // mask the lower 1% and upper 1% of data
        for (j = 0; j < beg; ++j) pt1.a[j].f = -1.0;
        for (j = end; j < pt1.n; ++j) pt1.a[j].f = -1.0;
        // sort pt by l
        qsort(pt1.a, pt1.n, sizeof(pt1_t), pt1_l_cmpfunc);
        // mask the lower 1% and upper 1% of data
        for (j = 0; j < beg; ++j) pt1.a[j].f = -1.0;
        for (j = end; j < pt1.n; ++j) pt1.a[j].f = -1.0;
        **/
        beg = floor(pt1.n * .05);
        end =  ceil(pt1.n * .95);
        // sort pt by l/c
        qsort(pt1.a, pt1.n, sizeof(pt1_t), pt1_f_cmpfunc);
        // estimate slope using the middle 90% data
        while (beg < end && pt1.a[beg].f < min_f)
            ++beg;
        
        if (end - beg < min_pt)
            break;

        if (verbose > 1) {
            for (j = beg; j < end; ++j)
                fprintf(stderr, "[M::%s] pt1: %lu %lu %u %u %.6f\n", __func__, i-2, j, pt1.a[j].c, pt1.a[j].l, pt1.a[j].f);
        }

        // do estimation
        double xy, x2, ybar;
        xy = x2 = ybar = .0;
        for (j = beg; j < end; ++j) {
            xy += (double) pt1.a[j].c * pt1.a[j].l;
            x2 += (double) pt1.a[j].c * pt1.a[j].c;
            ybar += pt1.a[j].l;
        }
        beta[i] = xy / x2;
        ybar /= end - beg;
        // calculate r2
        double res, tot, r;
        res = tot = .0;
        for (j = beg; j < end; ++j) {
            r = pt1.a[j].l - beta[i] * pt1.a[j].c;
            res += r * r;
            r = pt1.a[j].l - ybar;
            tot += r * r;
        }
        bse[i] = sqrt(res / x2 / (end - beg - 1));
        r2[i] = 1 - (tot == .0? .0 : res / tot);
        pt_n[i] = end - beg;
        n1 = i;
    }

#ifdef DEBUG_LINK_COVERAGE
    for (i = 0; i <= max_n; ++i)
        fprintf(stderr, "[DEBUG_LINK_COVERAGE::%s] rd_cnts: %lu %ld\n", __func__, i, rd_cnts[i]);
    for (i = 2; i < max_n; ++i)
        fprintf(stderr, "[DEBUG_LINK_COVERAGE::%s] rd_regs: %lu %u %.6f %.6f %.6f\n", __func__, 
                i-2, pt_n[i], beta[i], bse[i], r2[i]);
#endif

    if (verbose > 0) {
        for (i = 2; i < n1; ++i)
            fprintf(stderr, "[M::%s] G: %lu N: %u D: %ld coeff: %.6f bse: %.6f R2: %.6f\n", __func__, 
                    i-2, pt_n[i], rd_cnts[i], beta[i], bse[i], r2[i]);
    }

    if (n1 > 0) {
        // copy results
        if (_beta) {
            MYMALLOC(*_beta, n1 - 1);
            memcpy(*_beta, &beta[2], n1 - 1);
        }
        if (_bse) {
            MYMALLOC(*_bse, n1 - 1);
            memcpy(*_bse, &bse[2], n1 - 1);
        }
        if (_r2) {
            MYMALLOC(*_r2, n1 - 1);
            memcpy(*_r2, &r2[2], n1 - 1);
        }
    }

    kh_scm_destroy(a_cov);
    kh_kcn_destroy(k_cn);
    kv_destroy(pt1);
    free(scm_id);
    free(rl_cnts);
    free(rd_cnts);
    free(beta);
    free(bse);
    free(r2);
    free(pt_n);

    return n1 > 0? (n1 - 1) : 0;
}

