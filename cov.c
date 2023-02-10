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
 * 02/09/22 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include "kvec.h"

#include "misc.h"
#include "cov.h"

cov_t *cov_init(uint32_t n)
{
    cov_t *cov;
    MYMALLOC(cov, 1);
    cov->n = n;
    MYMALLOC(cov->p, n);
    uint32_t i;
    for (i = 0; i < n; ++i)
        kv_init(cov->p[i]);

    return cov;
}

void cov_destroy(cov_t* cov)
{
    uint32_t i;
    for (i = 0; i < cov->n; ++i)
        kv_destroy(cov->p[i]);
    free(cov->p);
    free(cov);

    return;
}

static int u64_cmpfunc(const void *a, const void *b)
{
    return (*(uint64_t *) a > *(uint64_t *) b) - (*(uint64_t *) a < *(uint64_t *) b);
}

uint64_t cov_pos_compression(cov_t *cov)
{
    uint32_t i, p, p1;
    int32_t c;
    size_t j, k, n;
    uint64_t m, *a;
    
    m = 0;
    for (i = 0; i < cov->n; ++i) {
        a = cov->p[i].a;
        n = cov->p[i].n;
        if (n == 0) continue;

        qsort(a, n, sizeof(uint64_t), u64_cmpfunc);
        p = (uint32_t) (a[0] >> 32);
        c = (int32_t) a[0];
        k = 0;
        for(j = 1; j < n; ++j) {
            p1 = (uint32_t) (a[j] >> 32);
            if (p != p1) {
                if (c != 0) a[k++] = (uint64_t) p << 32 | (uint32_t) c;
                p = p1;
                c = (int32_t) a[j];
            } else {
                c += (int32_t) a[j];
            }
        }
        if (c != 0) a[k++] = (uint64_t) p << 32 | (uint32_t) c;
        cov->p[i].n = k;
        m += k;
    }

    return m;
}

void print_cov_in_bed_format(cov_t *cov, char **sname, uint32_t *slen, FILE *fo)
{
    size_t i, j, s;
    int32_t c;
    uint32_t p, p1;
    uint64_t *a;

    for (i = 0; i < cov->n; ++i) {
        a = cov->p[i].a;
        s = cov->p[i].n;

        p = 0, c = 0;
        for (j = 0; j < s; ++j) {
            p1 = (uint32_t) (a[j] >> 32);
            if (p1 > p)
                fprintf(fo, "%s\t%u\t%u\t%d\n", sname[i], p, p1, c);
            p = p1;
            c += (int32_t) a[j];
        }

        if (p < slen[i])
            fprintf(fo, "%s\t%u\t%u\t%d\n", sname[i], p, slen[i], 0);
    }
}

