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
#include <stdlib.h>
#include <stdio.h>

#include "misc.h"
#include "sstream.h"

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

static kseq_stream_t *make_kseq_stream(char *file)
{
    kseq_stream_t *stream;
    MYCALLOC(stream, 1);
    int fd;
    // open query file
    stream->ko = kopen(file, &fd);
    if (stream->ko == 0) {
        fprintf(stderr, "[E::%s] fail to open file \"%s\"\n", __func__, file);
        exit(EXIT_FAILURE);
    }
    stream->fp = gzdopen(fd, "r");
    stream->ks = kseq_init(stream->fp);

    return stream;
}

static void kstream_close(kseq_stream_t *ks)
{
    kseq_destroy(ks->ks);
    gzclose(ks->fp);
    kclose(ks->ko);
    free(ks);
}

void sstream_close(sstream_t *ss)
{
    kstream_close(ss->s);
    free(ss);
}

sstream_t *sstream_open(char **files, int n_files)
{
    sstream_t *ss;
    MYCALLOC(ss, 1);
    ss->n_seq = 0;
    ss->files = files;
    ss->n_files = n_files;
    ss->n = 0;
    ss->s = make_kseq_stream(files[0]);

    return ss;
}

int sstream_read(sstream_t *ss)
{
    int l = kseq_read(ss->s->ks);
    if (l >= 0) {
        ++ss->n_seq;
        return l;
    }
    if (ss->n < ss->n_files - 1) {
        kstream_close(ss->s);
        ++ss->n;
        ss->s = make_kseq_stream(ss->files[ss->n]);

        l = kseq_read(ss->s->ks);
        if (l >= 0) {
            ++ss->n_seq;
            return l;
        }
    }
    return l;
}

