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
 * 20/01/23 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/

#ifndef __MISC_UTILITY_H__
#define __MISC_UTILITY_H__
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

extern double realtime0;

#define MIN(a, b) (((a)<(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))
#define SWAP(a, b) do { __typeof__(a) __tmp = (a); (a) = (b); (b) = __tmp; } while (0)

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif
#ifndef kroundup64
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif

#define MYMALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc(sizeof(*(ptr)) * (len)))
#define MYCALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MYREALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), sizeof(*(ptr)) * (len)))
#define MYBZERO(ptr, len) memset((ptr), 0, sizeof(*(ptr)) * (len))
#define MYBONE(ptr, len) memset((ptr), 0xff, sizeof(*(ptr)) * (len))
#define MYEXPAND(a, m) do { \
    ++(m); \
    kroundup64((m)); \
    MYREALLOC((a), (m)); \
} while (0)

#define array_reverse(a, n) do { \
    if ((n) > 0) { \
        size_t __i, __j; \
        __typeof__((a)[0]) __tmp; \
        for (__i = 0, __j = (n) - 1; __i < __j; ++__i, --__j) { \
            __tmp = (a)[__i]; (a)[__i] = (a)[__j]; (a)[__j] = __tmp; \
        } \
    } \
} while (0)

#ifdef __cplusplus
extern "C" {
#endif

double cputime(void);
long peakrss(void);
double realtime(void);
void liftrlimit(void);
void sys_init(void);
void sleep_ms(int ms);
void check_executable(char *exe);
int run_system_cmd(char *cmd, int retry);
int is_file(const char *path);
int is_dir(const char *path);
int is_fifo(const char *path);
char *make_tempfile(char *temp_dir, char *file_template, const char *suffix);
FILE *open_outstream(char *prefix, char *suffix);
void parse_pathname(char *path, char **_dirname, char **_basename);
#ifdef __cplusplus
}
#endif

#endif /*__MISC_UTILITY_H__*/



