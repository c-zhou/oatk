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
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

#include "misc.h"

double realtime0;

#if defined(WIN32) || defined(_WIN32)
#include <windows.h>

struct timezone
{
  __int32  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};

/*
 * gettimeofday.c
 *    Win32 gettimeofday() replacement
 *    taken from PostgreSQL, according to
 *    https://stackoverflow.com/questions/1676036/what-should-i-use-to-replace-gettimeofday-on-windows
 *
 * src/port/gettimeofday.c
 *
 * Copyright (c) 2003 SRA, Inc.
 * Copyright (c) 2003 SKC, Inc.
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose, without fee, and without a
 * written agreement is hereby granted, provided that the above
 * copyright notice and this paragraph and the following two
 * paragraphs appear in all copies.
 *
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE TO ANY PARTY FOR DIRECT,
 * INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
 * LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
 * DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE AUTHOR SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS
 * IS" BASIS, AND THE AUTHOR HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE,
 * SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 */

/* FILETIME of Jan 1 1970 00:00:00. */
static const unsigned __int64 epoch = ((unsigned __int64) 116444736000000000ULL);

/*
 * timezone information is stored outside the kernel so tzp isn't used anymore.
 *
 * Note: this function is not for Win32 high precision timing purpose. See
 * elapsed_time().
 */
static int gettimeofday(struct timeval * tp, struct timezone *tzp)
{
    FILETIME    file_time;
    SYSTEMTIME  system_time;
    ULARGE_INTEGER ularge;

    GetSystemTime(&system_time);
    SystemTimeToFileTime(&system_time, &file_time);
    ularge.LowPart = file_time.dwLowDateTime;
    ularge.HighPart = file_time.dwHighDateTime;

    tp->tv_sec = (long) ((ularge.QuadPart - epoch) / 10000000L);
    tp->tv_usec = (long) (system_time.wMilliseconds * 1000);

    return 0;
}

// taken from https://stackoverflow.com/questions/5272470/c-get-cpu-usage-on-linux-and-windows
double cputime(void)
{
    HANDLE hProcess = GetCurrentProcess();
    FILETIME ftCreation, ftExit, ftKernel, ftUser;
    SYSTEMTIME stKernel;
    SYSTEMTIME stUser;

    GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser);
    FileTimeToSystemTime(&ftKernel, &stKernel);
    FileTimeToSystemTime(&ftUser, &stUser);

    double kernelModeTime = ((stKernel.wHour * 60.) + stKernel.wMinute * 60.) + stKernel.wSecond * 1. + stKernel.wMilliseconds / 1000.;
    double userModeTime = ((stUser.wHour * 60.) + stUser.wMinute * 60.) + stUser.wSecond * 1. + stUser.wMilliseconds / 1000.;

    return kernelModeTime + userModeTime;
}

long peakrss(void) { return 0; }
#else
#include <sys/resource.h>
#include <sys/time.h>

double cputime(void)
{
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long peakrss(void)
{
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
    return r.ru_maxrss * 1024;
#else
    return r.ru_maxrss;
#endif
}

#endif /* WIN32 || _WIN32 */

double realtime(void)
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

void liftrlimit(void)
{
#ifdef __linux__
    struct rlimit r;
    getrlimit(RLIMIT_AS, &r);
    r.rlim_cur = r.rlim_max;
    setrlimit(RLIMIT_AS, &r);
#endif
}

void sys_init(void)
{
    liftrlimit();
    realtime0 = realtime();
}

#ifdef WIN32
#include <windows.h>
#elif _POSIX_C_SOURCE >= 199309L
#include <time.h>   // for nanosleep
#else
#include <unistd.h> // for usleep
#endif

void sleep_ms(int milliseconds) {
#ifdef WIN32
    Sleep(milliseconds);
#elif _POSIX_C_SOURCE >= 199309L
    struct timespec ts;
    ts.tv_sec = milliseconds / 1000;
    ts.tv_nsec = (milliseconds % 1000) * 1000000;
    nanosleep(&ts, NULL);
#else
    if (milliseconds >= 1000)
        sleep(milliseconds / 1000);
    usleep((milliseconds % 1000) * 1000);
#endif
}

char *make_tempfile(char *temp_dir, char *file_template, const char *suffix)
{
    sprintf(file_template, "%s/tmpXXXXXXXXXX%s", temp_dir, suffix);
    int ret = mkstemps(file_template, strlen(suffix));
    if (ret >= 0) {
        // need to close the file descriptor
        close(ret);
        return strdup(file_template);
    }
    return 0;
}

int run_system_cmd(char *cmd, int retry)
{
    int exit_code = system(cmd);
    --retry;
    if ((exit_code != -1 && !WEXITSTATUS(exit_code)) || !retry)
        return exit_code;
    return run_system_cmd(cmd, retry);
}

void check_executable(char *exe)
{
    char cmd[4096];
    sprintf(cmd, "command -v %s 1>/dev/null 2>/dev/null", exe);

    int exit_code = run_system_cmd(cmd, 1);
    if (exit_code == -1 || WEXITSTATUS(exit_code)) {
        fprintf(stderr, "[E::%s] executable %s is not available\n", __func__, exe);
        exit(EXIT_FAILURE);
    }
}

int is_file(const char *path) {
    struct stat path_stat;
    if (stat(path, &path_stat) == -1)
        return 0;  // stat error

    return S_ISREG(path_stat.st_mode);
}

int is_dir(const char *path) {
    struct stat path_stat;
    if (stat(path, &path_stat) == -1)
        return 0;  // stat error

    return S_ISDIR(path_stat.st_mode);
}

int is_fifo(const char *path) {
    if (*path == '-')
        return 1;

    struct stat path_stat;
    if (stat(path, &path_stat) == -1)
        return 0;  // stat error

    return S_ISFIFO(path_stat.st_mode);
}

FILE *open_outstream(char *prefix, char *suffix)
{
    FILE *fo;
    char *out;

    MYMALLOC(out, strlen(prefix) + strlen(suffix) + 1);
    sprintf(out, "%s%s", prefix, suffix);
    fo = fopen(out, "w");
    if (!fo) {
        fprintf(stderr, "[E::%s] failed to open file '%s' to write: %s\n", __func__, out, strerror(errno));
        exit(EXIT_FAILURE);
    }
    free(out);

    return fo;
}

void parse_pathname(char *path, char **_dirname, char **_basename)
{
    if (_dirname) *_dirname = 0;
    if (_basename) *_basename = 0;
    if (path == 0) return;
    int len, i, n_token;
    len = strlen(path);
    if (len == 0) return;
    if (path[len-1] == '/') {
        if (_dirname) *_dirname = strdup(path);
        return;
    }
    char *dirname, *basename, *token;
    dirname = strdup(path);
    token = strtok(dirname, "/");
    n_token = 0;
    do {
        ++n_token;
        basename = token;
        token = strtok(NULL, "/");
    } while (token);
    // necessary for duplicate '/'
    for (i = 0; i < len; ++i)
        if (dirname[i] == '/')
            dirname[i] = 0;
    len = basename - dirname - 1;
    for (i = 0; i < len; ++i)
        if (!dirname[i])
            dirname[i] = '/';
    if (dirname != basename) {
        if (len) {
            // here need to strdup basename
            if (_dirname) *_dirname = dirname;
            if (_basename) *_basename = strdup(basename);
        } else {
            // basename with point to dirname
            if (_dirname) *_dirname = strdup("/");
            if (_basename) *_basename = basename;
        }
    } else {
        if ((strlen(basename)==1 && basename[0] == '.') ||
                (strlen(basename)==2 && 
                 basename[0] == '.' && basename[1] == '.')) {
            if (_dirname) *_dirname = dirname;
        } else {
            if (_dirname) *_dirname = strdup(".");
            if (_basename) *_basename = basename;
        }
    }
    return;
}

