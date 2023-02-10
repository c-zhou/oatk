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
 * 03/08/22 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "ketopt.h"
#include "kvec.h"
#include "kseq.h"
#include "path.h"

#define PATHTOFASTA_VERSION "0.1"

KSTREAM_INIT(gzFile, gzread, 65536)

int VERBOSE = 0;
double realtime0;

static ko_longopt_t long_options[] = {
    { "linear",  ko_no_argument, 300 },
    { "version", ko_no_argument, 'V' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    const char *opt_str = "p:s:o:Vv:h";
    ketopt_t opt = KETOPT_INIT;
    int c, force_linear, ret = 0;
    asg_t *g;
    FILE *fp_help, *out_seqs;
    char *out, *seq_id, *path_file, *path_str;
    
    sys_init();

    fp_help = stderr;
    out_seqs = stdout;
    force_linear = 0;
    out = seq_id = path_file = path_str = 0;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >=0 ) {
        if (c == 's') seq_id = opt.arg;
        else if (c == 'p') path_file = opt.arg;
        else if (c == 300) force_linear = 1;
        else if (c == 'v') VERBOSE = atoi(opt.arg);
        else if (c == 'h') fp_help = stdout;
        else if (c == 'o') {
            if (strcmp(opt.arg, "-") != 0) {
                out = opt.arg;
            }
        }
        else if (c == 'V') {
            puts(PATHTOFASTA_VERSION);
            return 0;
        }
        else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
        else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }
    
    if (argc == opt.ind || fp_help == stdout) {
        fprintf(fp_help, "Usage: path_to_fasta [options] file[.gfa[.gz]] [path_str]\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "    -p STR        two-column path file\n");
        fprintf(fp_help, "    -s STR        output sequence id\n");
        fprintf(fp_help, "    -o FILE       output results to FILE [stdout]\n");
        fprintf(fp_help, "    -v INT        verbose level [%d]\n", VERBOSE);
        fprintf(fp_help, "    --linear      force linear output\n");
        fprintf(fp_help, "    --version     show version number\n");
        return fp_help == stdout? 0 : 1;
    }

    if (argc - opt.ind < 1) {
        fprintf(stderr, "[E::%s] missing input: please specify the GFA file\n", __func__);
        return 1;
    }

    if (argc - opt.ind > 1) {
        path_str = argv[opt.ind + 1];
        if (path_file) {
            fprintf(stderr, "[W::%s] path file %s ignored\n", __func__, path_file);
            fprintf(stderr, "[W::%s] using path string: %s\n", __func__, path_str);
        }
    } else {
        if (!path_file) {
            fprintf(stderr, "[E::%s] need a path file (-p) or path string\n", __func__);
            return 1;
        }
    }

    if (out) {
        out_seqs = fopen(out, "w");
        if (!out_seqs) {
            fprintf(stderr, "[E::%s] failed to open file %s to write\n", __func__, out);
            return 1;
        }
    }

    g = asg_read(argv[opt.ind]);
    if (g == 0) {
        fprintf(stderr, "[E::%s] failed to read the graph: %s\n", __func__, argv[opt.ind]);
        return 1;
    }

    path_v paths;
    kv_init(paths);
    if (path_str) {
        kv_push(path_t, paths, make_path_from_str(g, path_str, seq_id? seq_id : 0));
    } else {
        gzFile fp;
        kstring_t s = {0,0,0};
        kstream_t *ks;
        int dret;
        uint64_t lineno = 0;
        char *p, *q, *sid;

        fp = gzopen(path_file, "r");
        if (!fp) {
            fprintf(stderr, "[E::%s] failed to open file %s to read\n", __func__, path_file);
            return 1;
        }
        ks = ks_init(fp);
        while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
            ++lineno;
            
            p = q = s.s;
            while (isspace(*q) && *q != '\0') ++q;
            p = q;
            while (!isspace(*q) && *q != '\0') ++q;
            if (p == q) continue; // empty line
            MYMALLOC(sid, q - p + 1);
            memcpy(sid, p, q - p);
            sid[q-p] = '\0';

            while (isspace(*q) && *q != '\0') ++q;
            p = q;
            while (!isspace(*q) && *q != '\0') ++q;
            if (p == q) {
                // no second column found
                fprintf(stderr, "[E::%s] invalid line at line %ld: %s\n", __func__, (long) lineno, s.s);
                return 1;
            }
            
            kv_push(path_t, paths, make_path_from_str(g, p, sid));
            free(sid);
        }

        ks_destroy(ks);
        gzclose(fp);
    }

    size_t i;
    for (i = 0; i < paths.n; ++i)
        print_seq(g, &paths.a[i], out_seqs? out_seqs : stdout, i, force_linear, 60);
    
    for (i = 0; i < paths.n; ++i)
        path_destroy(&paths.a[i]);
    kv_destroy(paths);

    if (ret) {
        fprintf(stderr, "[E::%s] failed to analysis the GFA file\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (fflush(stdout) == EOF) {
        perror("[E::%s] failed to write the results");
        exit(EXIT_FAILURE);
    }

    if (out_seqs)
        fclose(out_seqs);

    asg_destroy(g);

    if (VERBOSE >= 0) {
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, PATHTOFASTA_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        int i;
        for (i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    return 0;
}

