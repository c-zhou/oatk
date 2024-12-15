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
 * 14/12/24 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>

#include "ketopt.h"
#include "kvec.h"
#include "kseq.h"
#include "khashl.h"
#include "sstream.h"
#include "version.h"
#include "misc.h"

KHASHL_MAP_INIT(KH_LOCAL, kh_str_t, kh_str, kh_cstr_t, uint64_t, kh_hash_str, kh_eq_str)

int VERBOSE = 0;

static char comp_table[] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
    32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
    48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', 123, 124, 125, 126, 127
};

static void print_seq(char *seq, int len, int rv, int pos, int line_wd, FILE *fo)
{
    char *s;
    if (!rv) {
        // forward
        for (s = seq, seq += len; s < seq; s++) {
            fputc(*s, fo);
            if ((++pos) % line_wd == 0)
                fputc('\n', fo);
        }
    } else {
         // reverse
        for (s = seq + len - 1; s >= seq; s--) {
            fputc(comp_table[(int) (*s)], fo);
            if ((++pos) % line_wd == 0)
                fputc('\n', fo);
        }
    }
}

static ko_longopt_t long_options[] = {
    { "verbose", ko_required_argument, 'v' },
    { "version", ko_no_argument,       'V' },
    { "help",    ko_no_argument,       'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    const char *opt_str = "l:s:o:rVv:h";
    ketopt_t opt = KETOPT_INIT;
    int c, line_width, absent, ret = 0;
    sstream_t *fin;
    kh_str_t *reg_list;
    khint32_t k;
    FILE *fp_help;
    char *rotate_file, *seq_id, *seq;
    int seq_pos, strand, ln;
    
    sys_init();

    fp_help = stderr;
    line_width = 60;
    rotate_file = seq_id = 0;
    seq_pos = 0;
    strand = 0;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >=0 ) {
        if (c == 's') rotate_file = opt.arg;
        else if (c == 'l') line_width = atoi(opt.arg);
        else if (c == 'r') strand = 1;
        else if (c == 'v') VERBOSE = atoi(opt.arg);
        else if (c == 'h') fp_help = stdout;
        else if (c == 'o') {
            if (strcmp(opt.arg, "-") != 0) {
                if (freopen(opt.arg, "wb", stdout) == NULL) {
                    fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m: %s\n", opt.arg, strerror(errno));
                    return 1;
                }
            }
        }
        else if (c == 'V') {
            puts(SEQROTATE_VERSION);
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
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Usage: rotate [options] <file>[.fa[.gz]] [seq pos]\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "    -s STR        two/three-column rotation file\n");
        fprintf(fp_help, "    -r STR        rotate in reverse strand (effective without -s)\n");
        fprintf(fp_help, "    -l INT        number of residues per line; 0 for 2^31-1 [%d]\n", line_width);
        fprintf(fp_help, "    -o FILE       output results to FILE [stdout]\n");
        fprintf(fp_help, "    -v INT        verbose level [%d]\n", VERBOSE);
        fprintf(fp_help, "    --version     show version number\n");
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Example: ./rotate asm.fa ctg1 100\n\n");
        return fp_help == stdout? 0 : 1;
    }

    if (argc - opt.ind < 1) {
        fprintf(stderr, "[E::%s] missing input: please specify the input FASTA file\n", __func__);
        return 1;
    }

    if (argc - opt.ind > 2) {
        seq_id  = argv[opt.ind + 1];
        seq_pos = atoll(argv[opt.ind + 2]);
        if (seq_pos <= 0) {
            fprintf(stderr, "[E::%s] rotate position must be positive: %d\n", __func__, seq_pos);
            exit(1);
        }
        if (rotate_file) {
            fprintf(stderr, "[W::%s] rotate file %s ignored\n", __func__, rotate_file);
            fprintf(stderr, "[W::%s] using positional parameters: %s %u\n", __func__, seq_id, seq_pos);
        }
    } else {
        if (!rotate_file) {
            fprintf(stderr, "[E::%s] need a file (-s) or two rotation parameters\n", __func__);
            return 1;
        }
    }

    if (line_width == 0) line_width = INT32_MAX;

    reg_list = kh_str_init();
    if (seq_id) {
        k = kh_str_put(reg_list, strdup(seq_id), &absent);
        kh_val(reg_list, k) = (uint64_t) seq_pos << 1 | strand;
    } else {
        gzFile fp;
        kstring_t s = {0,0,0};
        kstream_t *ks;
        int dret;
        int64_t pos, lineno = 0;
        char *p, *q, *sid;

        if (strand)
            fprintf(stderr, "[W::%s] the option '-r' is not effective with '-s' \n", __func__);

        fp = gzopen(rotate_file, "r");
        if (!fp) {
            fprintf(stderr, "[E::%s] failed to open file %s to read\n", __func__, rotate_file);
            exit(1);
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
            pos = atoll(p);
            if (p == q) {
                // no second column found
                fprintf(stderr, "[E::%s] invalid line at line %ld: %s\n", __func__, (long) lineno, s.s);
                fprintf(stderr, "[E::%s] need at least two columns\n", __func__);
                exit(1);
            }
            if (pos <= 0) {
                fprintf(stderr, "[E::%s] rotate position must be positive: %ld\n", __func__, pos);
                exit(1);
            }
            pos <<= 1;

            while (isspace(*q) && *q != '\0') ++q;
            p = q;
            while (!isspace(*q) && *q != '\0') ++q;
            
            if (p != q) {
                if (*p == '-') pos |= 1;
                else if (*p != '+') {
                    fprintf(stderr, "[E::%s] invalid line at line %ld: %s\n", __func__, (long) lineno, s.s);
                    fprintf(stderr, "[E::%s] the third column (strand) must be '+' or '-'\n", __func__);
                    exit(1);
                }
            }
            
            absent = 1;
            k = kh_str_put(reg_list, sid, &absent);
            if (absent)
                kh_val(reg_list, k) = (uint64_t) pos;
            else {
                fprintf(stderr, "[E::%s] invalid line at line %ld: %s\n", __func__, (long) lineno, s.s);
                fprintf(stderr, "[E::%s] duplicate sequence '%s'\n", __func__, sid);
                exit(1);
            }
        }

        ks_destroy(ks);
        gzclose(fp);
    }

    fin = sstream_open(&argv[opt.ind], 1);
    if (fin == 0) {
        fprintf(stderr, "[E::%s] failed to open files: %s\n", __func__, strerror(errno));
        exit(1);
    }
    while ((ln = sstream_read(fin)) >= 0) {
        seq_id = fin->s->ks->name.s;
        seq = fin->s->ks->seq.s;
        k = kh_str_get(reg_list, seq_id);
        fprintf(stdout, ">%s\n", seq_id);
        if (k != kh_end(reg_list)) {
            seq_pos = (kh_val(reg_list, k)) >> 1;
            strand = (kh_val(reg_list, k)) & 1;
            if (seq_pos > ln) {
                fprintf(stderr, "[E::%s] rotation position (%d) larger than sequence length (%d)\n", __func__, seq_pos, ln);
                exit(1);
            }
            if (strand) {
                print_seq(seq, seq_pos, 1, 0, line_width, stdout);
                print_seq(seq + seq_pos, ln - seq_pos, 1, seq_pos, line_width, stdout);
            } else {
                seq_pos -= 1;
                print_seq(seq + seq_pos, ln - seq_pos, 0, 0, line_width, stdout);
                print_seq(seq, seq_pos, 0, ln - seq_pos, line_width, stdout);
            }
            free((void *)kh_key(reg_list, k));
            kh_str_del(reg_list, k);
        } else
            print_seq(seq, ln, 0, 0, line_width, stdout);
        if (ln % line_width != 0) fputc('\n', stdout);
    }
    sstream_close(fin);

    for (k = (khint32_t) 0; k < kh_end(reg_list); k++) {
        if (kh_exist(reg_list, k)) {
            seq_id = (char *) kh_key(reg_list, k);
            fprintf(stderr, "[W::%s] sequence '%s' not found in the FASTA file\n", __func__, seq_id);
            free(seq_id);
        }
    }
    kh_str_destroy(reg_list);

    if (fflush(stdout) == EOF) {
        fprintf(stderr, "[E::%s] failed to write the results\n", __func__);
        exit(1);
    }

    if (VERBOSE >= 0 && !ret) {
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, SEQROTATE_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        int i;
        for (i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    return 0;
}

