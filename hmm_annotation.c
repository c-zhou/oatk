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
 * 04/08/22 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
#include <sys/stat.h>
#include <assert.h>

#include "kthread.h"
#include "kvec.h"
#include "kseq.h"

#include "misc.h"
#include "version.h"

KSTREAM_INIT(gzFile, gzread, 65536)

/***************************
 * Run nhmmscan annotation *
 ***************************/

typedef struct {
    uint32_t max_batch_size;
    uint32_t max_batch_num;
    char *nhmmscan;
    char *nhmmdb;
    kstream_t *ks;
    kstring_t *s;
    FILE *fo;
    char *tmpdir;
    int n_threads;
} annot_pipeline_t;

typedef struct {
    annot_pipeline_t *shared;
    int batch_num;
    char **temp_in;
    char **temp_out;
} annot_step_t;

static inline int parse_fheader(char *s)
{
    char *c = s;
    while (*c != 0 && !isspace(*c))
        ++c;
    *c = 0;
    if (c - s == 1) // empty
        return 1;
    return 0;
}

static inline int parse_gseq(char *s, char **seg, char **seq)
{
    if (*s != 'S') return 4;
    
    int i;
    char *p, *q;
    int c;
    for (i = 0, p = q = s + 2;; ++p) {
        if (*p == 0 || *p == '\t') {
            c = *p;
            *p = 0;
            if (i == 0)
                *seg = q;
            else if (i == 1)
                *seq = q;
            ++i, q = p + 1;
            if (i == 2 || c == 0)
                break;
        }
    }

    if (i != 2) return 5;
    if (**seq == '*') return 6;

    return 0;
}

static void annot_worker_for(void *_data, long i, int tid) // kt_for() callback
{
    annot_step_t *annot_s = (annot_step_t *) _data;
    int exit_code, wexit_st;

    char cmd[4096];
    sprintf(cmd, "%s --noali --cpu 1 -o /dev/null --tblout %s %s %s", annot_s->shared->nhmmscan, annot_s->temp_out[i], annot_s->shared->nhmmdb, annot_s->temp_in[i]);

    exit_code = run_system_cmd(cmd, 3);
    if (exit_code == -1) {
        fprintf(stderr, "[E::%s] failed to execute command: %s\n", __func__, cmd);
        exit(EXIT_FAILURE);
    }
    
    wexit_st = WEXITSTATUS(exit_code);
    if (wexit_st)
        fprintf(stderr, "[E::%s] command with non-zero exit code: %d\n", __func__, wexit_st);
    
    return;
}

static void *annot_worker_pipeline(void *shared, int step, void *in)
{
    annot_pipeline_t *p = (annot_pipeline_t *) shared;
    if (step == 0) { // read sequence data into bacthes        
        annot_step_t *annot_s;
        MYCALLOC(annot_s, 1);
        MYMALLOC(annot_s->temp_in, p->max_batch_num);
        MYMALLOC(annot_s->temp_out, p->max_batch_num);

        kstream_t *ks = p->ks;
        kstring_t *s = p->s;
        int dret, ret, is_fa, is_fq, is_gfa;
        uint32_t batch_size, batch_num, n_seq, l_seq;
        FILE *fo;
        char file_template[strlen(p->tmpdir) + 16];

        is_fa = is_fq = is_gfa = 0;
        batch_size = batch_num = 0;
        fo = 0;
        n_seq = l_seq = 0;

        // need the file names only and the file descriptors have already been closed
        annot_s->temp_in[batch_num]  = make_tempfile(p->tmpdir, file_template, ".fa");
        annot_s->temp_out[batch_num] = make_tempfile(p->tmpdir, file_template, ".out");
        
        fo = fopen(annot_s->temp_in[batch_num], "w");
        ++batch_num;

        while (s->l || ks_getuntil(ks, KS_SEP_LINE, s, &dret) >= 0) {
            if (s->l == 0) // empty line
                continue;
            
            if (batch_size >= p->max_batch_size && (is_gfa || (is_fa && s->s[0] == '>') || (is_fq && s->s[0] == '@'))) {
                fclose(fo);
                fo = 0;
                l_seq += batch_size;

                if (batch_num >= p->max_batch_num) {
                    break;
                } else {
                    // make a new batch
                    // need the file names only and the file descriptors have already been closed
                    annot_s->temp_in[batch_num]  = make_tempfile(p->tmpdir, file_template, ".fa");
                    annot_s->temp_out[batch_num] = make_tempfile(p->tmpdir, file_template, ".out");
                    fo = fopen(annot_s->temp_in[batch_num], "w");
                    ++batch_num;
                    batch_size = 0;
                }
            }

            ret = 0; 
            if (!is_gfa && s->s[0] == '>') { // FASTA header
                is_fa = 1;
                // parse header
                ret = parse_fheader(s->s);
                if (!ret)
                    fprintf(fo, "%s\n", s->s);
                ++n_seq;
            } else if (!is_gfa && s->s[0] == '@') { // FASTQ header
                is_fq = 1;
                // parse and write header
                ret = parse_fheader(s->s);
                if (!ret)
                    fprintf(fo, ">%s\n", s->s + 1);
                // parse sequence
                if(ks_getuntil(ks, KS_SEP_LINE, s, &dret) < 0)
                    ret = 2;
                if (!ret) {
                    fprintf(fo, "%s\n", s->s);
                    batch_size += s->l;
                }

                // skip quality score lines
                if (ks_getuntil(ks, KS_SEP_LINE, s, &dret) < 0 ||
                        ks_getuntil(ks, KS_SEP_LINE, s, &dret) < 0)
                    ret = 3;
                ++n_seq;
            } else if (is_fa) { // FASTA sequence line
                fprintf(fo, "%s\n", s->s);
                batch_size += s->l;
            } else {
                is_gfa = 1;
                if (s->s[0] == 'S') {
                    char *seg, *seq;
                    ret = parse_gseq(s->s, &seg, &seq);
                    if (!ret) {
                        fprintf(fo, ">%s\n", seg);
                        fprintf(fo, "%s\n", seq);
                        batch_size += strlen(seq);
                    }
                    ++n_seq;
                }
            }

            if (ret) {
                fprintf(stderr, "[E::%s] failed to parse %s file (error code: %d)\n", __func__, is_fa? "FASTA" : (is_fq? "FASTQ" : "GFA"), ret);
                exit(EXIT_FAILURE);
            }

            s->l = 0;
        }

        if (fo) {
            fclose(fo);
            l_seq += batch_size;
            if (batch_size == 0) {
                --batch_num;
                remove(annot_s->temp_in[batch_num]);
                remove(annot_s->temp_out[batch_num]);
                free(annot_s->temp_in[batch_num]);
                free(annot_s->temp_out[batch_num]);
            }
        }

        assert(is_fa + is_fq + is_gfa <= 1);

        if (batch_num == 0) {
            free(annot_s->temp_in);
            free(annot_s->temp_out);
            free(annot_s);
        } else {
            annot_s->batch_num = batch_num;
            annot_s->shared = p;
            fprintf(stderr, "[M::%s] %u sequences (%u bp) loaded in %u batc%s\n", __func__, n_seq, l_seq, batch_num, batch_num > 1? "hes" : "h");
            return annot_s;
        }
    } else if (step == 1) { // do nhmmscan annotation
        kt_for(p->n_threads, annot_worker_for, in, ((annot_step_t*)in)->batch_num);
        return in;
    } else if (step == 2) { // parse nhmmscan output
        annot_step_t *annot_s = (annot_step_t *) in;
        int i;
        FILE *fp;
        char c;

        for (i = 0; i < annot_s->batch_num; ++i) {
            fp = fopen(annot_s->temp_out[i], "r");
            while ((c = fgetc(fp)) != EOF)
                fputc(c, p->fo);
            fclose(fp);

            remove(annot_s->temp_in[i]);
            remove(annot_s->temp_out[i]);
            free(annot_s->temp_in[i]);
            free(annot_s->temp_out[i]);
        }
        free(annot_s->temp_in);
        free(annot_s->temp_out);
        free(annot_s);
    }

    return 0;
}

int hmm_annotate(char **file_in, int n_file, char *nhmmscan, char *nhmmdb, FILE *fo, uint32_t max_batch_size, uint32_t max_batch_num, int n_threads, char *tmpdir)
{
    annot_pipeline_t pl;
    MYBZERO(&pl, 1);

    pl.max_batch_size = max_batch_size;
    pl.max_batch_num = max_batch_num;
    pl.nhmmscan = nhmmscan? nhmmscan : "nhmmscan";
    pl.nhmmdb = nhmmdb;
    pl.fo = fo;
    pl.n_threads = n_threads;

    int rm_tmpdir = 0;
    if (tmpdir) {
        struct stat st = {0};
        if (stat(tmpdir, &st) == -1) {
            mkdir(tmpdir, 0700);
            rm_tmpdir = 1;
        }
        pl.tmpdir = tmpdir;
    } else {
        char file_template[] = "tmp_XXXXXXXXXX";
        pl.tmpdir = mkdtemp(file_template); // no need to free
        rm_tmpdir = 1;
    }

    int i;
    gzFile fp;
    kstream_t *ks;
    for (i = 0; i < n_file; ++i) {
        fp = file_in[i] && strcmp(file_in[i], "-")? gzopen(file_in[i], "r") : gzdopen(0, "r");
        if (fp == 0) {
            fprintf(stderr, "[E::%s] failed to open file %s to read\n", __func__, file_in[i]);
            exit(EXIT_FAILURE);
        }
        ks = ks_init(fp);
        pl.ks = ks;
        MYCALLOC(pl.s, 1);
        kt_pipeline(n_threads, annot_worker_pipeline, &pl, 3);
        ks_destroy(ks);
        if (pl.s->s) free(pl.s->s);
        gzclose(fp);
    }

    if (rm_tmpdir) rmdir(pl.tmpdir); // should be empty

    return 0;
}

#ifdef ANNOTATION_MAIN
#include "ketopt.h"

int VERBOSE = 0;

static ko_longopt_t long_options[] = {
    { "nhmmscan", ko_required_argument, 301 },
    { "threads",  ko_required_argument, 't' },
    { "version",  ko_no_argument,       'V' },
    { "help",     ko_no_argument,       'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    const char *opt_str = "t:b:T:o:Vv:h";
    ketopt_t opt = KETOPT_INIT;
    int c, ret = 0;
    int n_threads, batch_size, n_file;
    FILE *fp_help, *out_fp;
    char *out, *nhmmdb, *tmpdir;
    char **file_in;

    sys_init();
    
    fp_help = stderr;
    out_fp = stdout;
    out = nhmmdb = tmpdir = 0;
    file_in = 0;
    n_file = 0;
    batch_size = 1000000;
    n_threads = 4;
    char *nhmmscan = "nhmmscan";

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >=0 ) {
        if (c == 't') n_threads = atoi(opt.arg);
        else if (c == 'b') batch_size = atoi(opt.arg);
        else if (c == 'T') tmpdir = opt.arg;
        else if (c == 301) nhmmscan = opt.arg;
        else if (c == 'v') VERBOSE = atoi(opt.arg);
        else if (c == 'h') fp_help = stdout;
        else if (c == 'o') {
            if (strcmp(opt.arg, "-") != 0) {
                out = opt.arg;
            }
        }
        else if (c == 'V') {
            puts(HMMANNOTATION_VERSION);
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
        fprintf(fp_help, "Usage: hmm_annotation [options] <hmmdb> <query>[.gfa|.fa|.fq][.gz] [...]\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "    -b INT           batch size [%d]\n", batch_size);
        fprintf(fp_help, "    -t INT           number threads [4]\n");
        fprintf(fp_help, "    -T STR           temporary directory [NULL]\n");
        fprintf(fp_help, "    -o FILE          output results to FILE [stdout]\n");
        fprintf(fp_help, "    --nhmmscan STR   nhmmscan executable path [%s]\n", nhmmscan);
        fprintf(fp_help, "    -v INT           verbose level [%d]\n", VERBOSE);
        fprintf(fp_help, "    --version        show version number\n");
        fprintf(fp_help, "Example: ./hmm_annotation -t 8 --nhmmscan /usr/bin/nhmmscan -o asm_annot.mito.txt angiosperm_mito.fam asm.gfa\n");
        return fp_help == stdout? 0 : 1;
    }

    if (argc - opt.ind < 1) {
        fprintf(stderr, "[E::%s] missing input: please specify the hmmdb file\n", __func__);
        return 1;
    } else {
        nhmmdb = argv[opt.ind];
    }

    if (argc - opt.ind < 2) {
        fprintf(stderr, "[E::%s] missing input: please specify at least one query file\n", __func__);
        return 1;
    } else {
        file_in = argv + opt.ind + 1;
        n_file = argc - opt.ind - 1;
    }

    // check if nhmmscan executable is available
    check_executable(nhmmscan);

    if (out) out_fp = fopen(out, "w");

    ret = hmm_annotate(file_in, n_file, nhmmscan, nhmmdb, out_fp, batch_size, n_threads * 5, n_threads, tmpdir);

    if (out) fclose(out_fp);

    if (ret) {
        fprintf(stderr, "[E::%s] annotation failed\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (fflush(stdout) == EOF) {
        fprintf(stderr, "[E::%s] failed to write the results\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (VERBOSE >= 0) {
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, HMMANNOTATION_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        int i;
        for (i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    return 0;
}
#endif

