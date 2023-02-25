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
 * 24/02/23 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <sys/stat.h>

#include "kvec.h"
#include "khashl.h"
#include "ketopt.h"

#include "misc.h"
#include "path.h"
#include "graph.h"
#include "syncasm.h"
#include "hmmannot.h"

#define OATK_VERSION "0.1"

int VERBOSE = 0;
double realtime0;

int syncasm(char **file_in, int n_file, int k, int s, int bubble_size, int tip_size, int min_k_cov, int min_a_cov,
        double weak_cross, int do_unzip, int n_threads, char *out, scg_meta_t *meta, int VERBOSE);

int hmm_annotate(char **file_in, int n_file, char *nhmmscan, char *nhmmdb, FILE *fo, uint32_t max_batch_size, 
        uint32_t max_batch_num, int n_threads, char *tmpdir);

int pathfinder(char *asg_file, char *mito_annot, char *pltd_annot, int n_core, int min_len, int min_ex_g, int max_d_len,
        int max_copy, double max_eval, double min_score, double min_cf, double seq_cf, int no_trn, int do_graph_clean, 
        int bubble_size, int tip_size, double weak_cross, int out_s, char *out_pref, int VERBOSE);

static int make_dir(char *dir)
{
    struct stat st = {0};
    if (stat(dir, &st) == -1) {
        int ret = mkdir(dir, 0700);
        if (ret) {
            fprintf(stderr, "[W::%s] create output directory '%s' failed: %s\n", __func__, dir, strerror(errno));
            exit(EXIT_FAILURE);
        }
        return 1;
    }
    return 0;
}

static ko_longopt_t long_options[] = {
    { "mini-circle",    ko_no_argument,       301 },
    { "max-bubble",     ko_required_argument, 302 },
    { "max-tip",        ko_required_argument, 303 },
    { "weak-cross",     ko_required_argument, 304 },
    { "no-unzip",       ko_no_argument,       305 },
    { "nhmmscan",       ko_required_argument, 306 },
    { "longest",        ko_no_argument,       307 },
    { "circular",       ko_no_argument,       308 },
    { "all",            ko_no_argument,       309 },
    { "core-gene",      ko_required_argument, 310 },
    { "min-score",      ko_required_argument, 311 },
    { "min-gain",       ko_required_argument, 312 },
    { "max-eval",       ko_required_argument, 313 },
    { "min-s-length",   ko_required_argument, 314 },
    { "max-d-length",   ko_required_argument, 315 },
    { "min-s-cov",      ko_required_argument, 316 },
    { "max-copy",       ko_required_argument, 317 },
    { "no-graph-clean", ko_no_argument,       318 },
    { "include-trn",    ko_no_argument,       319 },
    { "mito-db",        ko_required_argument, 'm' },
    { "pltd-db",        ko_required_argument, 'p' },
    { "threads",        ko_required_argument, 't' },
    { "version",        ko_no_argument,       'V' },
    { "help",           ko_no_argument,       'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    const char *opt_str = "k:s:c:m:p:b:T:f:o:t:Vv:h";
    ketopt_t opt = KETOPT_INIT;
    int k, s, bubble_size, tip_size, do_unzip, min_k_cov, min_a_cov, batch_size;
    int out_s, out_c, n_db, max_copy, no_trn, n_core, min_len, min_ex_g, max_d_len, do_graph_clean;
    int mini_circle, n_threads;
    double weak_cross, max_eval, min_score, min_cf, seq_cf;
    FILE *fp_help;
    char *out, *nhmmscan, *mito_db, *pltd_db, *tmpdir;
    int c, ret = 0;

    sys_init();

    // public
    n_threads = 1;
    mini_circle = 0;
    out = "./oatk.asm";
    fp_help = stderr;
    // syncasm parameters
    k = 1001;
    s = 31;
    min_k_cov = 3;
    min_a_cov = 1;
    do_unzip = 1;
    bubble_size = 100000;
    tip_size = 10000;
    weak_cross = 0.3;
    // hmm_annotaton parameters
    mito_db = 0;
    pltd_db = 0;
    n_db = 0;
    batch_size = 1000000;
    nhmmscan = "nhmmscan";
    tmpdir = 0;
    // path_finder parameters
    out_s = -1;
    out_c = 0;
    max_copy = 10;
    no_trn = 1;
    do_graph_clean = 1;
    n_core = 10;
    max_eval = 1e-12;
    min_len = 10000;
    min_score = 100;
    min_ex_g = 3;
    max_d_len = 50000;
    seq_cf = .95;
    min_cf = .20;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >=0 ) {
        if (c == 'k') k = atoi(opt.arg);
        else if (c == 's') s = atoi(opt.arg);
        else if (c == 'c') min_k_cov = atoi(opt.arg);
        else if (c == 'm') mito_db = opt.arg, ++n_db;
        else if (c == 'p') pltd_db = opt.arg, ++n_db;
        else if (c == 'b') batch_size = atoi(opt.arg);
        else if (c == 'T') tmpdir = opt.arg;
        else if (c == 'f') seq_cf = atof(opt.arg);
        else if (c == 'o') {
            if (strcmp(opt.arg, "-") != 0)
                out = opt.arg;
        }
        else if (c == 't') n_threads = atoi(opt.arg);
        else if (c == 301) mini_circle = 1;
        else if (c == 302) bubble_size = atoi(opt.arg);
        else if (c == 303) tip_size = atoi(opt.arg);
        else if (c == 304) weak_cross = atof(opt.arg);
        else if (c == 305) do_unzip = 0;
        else if (c == 306) nhmmscan = opt.arg;
        else if (c == 307) out_s = 0, ++out_c;
        else if (c == 308) out_s = 1, ++out_c;
        else if (c == 309) out_s = 2, ++out_c;
        else if (c == 310) n_core = atoi(opt.arg);
        else if (c == 311) min_score = atof(opt.arg);
        else if (c == 312) min_ex_g = atoi(opt.arg);
        else if (c == 313) max_eval = atof(opt.arg);
        else if (c == 314) min_len = atoi(opt.arg);
        else if (c == 315) max_d_len = atoi(opt.arg);
        else if (c == 316) min_cf = atof(opt.arg);
        else if (c == 317) max_copy = atoi(opt.arg);
        else if (c == 318) do_graph_clean = 0;
        else if (c == 319) no_trn = 0;
        else if (c == 'v') VERBOSE = atoi(opt.arg);
        else if (c == 'h') fp_help = stdout;
        else if (c == 'V') {
            puts(OATK_VERSION);
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
        fprintf(fp_help, "Usage: oatk [options] <target.fa[stq][.gz]> [...]\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "  Input/Output:\n");
        fprintf(fp_help, "    -o FILE              prefix of output files [%s]\n", out);
        fprintf(fp_help, "    -t INT               number of threads [%d]\n", n_threads);
        fprintf(fp_help, "    --mini-circle        run mini circle mode such as animal mitochondria and plasmid\n");
        fprintf(fp_help, "    -v INT               verbose level [%d]\n", VERBOSE);
        fprintf(fp_help, "    --version            show version number\n");
        fprintf(fp_help, "  Syncasm:\n");
        fprintf(fp_help, "    -k INT               kmer size [%d]\n", k);
        fprintf(fp_help, "    -s INT               smer size (no larger than 31) [%d]\n", s);
        fprintf(fp_help, "    -c INT               minimum kmer coverage [%d]\n", min_k_cov);
        fprintf(fp_help, "    --max-bubble INT     maximum bubble size for assembly graph clean [%d]\n", bubble_size);
        fprintf(fp_help, "    --max-tip    INT     maximum tip size for assembly graph clean [%d]\n", tip_size);
        fprintf(fp_help, "    --weak-cross FLOAT   maximum relative edge coverage for weak crosslink clean [%.2f]\n", weak_cross);
        fprintf(fp_help, "    --no-unzip           do not run assembly graph unzipping\n");
        fprintf(fp_help, "  Annotation:\n");
        fprintf(fp_help, "    -m FILE              mitochondria gene annotation HMM profile database [NULL]\n");
        fprintf(fp_help, "    -p FILE              plastid gene annotation HMM profile database [NULL]\n");
        fprintf(fp_help, "    -b INT               batch size [%d]\n", batch_size);
        fprintf(fp_help, "    -T STR               temporary directory [NULL]\n");
        fprintf(fp_help, "    --nhmmscan STR       nhmmscan executable path [nhmmscan]\n");
        fprintf(fp_help, "  Pathfinder:\n");
        fprintf(fp_help, "    -f FLOAT             prefer circular path to longest if >= FLOAT sequence covered [%.3f]\n", seq_cf);
        fprintf(fp_help, "    --longest            output only the longest path [default]\n");
        fprintf(fp_help, "    --circular           output only the longest circular path\n");
        fprintf(fp_help, "    --all                output all best paths\n");
        fprintf(fp_help, "    --core-gene INT      number of top core gene annotations to consider [%d]\n", n_core);
        fprintf(fp_help, "    --min-score FLOAT    minimum annotation score of a subgraph [%.1f]\n", min_score);
        fprintf(fp_help, "    --min-gain INT       minimum number of addtional core genes to include a sequence [%d]\n", min_ex_g);
        fprintf(fp_help, "    --max-eval  FLOAT    maximum E-value of a core gene [%.3e]\n", max_eval);
        fprintf(fp_help, "    --min-s-length INT   minimum length of a singleton sequence to keep [%d]\n", min_len);
        fprintf(fp_help, "    --max-d-length INT   maximum length of a singleton sequence to delete [%d]\n", max_d_len);
        fprintf(fp_help, "    --min-s-cov FLOAT    minimum coverage of a sequence compared to the subgraph average [%.3f]\n", min_cf);
        fprintf(fp_help, "    --max-copy INT       maximum copy number to consider [%d]\n", max_copy);
        fprintf(fp_help, "    --no-graph-clean     do not do assembly graph clean\n");
        fprintf(fp_help, "    --include-trn        include TRN type genes\n");
        fprintf(fp_help, "Example: ./oatk -o oatk.asm -t 16 -m angiosperm_mito.fam -p angiosperm_pltd.fam hifi.fa.gz\n");
        return fp_help == stdout? 0 : 1;
    }

    if (out_c > 1) {
        fprintf(stderr, "[E::%s] options --longest, --circular and --all are mutually exclusive\n", __func__);
        exit(EXIT_FAILURE);
    }
    if (out_s < 0) out_s = 0;

    if (n_db == 0) {
        fprintf(stderr, "[E::%s] provide at least one HMM profile database (-m and/or -p)\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (mini_circle && n_db > 1) {
        fprintf(stderr, "[E::%s] only one HMM profile database (-m or -p) allowed for mini-circle mode\n", __func__);
        exit(EXIT_FAILURE);
    }

    /*** parse output dirname and basename ***/
    char *outdir, *outname;
    parse_pathname(out, &outdir, &outname);
    if (!outdir) {
        outdir = ".";
        fprintf(stderr, "[W::%s] invalid output directory - using '%s' instead\n", __func__, outdir);
    } else {
        make_dir(outdir);
    }
    if (!outname) {
        outname = "oatk.asm";
        fprintf(stderr, "[W::%s] invalid output name prefix - using '%s' instead\n", __func__, outname);
    }
    int outlen = strlen(outdir) + strlen(outname) + 1;
    char *outpref;
    MYMALLOC(outpref, outlen);
    sprintf(outpref, "%s/%s", outdir, outname);

    /*** syncasm assembly ***/
    scg_meta_t *scg_meta;
    MYMALLOC(scg_meta, 1);
    ret = syncasm(argv + opt.ind, argc - opt.ind, k, s, bubble_size, tip_size, min_k_cov, min_a_cov, weak_cross, do_unzip, n_threads, outpref, scg_meta, VERBOSE);
    if (ret) {
        fprintf(stderr, "[E::%s] syncasm assembly program failed\n", __func__);
        exit(EXIT_FAILURE);
    }
    char *asg_file;
    MYMALLOC(asg_file, outlen + 16);
    sprintf(asg_file, "%s_utg_final.gfa", outpref);

    /*** hmm annotation ***/
    char *mito_annot, *pltd_annot;
    mito_annot = pltd_annot = 0;
    check_executable(nhmmscan); // check if nhmmscan executable is available
    int rm_tmpdir, free_tmpdir; // make temp dir if not provided or existed
    rm_tmpdir = free_tmpdir = 0;
    if (tmpdir) {
        rm_tmpdir = make_dir(tmpdir);
    } else {
        char *file_template;
        MYMALLOC(file_template, outlen + 16);
        sprintf(file_template, "%s/%s", outdir, "tmp_XXXXXXXXXX");
        tmpdir = mkdtemp(file_template);
        rm_tmpdir = 1;
        free_tmpdir = 1;
    }
    if (mito_db) {    
        FILE *fo;
        MYMALLOC(mito_annot, outlen + 16);
        sprintf(mito_annot, "%s.annot_mito.txt", outpref);
        fo = fopen(mito_annot, "w");
        ret  = hmm_annotate(&asg_file, 1, nhmmscan, mito_db, fo, batch_size, n_threads * 5, n_threads, tmpdir);
        fclose(fo);
    }
    if (pltd_db) {
        FILE *fo;
        MYMALLOC(pltd_annot, outlen + 16);
        sprintf(pltd_annot, "%s.annot_pltd.txt", outpref);
        fo = fopen(pltd_annot, "w");
        ret |= hmm_annotate(&asg_file, 1, nhmmscan, pltd_db, fo, batch_size, n_threads * 5, n_threads, tmpdir);
        fclose(fo);
    }
    if (ret) {
        fprintf(stderr, "[E::%s] annotation program failed\n", __func__);
        exit(EXIT_FAILURE);
    }
    
    /*** pathfinder ***/
    if (mini_circle) // pathfinder in mini-circle mode
        return 1;
    else // pathfinder in normal mode
        ret = pathfinder(asg_file, mito_annot, pltd_annot, n_core, min_len, min_ex_g, max_d_len, max_copy,
                max_eval, min_score, min_cf, seq_cf, no_trn, do_graph_clean, bubble_size, tip_size, weak_cross,
                out_s, outpref, VERBOSE);
    if (ret) {
        fprintf(stderr, "[E::%s] pathfinder program failed\n", __func__);
        exit(EXIT_FAILURE);
    }

    /*** final clean ***/
    if (rm_tmpdir) rmdir(tmpdir); // should be empty
    if (free_tmpdir) free(tmpdir);
    free(mito_annot);
    free(pltd_annot);
    free(asg_file);
    free(outpref);
    scg_meta_destroy(scg_meta);

    if (ret) {
        fprintf(stderr, "[E::%s] oatk program halted\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (fflush(stdout) == EOF) {
        fprintf(stderr, "[E::%s] failed to write the results\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (VERBOSE >= 0) {
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, OATK_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        int i;
        for (i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    return 0;
}

