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
#include <string.h>
#include <errno.h>

#include "kvec.h"
#include "kstring.h"

#include "misc.h"
#include "sstream.h"
#include "syncasm.h"
#include "syncmer.h"
#include "graph.h"
#include "version.h"

#undef DEBUG_SYNCMER_SEQ
#undef DEBUG_SYNCMER_GRAPH
#undef DEBUG_GRAPH_MULTIPLEX
#undef DEBUG_GRAPH_ALIGNMENT

int syncasm(char **file_in, int n_file, int k, int s, int bubble_size, int tip_size, int min_k_cov, int min_a_cov,
        double weak_cross, int do_unzip, int n_threads, char *out, scg_meta_t *meta, int VERBOSE)
{
    FILE *fo;
    sstream_t *sr_rdr;
    scg_t *scg;
    sr_v *sr;
    scg_ra_v *ra;
    int ret = 0;

    scg = 0;
    sr = 0;
    ra = 0;

    sr_rdr = sstream_open(file_in, n_file);
    if (sr_rdr == 0) {
        fprintf(stderr, "[E::%s] failed to open files: %s\n", __func__, strerror(errno));
        ret = 1;
        goto do_clean;
    }

    MYCALLOC(sr, 1);
    sr_read(sr_rdr, sr, s, k, n_threads);
    fprintf(stderr, "[M::%s] collected syncmers from %lu target sequence(s)\n", __func__, sr_rdr->n_seq);
    if (sr_validate(sr)) {
        ret = 1;
        goto do_clean;
    }

    sr_stat_t *stats;
    MYCALLOC(stats, 1);
    sr_stat(sr, stats, k, stderr, VERBOSE);

    if (min_k_cov == 0) {
        min_k_cov = stats->kmer_peak_het > 0? (stats->kmer_peak_het * 10) : (stats->kmer_peak_hom * 10);
        min_a_cov = MAX(min_k_cov / 3, 1);
        fprintf(stderr, "[M::%s] set minimum kmer coverage as %d\n", __func__, min_k_cov);
        fprintf(stderr, "[M::%s] set minimum edge coverage as %d\n", __func__, min_a_cov);
    }

    free(stats);
    sstream_close(sr_rdr);

#ifdef DEBUG_SYNCMER_SEQ
    uint32_t i;
    fo = open_outstream(out, "_syncmer_debug.fa");
    for (i = 0; i < sr->n; ++i) print_all_syncmers_on_seq(&sr->a[i], s, k, fo);
    fclose(fo);
#endif

    // make syncmer graph
    fprintf(stderr, "[M::%s] make syncmer graph\n", __func__);
    scg = make_syncmer_graph(sr, min_k_cov, min_a_cov);
    if (!scg || scg_is_empty(scg)) {
        fprintf(stderr, "[E::%s] empty syncmer graph\n", __func__);
        ret = 1;
        goto do_clean;
    }
    fprintf(stderr, "[M::%s] syncmer graph stats\n", __func__);
    scg_stat(scg, stderr, 0);
#ifdef DEBUG_SYNCMER_GRAPH
    fo = open_outstream(out, "_syncmer.gfa");
    scg_consensus(sr, scg, k, 0, fo);
    asmg_clean_consensus(scg->utg_asmg);
    fclose(fo);
#endif

    // make unitigs
    fprintf(stderr, "[M::%s] syncmer graph unitigging\n", __func__);
    process_mergeable_unitigs(scg);
    fprintf(stderr, "[M::%s] syncmer graph stats after unitigging\n", __func__);
    scg_stat(scg, stderr, 0);
    fo = open_outstream(out, "_utg.gfa");
    scg_consensus(sr, scg, k, 0, fo);
    fclose(fo);

    // do basic cleanup
    // already have consensus information
    asmg_drop_tip(scg->utg_asmg, 1, tip_size, 0, VERBOSE);
    asmg_pop_bubble(scg->utg_asmg, bubble_size, 0, 0, 1, 0, VERBOSE);
    asmg_remove_weak_crosslink(scg->utg_asmg, weak_cross, 0, VERBOSE);
    asmg_clean_consensus(scg->utg_asmg); // need to clean consensus before unitigging
    process_mergeable_unitigs(scg);
    fprintf(stderr, "[M::%s] syncmer graph stats after cleanup\n", __func__);
    scg_stat(scg, stderr, 0);
    fo = open_outstream(out, do_unzip? "_utg_clean.gfa" : "_utg_final.gfa");
    scg_consensus(sr, scg, k, 0, fo);
    fclose(fo);

    // do read threading
    MYCALLOC(ra, 1);
    if (do_unzip) {
        fprintf(stderr, "[M::%s] assembly graph unzipping\n", __func__);
        int round, updated;

        // do multiplexing
        round = 0;
        updated = 1;
        asmg_clean_consensus(scg->utg_asmg); // need to clean consensus information
        while (updated != 0 && round < 10) {
            scg_read_alignment(sr, ra, scg, n_threads, 1);
            updated = scg_multiplex(scg, ra, 3);
            if (VERBOSE > 0) {
                fprintf(stderr, "[M::%s] syncmer graph stats after multiplexing round %d\n", __func__, ++round);
                scg_stat(scg, stderr, 0);
            }
#ifdef DEBUG_GRAPH_MULTIPLEX
            scg_print_unitig_syncmer_list(scg, stderr);
#endif
        }

        // do demultiplexing
        scg_demultiplex(scg);
        // the arc coverage is lost
        scg_arc_coverage(scg, sr);
        fprintf(stderr, "[M::%s] syncmer graph stats after unzipping\n", __func__);
        scg_stat(scg, stderr, 0);
        fo = open_outstream(out, "_utg_unzip.gfa");
        scg_consensus(sr, scg, k, 0, fo);
        fclose(fo);

        // do basic cleanup
        // already have consensus information
        asmg_drop_tip(scg->utg_asmg, 1, tip_size, 0, VERBOSE);
        asmg_pop_bubble(scg->utg_asmg, bubble_size, 0, 0, 1, 0, VERBOSE);
        asmg_remove_weak_crosslink(scg->utg_asmg, weak_cross, 0, VERBOSE);
        asmg_clean_consensus(scg->utg_asmg); // need to clean consensus before unitigging
        process_mergeable_unitigs(scg);
        fprintf(stderr, "[M::%s] syncmer graph stats after final cleanup\n", __func__);
        scg_stat(scg, stderr, 0);
        fo = open_outstream(out, "_utg_final.gfa");
        scg_consensus(sr, scg, k, 0, fo);
        fclose(fo);

#ifdef DEBUG_GRAPH_ALIGNMENT
        asmg_clean_consensus(scg->utg_asmg); // need to clean consensus before alignment
        scg_read_alignment(sr, ra, scg, n_threads, 0);
        fprintf(stderr, "[DEBUG_GRAPH_ALIGNMENT::%s] read alignment\n", __func__);
        scg_rv_print(ra, stderr);
#endif
    }

do_clean:
    if (meta) {
        scg_meta_clean(meta);
        meta->k = k;
        meta->s = s;
        meta->scg = scg;
        meta->sr = sr;
        meta->ra = ra;
    } else {
        scg_destroy(scg);
        sr_v_destroy(sr);
        scg_ra_v_destroy(ra);
    }

    return ret;
}

#ifdef SYNCASM_MAIN
#include "ketopt.h"

int VERBOSE = 0;

static ko_longopt_t long_options[] = {
    { "max-bubble", ko_required_argument, 301 },
    { "max-tip",    ko_required_argument, 302 },
    { "weak-cross", ko_required_argument, 303 },
    { "no-unzip",   ko_no_argument,       304 },
    { "threads",    ko_required_argument, 't' },
    { "version",    ko_no_argument,       'V' },
    { "help",       ko_no_argument,       'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    const char *opt_str = "k:s:c:t:v:o:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c, k, s, bubble_size, tip_size, do_unzip, min_k_cov, min_a_cov, n_threads;
    double weak_cross;
    char *out;
    FILE *fp_help = stderr;
    int ret = 0;

    sys_init();

    k = 1001;
    s = 31;
    min_k_cov = 3;
    min_a_cov = 1;
    n_threads = 1;
    do_unzip = 1;
    bubble_size = 100000;
    tip_size = 10000;
    weak_cross = 0.3;
    out = "syncasm.asm";

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0) {
        if (c == 'k') k = atoi(opt.arg);
        else if (c == 's') s = atoi(opt.arg);
        else if (c == 'c') min_k_cov = atoi(opt.arg), min_a_cov = MAX(min_k_cov / 3, 1);
        else if (c == 't') n_threads = atoi(opt.arg);
        else if (c == 301) bubble_size = atoi(opt.arg);
        else if (c == 302) tip_size = atoi(opt.arg);
        else if (c == 303) weak_cross = atof(opt.arg);
        else if (c == 304) do_unzip = 0;
        else if (c == 'o') {
            if (strcmp(opt.arg, "-") != 0)
                out = opt.arg;
        }
        else if (c == 'v') VERBOSE = atoi(opt.arg);
        else if (c == 'V') {
            puts(SYNCASM_VERSION);
            return 0;
        }
        else if (c == 'h') fp_help = stdout;
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
        fprintf(fp_help, "Usage: syncasm [options] <target.fa[stq][.gz]> [...]\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "    -k INT               kmer size [%d]\n", k);
        fprintf(fp_help, "    -s INT               smer size (no larger than 31) [%d]\n", s);
        fprintf(fp_help, "    -c INT               minimum kmer coverage [%d]\n", min_k_cov);
        fprintf(fp_help, "    -t INT               number of threads [%d]\n", n_threads);
        fprintf(fp_help, "    -o FILE              prefix of output files [%s]\n", out);
        fprintf(fp_help, "    --max-bubble INT     maximum bubble size for assembly graph clean [%d]\n", bubble_size);
        fprintf(fp_help, "    --max-tip    INT     maximum tip size for assembly graph clean [%d]\n", tip_size);
        fprintf(fp_help, "    --weak-cross FLOAT   maximum relative edge coverage for weak crosslink clean [%.2f]\n", weak_cross);
        fprintf(fp_help, "    --no-unzip           do not run assembly graph unzipping\n");
        fprintf(fp_help, "    -v INT               verbose level [%d]\n", VERBOSE);
        fprintf(fp_help, "    --version            show version number\n");
        fprintf(fp_help, "Example: ./syncasm -k 1001 -c 50 -t 8 -o syncasm.asm hifi.fa.gz\n");
        return fp_help == stdout? 0 : 1;
    }

    ret = syncasm(argv + opt.ind, argc - opt.ind, k, s, bubble_size, tip_size, min_k_cov, min_a_cov, weak_cross, do_unzip, n_threads, out, 0, VERBOSE);

    if (ret) {
        fprintf(stderr, "[E::%s] failed to constrcut assembly\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (fflush(stdout) == EOF) {
        fprintf(stderr, "[E::%s] failed to write the results\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (VERBOSE >= 0) {
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, SYNCASM_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        int i;
        for (i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    return 0;
}

#endif
