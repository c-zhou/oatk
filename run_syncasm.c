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
#include <math.h>
#include <errno.h>
#include <limits.h>

#include "kvec.h"
#include "kstring.h"

#include "misc.h"
#include "sstream.h"
#include "syncmer.h"
#include "syncasm.h"
#include "graph.h"
#include "version.h"

#undef DEBUG_SYNCMER_SEQ
#undef DEBUG_SYNCMER_GRAPH
#undef DEBUG_GRAPH_UNITIG
#undef DEBUG_GRAPH_MULTIPLEX
#undef DEBUG_GRAPH_ALIGNMENT
#undef DEBUG_GRAPH_ERROR_CORRECTION

void read_error_correction(sr_db_t *sr_db, scg_t *g, double max_edist, uint32_t err_mer_c, uint32_t max_err_c,
        uint32_t err_arc_c, double max_arc_f, int threads, FILE *fo, int verbose);

int syncasm(char **file_in, int n_file, size_t m_data, int k, int s, int bubble_size, int tip_size, int min_k_cov, double min_a_cov_f,
        double weak_cross, int do_ec, int do_unzip, int n_threads, char *out, scg_meta_t *meta, int VERBOSE)
{
    FILE *fo;
    sstream_t *sr_rdr;
    scg_t *scg;
    sr_db_t *sr_db;
    syncmer_db_t *scm_db;
    scg_ra_v *ra_db;
    int ret = 0;

    scg = 0;
    sr_db = 0;
    scm_db = 0;
    ra_db = 0;

    sr_rdr = sstream_open(file_in, n_file);
    if (sr_rdr == 0) {
        fprintf(stderr, "[E::%s] failed to open files: %s\n", __func__, strerror(errno));
        ret = 1;
        goto do_clean;
    }

    MYMALLOC(sr_db, 1);
    sr_db_init(sr_db, k, s);
    sr_read(sr_rdr, sr_db, m_data, n_threads);
    fprintf(stderr, "[M::%s] collected syncmers from %lu target sequence(s)\n", __func__, sr_rdr->n_seq);
    sstream_close(sr_rdr);
    if (sr_db_validate(sr_db)) {
        ret = 1;
        goto do_clean;
    }
    sr_db_stat(sr_db, stderr, VERBOSE);

    if (min_k_cov == 0) {
        min_k_cov = sr_db->stats->kmer_peak_het > 0? (sr_db->stats->kmer_peak_het * 10) : (sr_db->stats->kmer_peak_hom * 10);
        fprintf(stderr, "[M::%s] set minimum kmer coverage as %d\n", __func__, min_k_cov);
    }

#ifdef DEBUG_SYNCMER_SEQ
    uint32_t i;
    fo = open_outstream(out, "_syncmer_debug.fa");
    for (i = 0; i < sr_db->n; ++i) print_all_syncmers_on_seq(&sr_db->a[i], s, k, fo);
    fclose(fo);
#endif

    // make syncmer database
    scm_db = collect_syncmer_from_reads(sr_db);
    
    // syncmer_link_coverage_analysis(sr_db, scm_db, min_k_cov, 30, 30, .7, 0, 0, 0, VERBOSE);
    
    if (do_ec) {
        // make syncmer graph with all syncmers for error correction
        scg = make_syncmer_graph(sr_db, scm_db, 0, 0.);
#ifdef DEBUG_GRAPH_ERROR_CORRECTION
        // save sequence in graph
        fo = open_outstream(out, "_syncmer_hoco.gfa");
        scg_consensus(sr_db, scg, 1, 1, fo);
        fclose(fo);
#endif
        // do consensus in hoco space
        scg_consensus(sr_db, scg, 1, 1, 0);
        // naive error finding without checking arc coverage
#ifdef DEBUG_GRAPH_ERROR_CORRECTION
        fo = open_outstream(out, ".ec.fa");
#else
        fo = 0;
#endif
        read_error_correction(sr_db, scg, 0.02, min_k_cov, min_k_cov * 10, min_k_cov, min_a_cov_f, n_threads, fo, VERBOSE);
#ifdef DEBUG_GRAPH_ERROR_CORRECTION
        fclose(fo);
        fo = open_outstream(out, "_syncmer_hoco.noerr.gfa");
        scg_print(scg, fo, 0);
        fclose(fo);
#endif
        sr_db_stat(sr_db, stderr, VERBOSE);
        scg_destroy(scg); scg = 0;
        // goto do_clean;
    }

    // make syncmer graph
    fprintf(stderr, "[M::%s] make syncmer graph\n", __func__);
    scg = make_syncmer_graph(sr_db, scm_db, min_k_cov, min_a_cov_f);
    if (!scg || scg_is_empty(scg)) {
        fprintf(stderr, "[E::%s] empty syncmer graph\n", __func__);
        ret = 1;
        goto do_clean;
    }
    fprintf(stderr, "[M::%s] syncmer graph stats\n", __func__);
    scg_stat(scg, stderr, 0);
    if (VERBOSE > 1) scg_subgraph_stat(scg, stderr);

#ifdef DEBUG_SYNCMER_GRAPH
    fo = open_outstream(out, "_syncmer.gfa");
    scg_consensus(sr_db, scg, 0, 0, fo);
    fclose(fo);
    MYCALLOC(ra_db, 1);
    scg_read_alignment(sr_db, ra_db, scg, n_threads, 0);
    fprintf(stderr, "[DEBUG_SYNCMER_GRAPH::%s] read alignment\n", __func__);
    scg_rv_print(ra_db, stderr);
    scg_ra_v_destroy(ra_db);
#endif

    // make unitigs
    fprintf(stderr, "[M::%s] syncmer graph unitigging\n", __func__);
    process_mergeable_unitigs(scg);
    fprintf(stderr, "[M::%s] syncmer graph stats after unitigging\n", __func__);
    scg_stat(scg, stderr, 0);
    fo = open_outstream(out, ".utg.gfa");
    scg_consensus(sr_db, scg, 0, 0, fo);
    fclose(fo);

    if (VERBOSE > 1) scg_subgraph_stat(scg, stderr);

#ifdef DEBUG_GRAPH_UNITIG
    scg_print_unitig_syncmer_list(scg, stderr);
    MYCALLOC(ra_db, 1);
    scg_read_alignment(sr_db, ra_db, scg, n_threads, 0);
    fprintf(stderr, "[DEBUG_GRAPH_UNITIG::%s] read alignment\n", __func__);
    scg_rv_print(ra_db, stderr);
    scg_ra_v_destroy(ra_db);
#endif

    // do basic cleanup
    // already have consensus information
    fprintf(stderr, "[M::%s] syncmer graph cleanup\n", __func__);
    uint64_t cleaned = 1;
    while (cleaned) {
        // do not do bubble popping before unzipping to avoid removing haplotypes
        cleaned = 0;
        if (do_unzip <= 0) {
            cleaned += asmg_pop_bubble(scg->utg_asmg, bubble_size, 0, 0, 1, 0, VERBOSE);
            cleaned += asmg_remove_weak_crosslink(scg->utg_asmg, weak_cross, 10, 0, VERBOSE);
        }
        cleaned += asmg_drop_tip(scg->utg_asmg, INT_MAX, tip_size, 1, 0, VERBOSE);
    }
    process_mergeable_unitigs(scg);
    
#ifdef DEBUG_GRAPH_MULTIPLEX
    fprintf(stderr, "[M::%s] syncmer graph stats after cleanup\n", __func__);
    scg_stat(scg, stderr, 0);
    fo = open_outstream(out, ".utg.clean.gfa");
    scg_consensus(sr_db, scg, 0, 0, fo);
    fclose(fo);
#endif

#ifdef DEBUG_GRAPH_MULTIPLEX
    scg_print_unitig_syncmer_list(scg, stderr);
#endif

    MYCALLOC(ra_db, 1);
    // do read threading
    if (do_unzip > 0) {
        fprintf(stderr, "[M::%s] assembly graph unzipping\n", __func__);
        int round, updated;
        uint32_t max_n_scm;

        // set max_n_scm for maximum repeats of ~15Kb which reach the HiFi read length limit
        max_n_scm = ceil(30000.0 / k);

        // do multiplexing
        round = 0;
        updated = 1;
        while (updated != 0 && round < do_unzip) {
            ++round;
            scg_read_alignment(sr_db, ra_db, scg, n_threads, 1);
            // scg_rv_print(ra_db, stderr);
            scg_update_utg_cov(scg);
            updated = scg_multiplex(scg, ra_db, max_n_scm, 10, .3);
            if (VERBOSE > 0) {
                fprintf(stderr, "[M::%s] syncmer graph stats after multiplexing round %d\n", __func__, round);
                scg_stat(scg, stderr, 0);
            }
#ifdef DEBUG_GRAPH_MULTIPLEX
            char *out1;
            MYMALLOC(out1, strlen(out) + 36);
            sprintf(out1, "%s.utg.unzip.r%02d.gfa", out, round);
            fo = open_outstream(out1, "");
            scg_consensus(sr_db, scg, 0, 0, fo);
            fclose(fo);
            free(out1);
            scg_print_unitig_syncmer_list(scg, stderr);
#endif
        }
        
        // arc coverage estimation from aligned reads
        // to remove weak cross arcs
        // only arc coverage is required
        scg_read_alignment(sr_db, ra_db, scg, n_threads, 1);
        scg_ra_arc_coverage(scg, sr_db, ra_db, 0, VERBOSE);
        asmg_remove_weak_crosslink(scg->utg_asmg, weak_cross, 10, 0, VERBOSE);

#ifdef DEBUG_GRAPH_MULTIPLEX
        fprintf(stderr, "[M::%s] syncmer graph stats after multiplexing\n", __func__);
        scg_stat(scg, stderr, 0);
        fo = open_outstream(out, ".utg.multiplex.gfa");
        scg_consensus(sr_db, scg, 0, 0, fo);
        fclose(fo);
#endif

        // do demultiplexing
        scg_demultiplex(scg);
        // the arc coverage is lost
        scg_read_alignment(sr_db, ra_db, scg, n_threads, 0);
        scg_ra_utg_coverage(scg, sr_db, ra_db, VERBOSE);
        scg_ra_arc_coverage(scg, sr_db, ra_db, 1, VERBOSE);
        
#ifdef DEBUG_GRAPH_MULTIPLEX
        fprintf(stderr, "[M::%s] syncmer graph stats after unzipping\n", __func__);
        scg_stat(scg, stderr, 0);
        fo = open_outstream(out, ".utg.unzip.gfa");
        scg_consensus(sr_db, scg, 0, 0, fo);
        fclose(fo);
#else
        // consensus infomration is required for basic cleanup
        scg_consensus(sr_db, scg, 0, 0, 0);
#endif

        // do basic cleanup
        cleaned = 1;
        while (cleaned) {
            cleaned = 0;
            cleaned += asmg_pop_bubble(scg->utg_asmg, bubble_size, 0, 0, 1, 0, VERBOSE);
            cleaned += asmg_remove_weak_crosslink(scg->utg_asmg, weak_cross, 10, 0, VERBOSE);
            cleaned += asmg_drop_tip(scg->utg_asmg, INT_MAX, tip_size, 1, 0, VERBOSE);
        }
        process_mergeable_unitigs(scg);

#ifdef DEBUG_GRAPH_UNITIG
        scg_print_unitig_syncmer_list(scg, stderr);
#endif

#ifdef DEBUG_GRAPH_ALIGNMENT
        fprintf(stderr, "[DEBUG_GRAPH_ALIGNMENT::%s] read alignment\n", __func__);
        scg_rv_print(ra_db, stderr);
#endif
    }

    // unitig and arc coverage estimation
    scg_read_alignment(sr_db, ra_db, scg, n_threads, 0);
    scg_ra_utg_coverage(scg, sr_db, ra_db, VERBOSE);
    scg_ra_arc_coverage(scg, sr_db, ra_db, 1, VERBOSE);

    fprintf(stderr, "[M::%s] syncmer graph stats after final processing\n", __func__);
    scg_stat(scg, stderr, 0);
    fo = open_outstream(out, ".utg.final.gfa");
    scg_consensus(sr_db, scg, 0, 0, fo);
    fclose(fo);

do_clean:
    if (meta) {
        scg_meta_clean(meta);
        meta->k = k;
        meta->s = s;
        meta->scg = scg;
        meta->scm_db = scm_db;
        meta->sr_db = sr_db;
        meta->ra_db = ra_db;
    } else {
        scg_destroy(scg);
        syncmer_db_destroy(scm_db);
        sr_db_destroy(sr_db);
        scg_ra_v_destroy(ra_db);
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
    { "unzip-round",ko_required_argument, 304 },
    { "no-read-ec", ko_no_argument,       305 },
    { "threads",    ko_required_argument, 't' },
    { "verbose",    ko_required_argument, 'v' },
    { "version",    ko_no_argument,       'V' },
    { "help",       ko_no_argument,       'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    const char *opt_str = "k:s:c:a:D:t:v:o:Vh";
    ketopt_t opt = KETOPT_INIT;
    int c, k, s, bubble_size, tip_size, min_k_cov, n_threads;
    size_t m_data;
    double min_a_cov_f, weak_cross;
    char *out;
    int do_ec, do_unzip;
    FILE *fp_help = stderr;
    int ret = 0;

    sys_init();

    k = 1001;
    s = 31;
    min_k_cov = 3;
    min_a_cov_f = .35;
    n_threads = 1;
    bubble_size = 100000;
    tip_size = 10000;
    weak_cross = 0.3;
    m_data = 0;
    do_ec = 1;
    do_unzip = 3;
    out = "syncasm.asm";

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0) {
        if (c == 'k') k = atoi(opt.arg);
        else if (c == 's') s = atoi(opt.arg);
        else if (c == 'c') min_k_cov = atoi(opt.arg);
        else if (c == 'a') min_a_cov_f = atof(opt.arg);
        else if (c == 'D') {
            char *q;
            m_data = strtol(opt.arg, &q, 0);
            if (*q == 'k' || *q == 'K') m_data <<= 10;
            else if (*q == 'm' || *q == 'M') m_data <<= 20;
            else if (*q == 'g' || *q == 'G') m_data <<= 30;
        }
        else if (c == 't') n_threads = atoi(opt.arg);
        else if (c == 301) bubble_size = atoi(opt.arg);
        else if (c == 302) tip_size = atoi(opt.arg);
        else if (c == 303) weak_cross = atof(opt.arg);
        else if (c == 304) do_unzip = atoi(opt.arg);
        else if (c == 305) do_ec = 0;
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
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Usage: syncasm [options] <target.fa[stq][.gz]> [...]\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "    -k INT               kmer size [%d]\n", k);
        fprintf(fp_help, "    -s INT               smer size (no larger than 31) [%d]\n", s);
        fprintf(fp_help, "    -c INT               minimum kmer coverage [%d]\n", min_k_cov);
        fprintf(fp_help, "    -a FLOAT             minimum arc coverage [%.2f]\n", min_a_cov_f);
        fprintf(fp_help, "    -D INT               maximum amount of data to use; suffix K/M/G recognized [%lu]\n", m_data);
        fprintf(fp_help, "    -t INT               number of threads [%d]\n", n_threads);
        fprintf(fp_help, "    -o FILE              prefix of output files [%s]\n", out);
        fprintf(fp_help, "    --max-bubble  INT    maximum bubble size for assembly graph clean [%d]\n", bubble_size);
        fprintf(fp_help, "    --max-tip     INT    maximum tip size for assembly graph clean [%d]\n", tip_size);
        fprintf(fp_help, "    --weak-cross  FLOAT  maximum relative edge coverage for weak crosslink clean [%.2f]\n", weak_cross);
        fprintf(fp_help, "    --unzip-round INT    maximum round of assembly graph unzipping [%d]\n", do_unzip);
        fprintf(fp_help, "    --no-read-ec         do not do read error correction\n");
        fprintf(fp_help, "    -v INT               verbose level [%d]\n", VERBOSE);
        fprintf(fp_help, "    --version            show version number\n");
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Example: ./syncasm -k 1001 -c 50 -t 8 -o syncasm.asm hifi.fa.gz\n\n");
        return fp_help == stdout? 0 : 1;
    }

    ret = syncasm(argv + opt.ind, argc - opt.ind, m_data, k, s, bubble_size, tip_size, min_k_cov, min_a_cov_f, weak_cross, do_ec, do_unzip, n_threads, out, 0, VERBOSE);

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
