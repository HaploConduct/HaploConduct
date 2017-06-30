#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random
import subprocess
from time import clock


__author__ = "Jasmijn Baaijens"

usage = """

Pipeline for de novo viralquasispecies assembly.

"""
# fixed settings
#tmp_path = sys.path[0].split('/')
tmp_path = os.path.dirname(os.path.abspath(__file__)).split('/')
selfpath = '/'.join(tmp_path[0:len(tmp_path)-1])
viralquasispecies = selfpath + "/bin/ViralQuasispecies"
COPYFILES = False

def get_original_readcount(fastq):
    count = 0;
    with open(fastq, 'r') as f1:
        for line in f1:
            count += 1
    assert count >= 0
    assert count % 4 == 0
    return int(count/4)

def get_max_subread_id(subreads):
    max_ID = 0
    with open(subreads, 'r') as f:
        for line in f:
            splitline = line.strip('\n').split('\t')
            for subread in splitline[1:]:
                ID = subread.split(':')[0]
                max_ID = max(max_ID, int(ID))
    return max_ID

# GLOBALS
threads = 1
original_readcount = 0
max_read_lengths = []
max_coverages = []
read_counts = []
singles_counts = []
paired_counts = []
overlap_counts = []
edge_counts = []
iteration = 0
verbose = False
stage_a = False
min_read_len = 0
diploid = "false"
max_tip_len = 150
separate_tips = "true"
remove_inclusions = "true"


def main():
    print "pipeline_per_stage.py"

    parser = ArgumentParser(description=usage)
    parser.add_argument('--stage', dest='stage', type=str, required=True,
                            help='specify the algorithm stage (a/b/c)')
    parser.add_argument('--min_overlap_len', dest='min_overlap_len', type=int, default=150)
    parser.add_argument('--min_overlap_perc', dest='min_overlap_perc', type=int, default=0)
    parser.add_argument('--edge_threshold', dest='edge_threshold', type=float, default=0.995)
    parser.add_argument('--merge_contigs', dest='merge_contigs', type=float, default=0)
    parser.add_argument('--fastq', dest='fastq', required=True, type=str)
    parser.add_argument('--overlaps', dest='overlaps', required=True, type=str)
    parser.add_argument('--no_error_correction', dest='error_correction', action='store_false',
                            help='skip error correction in stage a')
    parser.add_argument('--use_subreads', dest='use_subreads', action='store_true')
    parser.add_argument('--num_threads', dest='num_threads', type=int, default=1)
    parser.add_argument('--remove_branches', dest='remove_branches', type=str, default='false')
    parser.add_argument('--min_read_len', dest='min_read_len', type=int, default=0)
    parser.add_argument('--diploid', dest='diploid', action='store_true')
    parser.add_argument('--max_tip_len', dest='max_tip_len', type=int, required=True)
    parser.add_argument('--clique_size_EC', dest='clique_size_EC', type=int, default=4)
    parser.add_argument('--min_overlap_len_EC', dest='min_overlap_len_EC', type=int)
    parser.add_argument('--verbose', dest='verbose', action='store_true')
    args = parser.parse_args()

    FNULL = open(os.devnull, 'w')

    global iteration, original_readcount, max_read_lengths, max_coverages, overlap_counts
    global edge_counts, read_counts, threads, verbose, stage_a, min_read_len, diploid
    global max_tip_len, separate_tips, remove_inclusions, singles_counts, paired_counts

    if args.use_subreads:
        original_readcount = get_max_subread_id("subreads.txt") + 1
    elif args.stage == 'a':
        original_readcount = get_original_readcount(args.fastq + '/paired1.fastq') + get_original_readcount(args.fastq + '/singles.fastq')
    else:
        original_readcount = get_original_readcount(args.fastq + '/singles.fastq')
    if original_readcount == 0:
        sys.stderr.write("Given fastq files are empty. Exiting.\n")
        sys.sterr.flush()
        sys.exit(1)
    if args.stage == 'a':
        read_counts = [original_readcount]
    else:
        read_counts = []
    original_overlaps = analyze_overlaps(args.overlaps)
    overlap_counts = [original_overlaps]
    first_it = "false" if args.use_subreads else "true"
    threads = args.num_threads
    verbose = 'true' if args.verbose else 'false'
    stage_a = True if args.stage == 'a' else False
    min_read_len = args.min_read_len
    diploid = "true" if args.diploid else "false"
    max_tip_len = args.max_tip_len
    separate_tips = "false" if stage_a else "true"
    remove_inclusions = "false" if stage_a else "true" # only keep inclusions during error correction

    # create a global log file; after every iteration the log file is appended to this global log file
    subprocess.call(["rm", "pipeline.log"], stdout=FNULL, stderr=FNULL)
    subprocess.call(["touch", "pipeline.log"])
    # remove existing stats file
    subprocess.call(["rm", "stats.txt"], stdout=FNULL, stderr=FNULL)
    subprocess.call(["touch", "stats.txt"])
    # remove existing tips file
    if not stage_a:
        subprocess.call(["rm", "removed_tip_sequences.fastq"], stdout=FNULL, stderr=FNULL)
        subprocess.call(["touch", "removed_tip_sequences.fastq"])

    min_overlap_len = args.min_overlap_len
    if args.min_overlap_len_EC:
        min_overlap_len_EC = args.min_overlap_len_EC
    else:
        min_overlap_len_EC = args.min_overlap_len
    const_read_its = 0

    if args.stage == 'a':
        # Stage a
        if args.error_correction:
            run_error_correction(args.fastq, args.overlaps, args.edge_threshold, args.min_overlap_perc, min_overlap_len_EC, args.merge_contigs, first_it, args.clique_size_EC)
        else:
            run_first_it_noEC(args.fastq, args.overlaps, args.edge_threshold, args.min_overlap_perc, min_overlap_len, args.merge_contigs, first_it)
        remove_inclusions = "true"
        separate_tips = "true"
        while overlap_counts[-1] > 0 and edge_counts[-1] > 0 and const_read_its < 2:
            while overlap_counts[-1] > 0 and edge_counts[-1] > 0 and const_read_its < 2:
                # merge simple paths
                run_merging_it(args.edge_threshold, args.min_overlap_perc, args.min_overlap_len, 0)
                if read_counts[-1] == read_counts[-2]:
                    const_read_its += 1
                else:
                    const_read_its = 0
            # build super-reads from cliques
            if args.remove_branches == 'false':
                run_clique_it(args.edge_threshold, args.min_overlap_perc, args.min_overlap_len, 0)
                if read_counts[-1] == read_counts[-2]:
                    const_read_its += 1
                else:
                    const_read_its = 0
    #
    elif args.stage == 'b':
        # Stage b
        run_first_it_merge(args.fastq, args.overlaps, args.edge_threshold, args.min_overlap_perc, min_overlap_len, args.merge_contigs, first_it)
        while overlap_counts[-1] > 0 and edge_counts[-1] > 0 and const_read_its < 2:
            while overlap_counts[-1] > 0 and edge_counts[-1] > 0 and const_read_its < 2:
                # merge simple paths
                run_merging_it(args.edge_threshold, args.min_overlap_perc, min_overlap_len, 0)
                if read_counts[-1] == read_counts[-2]:
                    const_read_its += 1
                else:
                    const_read_its = 0
            # merge along branches
            if args.remove_branches == 'false':
                run_clique_it(args.edge_threshold, args.min_overlap_perc, min_overlap_len, 0)
                if read_counts[-1] == read_counts[-2]:
                    const_read_its += 1
                else:
                    const_read_its = 0
    #
    elif args.stage == 'c':
        # Stage c
        run_first_it_merge(args.fastq, args.overlaps, args.edge_threshold, args.min_overlap_perc, min_overlap_len, args.merge_contigs, first_it)
        while overlap_counts[-1] > 0 and edge_counts[-1] > 0 and const_read_its < 2:
            while overlap_counts[-1] > 0 and edge_counts[-1] > 0 and const_read_its < 2:
                # merge simple paths
                run_merging_it(args.edge_threshold, args.min_overlap_perc, min_overlap_len, args.merge_contigs)
                if read_counts[-1] == read_counts[-2]:
                    const_read_its += 1
                else:
                    const_read_its = 0
            # merge along branches
            if args.remove_branches == 'false':
                run_clique_it(args.edge_threshold, args.min_overlap_perc, min_overlap_len, args.merge_contigs)
                if read_counts[-1] == read_counts[-2]:
                    const_read_its += 1
                else:
                    const_read_its = 0
    #
    else:
        sys.stderr.write("ERROR: algorithm stage not properly specified; choose stage a, b, or c.\n")
        sys.stderr.flush()
        sys.exit(1)

    print "Stage %s done in %d iterations" %(args.stage, iteration)
    print "Maximum read length per iteration: \t", max_read_lengths
#    print "Maximum # subreads per iteration: \t", max_coverages
    print "Number of contigs per iteration: \t", singles_counts
    if max(paired_counts) > 0:
        print "Number of paired reads per iteration: \t", paired_counts
    print "Number of overlaps per iteration: \t", overlap_counts


def run_first_it_merge(fastq, overlaps, edge_threshold, min_overlap_perc, min_overlap_len, error_rate, first_it, remove_branches='true'):
    global iteration, max_read_lengths, max_coverages, overlap_counts, edge_counts
    global read_counts, singles_counts, paired_counts
    iteration += 1
    keep_singletons = max(min_overlap_len, min_read_len)
    if verbose == 'true':
        print "\n**************************************"
        print "**** Iteration %d = first_it_merge ****" %iteration
        print "**************************************"
    subprocess.check_call([viralquasispecies,
        "--singles=%s/singles.fastq" %fastq,
        "--overlaps=%s" %overlaps,
        "--threads=%d" %threads,
        "--edge_threshold=%f" %edge_threshold,
        "--first_it=%s" %first_it,
        "--min_clique_size=2",
        "--keep_singletons=%d" %keep_singletons,
        "--remove_branches=%s" %remove_branches,
        "--min_overlap_perc=%d" %min_overlap_perc,
        "--min_overlap_len=%d" %min_overlap_len,
        "--merge_contigs=%f" %error_rate,
        "--FNO=1",
        "--original_readcount=%d" %original_readcount,
        "--error_correction=false",
        "--remove_trans=1",
        "--optimize=false",
        "--verbose=%s" % verbose,
        "--diploid=%s" % diploid,
        "--base_path=%s" % selfpath,
        "--min_read_len=%s" % min_read_len,
        "--max_tip_len=%s" % max_tip_len,
        "--separate_tips=%s" % separate_tips,
        "--ignore_inclusions=%s" % remove_inclusions
    ])
    if COPYFILES:
        copy_files(iteration)
    copy_log()
    [singles_count, paired_count, n_overlaps] = analyze_results()
    readcount = singles_count + paired_count
    read_counts.append(readcount)
    singles_counts.append(singles_count)
    paired_counts.append(paired_count)
    overlap_counts.append(n_overlaps)
    n_edges = get_edge_count()
    edge_counts.append(n_edges)
    if verbose == 'true':
        print "***"


def run_first_it_noEC(fastq, overlaps, edge_threshold, min_overlap_perc, min_overlap_len, error_rate, first_it, remove_branches='true'):
    global iteration, max_read_lengths, max_coverages, overlap_counts, edge_counts
    global read_counts, singles_counts, paired_counts
    iteration += 1
    keep_singletons = max(min_overlap_len, min_read_len)
    if verbose == 'true':
        print "\n**************************************"
        print "**** Iteration %d = first_it_merge ****" %iteration
        print "**************************************"
    subprocess.check_call([viralquasispecies,
        "--singles=%s/singles.fastq" %fastq,
        "--paired1=%s/paired1.fastq" %fastq,
        "--paired2=%s/paired2.fastq" %fastq,
        "--overlaps=%s" %overlaps,
        "--threads=%d" %threads,
        "--edge_threshold=%f" %edge_threshold,
        "--first_it=%s" %first_it,
        "--min_clique_size=2",
        "--keep_singletons=0",
        "--remove_branches=%s" %remove_branches,
        "--min_overlap_perc=%d" %min_overlap_perc,
        "--min_overlap_len=%d" %min_overlap_len,
        "--merge_contigs=%f" %error_rate,
        "--FNO=1",
        "--original_readcount=%d" %original_readcount,
        "--error_correction=false",
        "--remove_trans=1",
        "--optimize=false",
        "--verbose=%s" % verbose,
        "--diploid=%s" % diploid,
        "--base_path=%s" % selfpath,
        "--min_read_len=%s" % min_read_len,
        "--max_tip_len=%s" % max_tip_len,
        "--separate_tips=%s" % separate_tips,
        "--ignore_inclusions=%s" % remove_inclusions
    ])
    if COPYFILES:
        copy_files(iteration)
    copy_log()
    [singles_count, paired_count, n_overlaps] = analyze_results()
    readcount = singles_count + paired_count
    read_counts.append(readcount)
    singles_counts.append(singles_count)
    paired_counts.append(paired_count)
    overlap_counts.append(n_overlaps)
    n_edges = get_edge_count()
    edge_counts.append(n_edges)
    if verbose == 'true':
        print "***"


def run_merging_it(edge_threshold, min_overlap_perc, min_overlap_len, error_rate, remove_branches='true'):
    global iteration, max_read_lengths, max_coverages, overlap_counts, edge_counts
    global read_counts, singles_counts, paired_counts
    iteration += 1
    if verbose == 'true':
        print "\n*******************************"
        print "**** Iteration %d = merging ****" %iteration
        print "*******************************"
    if stage_a:
        paired1 = "paired1.fastq"
        paired2 = "paired2.fastq"
        fno = 1
    else:
        paired1 = "None"
        paired2 = "None"
        fno = 1
    keep_singletons = max(min_overlap_len, min_read_len)
    subprocess.check_call([viralquasispecies,
        "--singles", "singles.fastq",
        "--paired1=%s" %paired1,
        "--paired2=%s" %paired2,
        "--overlaps=%s" %"overlaps.txt",
        "--threads=%d" %threads,
        "--edge_threshold=%f" %edge_threshold,
        "--first_it=false",
        "--keep_singletons=%d" %keep_singletons,
        "--min_clique_size=2",
        "--remove_branches=%s" %remove_branches,
        "--min_overlap_perc=%d" %min_overlap_perc,
        "--min_overlap_len=%d" %min_overlap_len,
        "--merge_contigs=%f" %error_rate,
        "--FNO=%d" %fno,
        "--original_readcount=%d" %original_readcount,
        "--error_correction=false",
        "--remove_trans=1",
        "--optimize=false",
        "--verbose=%s" % verbose,
        "--diploid=%s" % diploid,
        "--base_path=%s" % selfpath,
        "--min_read_len=%s" % min_read_len,
        "--max_tip_len=%s" % max_tip_len,
        "--separate_tips=%s" % separate_tips,
        "--ignore_inclusions=%s" % remove_inclusions
    ])
    if COPYFILES:
        copy_files(iteration)
    copy_log()
    [singles_count, paired_count, n_overlaps] = analyze_results()
    readcount = singles_count + paired_count
    read_counts.append(readcount)
    singles_counts.append(singles_count)
    paired_counts.append(paired_count)
    overlap_counts.append(n_overlaps)
    n_edges = get_edge_count()
    edge_counts.append(n_edges)
    if verbose == 'true':
        print "***"


def run_error_correction(fastq, overlaps, edge_threshold, min_overlap_perc, min_overlap_len, error_rate, first_it, min_clique_size):
    global iteration, max_read_lengths, max_coverages, overlap_counts, edge_counts
    global read_counts, singles_counts, paired_counts
    iteration += 1
    if verbose == 'true':
        print "\n****************************************"
        print "**** Iteration %d = first_it_cliques ****" %iteration
        print "****************************************"
    subprocess.check_call([viralquasispecies,
        "--singles", "%s/singles.fastq" %fastq,
        "--paired1", "%s/paired1.fastq" %fastq,
        "--paired2", "%s/paired2.fastq" %fastq,
        "--overlaps=%s" %overlaps,
        "--threads=%d" %threads,
        "--edge_threshold=%f" %edge_threshold,
        "--first_it=%s" %first_it,
        "--cliques=true",
        "--error_correction=true",
        "--keep_singletons=1000",
        "--min_clique_size=%d" %min_clique_size,
        "--remove_branches=false",
        "--min_overlap_perc=%d" %min_overlap_perc,
        "--min_overlap_len=%d" %min_overlap_len,
        "--merge_contigs=%f" %error_rate,
        "--FNO=3",
        "--original_readcount=%d" %original_readcount,
        "--remove_trans=2",
        "--optimize=false",
        "--verbose=%s" %verbose,
        "--base_path=%s" % selfpath,
        "--min_read_len=%s" % min_read_len,
        "--max_tip_len=%s" % max_tip_len,
        "--separate_tips=%s" % separate_tips,
        "--ignore_inclusions=%s" % remove_inclusions
    ])
    if COPYFILES:
        copy_files(iteration)
    copy_log()
    [singles_count, paired_count, n_overlaps] = analyze_results()
    readcount = singles_count + paired_count
    read_counts.append(readcount)
    singles_counts.append(singles_count)
    paired_counts.append(paired_count)
    overlap_counts.append(n_overlaps)
    n_edges = get_edge_count()
    edge_counts.append(n_edges)
    if verbose == 'true':
        print "***"


def run_clique_it(edge_threshold, min_overlap_perc, min_overlap_len, error_rate):
    global iteration, max_read_lengths, max_coverages, overlap_counts, edge_counts
    global read_counts, singles_counts, paired_counts
    iteration += 1
    if verbose == 'true':
        print "\n*******************************"
        print "**** Iteration %d = cliques ****" %iteration
        print "*******************************"
    if stage_a:
        paired1 = "paired1.fastq"
        paired2 = "paired2.fastq"
    else:
        paired1 = "None"
        paired2 = "None"
    keep_singletons = max(min_overlap_len, min_read_len)
    subprocess.check_call([viralquasispecies,
        "--singles", "singles.fastq",
        "--paired1=%s" %paired1,
        "--paired2=%s" %paired2,
        "--overlaps=%s" %"overlaps.txt",
        "--threads=%d" %threads,
        "--edge_threshold=%f" %edge_threshold,
        "--first_it=false",
        "--cliques=true",
        "--error_correction=false",
        "--keep_singletons=%d" %keep_singletons,
        "--min_clique_size=2",
        "--remove_branches=false",
        "--min_overlap_perc=%d" %min_overlap_perc,
        "--min_overlap_len=%d" %min_overlap_len,
        "--merge_contigs=%f" %error_rate,
        "--FNO=3",
        "--original_readcount=%d" %original_readcount,
        "--remove_trans=1",
        "--optimize=false",
        "--verbose=%s" %verbose,
        "--diploid=%s" % diploid,
        "--base_path=%s" % selfpath,
        "--min_read_len=%s" % min_read_len,
        "--max_tip_len=%s" % max_tip_len,
        "--separate_tips=%s" % separate_tips,
        "--ignore_inclusions=%s" % remove_inclusions
    ])
    if COPYFILES:
        copy_files(iteration)
    copy_log()
    [singles_count, paired_count, n_overlaps] = analyze_results()
    readcount = singles_count + paired_count
    read_counts.append(readcount)
    singles_counts.append(singles_count)
    paired_counts.append(paired_count)
    overlap_counts.append(n_overlaps)
    n_edges = get_edge_count()
    edge_counts.append(n_edges)
    if verbose == 'true':
        print "***"


def copy_files(it):
    subprocess.call(["cp", "singles.fastq", "it%d_singles.fastq" %it])
    subprocess.call(["cp", "overlaps.txt", "it%d_overlaps.txt" %it])
    subprocess.call(["cp", "subreads.txt", "it%d_subreads.txt" %it])
    subprocess.call(["cp", "graph.gfa", "it%d_graph.gfa" %it])
    if stage_a:
        subprocess.call(["cp", "paired1.fastq", "it%d_paired1.fastq" %it])
        subprocess.call(["cp", "paired2.fastq", "it%d_paired2.fastq" %it])

def copy_log():
    subprocess.call("cat viralquasispecies.log >> pipeline.log", shell=True)

def get_edge_count():
    edge_count = -2
    graphfile = 'graph.txt'
    if os.path.isfile(graphfile):
        with open(graphfile, 'r') as f:
            for line in f:
                edge_count +=1
    return edge_count

def analyze_results(cliques=False):
#    max_cov = analyze_coverage()
#    max_coverages.append(max_cov)
    [singles_count, max_len] = analyze_fastq("singles.fastq")
    if os.path.isfile('paired1.fastq'):
        [paired_count, paired_len] = analyze_fastq("paired1.fastq")
    else:
        paired_count = 0
    max_read_lengths.append(max_len)
    if cliques:
        analyze_cliques()

    if os.path.isfile('overlaps.txt'):
        n_overlaps = analyze_overlaps('overlaps.txt')
    else:
        n_overlaps = 0
    return [singles_count, paired_count, n_overlaps]


def analyze_coverage():
    cov_counts = [0 for i in xrange(original_readcount)]
    max_cov = 0

    infile = "subreads.txt"
    if os.path.isfile(infile):
        with open(infile, 'r') as f:
            for line in f:
                reads = line.split()
                cov = len(reads)-1
                cov_counts[cov-1] += 1
                if cov > max_cov:
                    max_cov = cov

#        print "cov\tcount (top 5)"
#        for i in range(max_cov-5, max_cov):
#            if cov_counts[i] != 0:
#                print "%d\t%d" % (i+1, cov_counts[i])
#        print "\n"
        return max_cov
    else:
        return 0


def analyze_fastq(infile):
    len_counts = [0 for i in xrange(100000)]
    max_len = 0

    if os.path.isfile(infile):
        with open(infile, 'r') as f:
            c = 0
            for line in f:
                c += 1
                if c%4 != 2:
                    continue
                seq = line.strip()
                l = len(seq)
                len_counts[l-1] += 1

                if l > max_len:
                    max_len = l
        if verbose == 'true':
            print "longest read: ", max_len
            print "\n"
        return [int(c/4), max_len]
    else:
        return [0, 0]


def analyze_cliques():
    clique_size_counts = [0 for i in xrange(10000)]
    max_size = 0

    infile = "cliques.txt"
    with open(infile, 'r') as f:
        for line in f:
            clique = line.split()
            s = len(clique)
            clique_size_counts[s-1] += 1
            if s > max_size:
                max_size = s
    if verbose:
        print "clique size top 5:"
        print "size\tcount"
        for i in range(max_size-5, max_size):
            if clique_size_counts[i] != 0:
                print "%d\t%d" % (i+1, clique_size_counts[i])
        print "\n"
    return clique_size_counts


def analyze_overlaps(filename):
    pp_count = [0 for i in xrange(4)]
    ps_count = [0 for i in xrange(4)]
    sp_count = [0 for i in xrange(4)]
    ss_count = [0 for i in xrange(4)]
    c = 0
    with open(filename) as f:
        for line in f:
            c += 1
            line = line.strip().split('\t')
            if line[11] == 'p' and line[12] == 'p':
                if line[5] == '-' and line[6] == '+':
                    pp_count[0] += 1
                elif line[5] == '+' and line[6] == '-':
                    pp_count[1] += 1
                elif line[5] == '+' and line[6] == '+':
                    pp_count[2] += 1
                elif line[5] == '-' and line[6] == '-':
                    pp_count[3] += 1
                else:
                    print 'orientation not found...'
            elif line[11] == 'p' and line[12] == 's':
                if line[5] == '-' and line[6] == '+':
                    ps_count[0] += 1
                elif line[5] == '+' and line[6] == '-':
                    ps_count[1] += 1
                elif line[5] == '+' and line[6] == '+':
                    ps_count[2] += 1
                elif line[5] == '-' and line[6] == '-':
                    ps_count[3] += 1
                else:
                    print 'orientation not found...'
            elif line[11] == 's' and line[12] == 'p':
                if line[5] == '-' and line[6] == '+':
                    sp_count[0] += 1
                elif line[5] == '+' and line[6] == '-':
                    sp_count[1] += 1
                elif line[5] == '+' and line[6] == '+':
                    sp_count[2] += 1
                elif line[5] == '-' and line[6] == '-':
                    sp_count[3] += 1
                else:
                    print 'orientation not found...'
            elif line[11] == 's' and line[12] == 's':
                if line[5] == '-' and line[6] == '+':
                    ss_count[0] += 1
                elif line[5] == '+' and line[6] == '-':
                    ss_count[1] += 1
                elif line[5] == '+' and line[6] == '+':
                    ss_count[2] += 1
                elif line[5] == '-' and line[6] == '-':
                    ss_count[3] += 1
                else:
                    print 'orientation not found...'
            else:
                print 'read types not recognized...'

    total = sum(pp_count) + sum(ps_count) + sum(sp_count) + sum(ss_count)
    if verbose == 'true':
        print "Overlaps:"
        print "[-+, +-, ++, --]"
        print "p-p: ", pp_count
        print "p-s: ", ps_count
        print "s-p: ", sp_count
        print "s-s: ", ss_count
        print total
        print "# lines: ", c
        print "\n"
    return total


if __name__ == '__main__':
    sys.exit(main())
