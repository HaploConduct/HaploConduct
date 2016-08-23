#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random
import subprocess
from time import clock


__author__ = "Jasmijn Baaijens"

usage = """%prog [options]

Pipeline for de novo viralquasispecies assembly. 

"""
# fixed settings
THREADS = 12
tmp_path = sys.path[0].split('/')
selfpath = '/'.join(tmp_path[0:len(tmp_path)-1])
viralquasispecies = selfpath + "/bin/ViralQuasispecies"

def get_original_readcount(fastq):
    count = 0;
    with open("%s/singles.fastq" %fastq, 'r') as f1:
        for line in f1:
            count += 1
    with open("%s/paired1.fastq" %fastq, 'r') as f2:
        for line in f2:
            count += 1
    return int(count/4)
     

# GLOBALS
ORIGINAL_READCOUNT = 1000000
max_read_lengths = []
max_coverages = []
comp_times = [] 
overlap_counts = []
iteration = 0
outputdir = "."
transitive_edges = 0


def main():
    print "Pipeline_stage_a.py"
    
    parser = ArgumentParser(description=usage)
    parser.add_argument('--error_correction', dest='error_correction', action='store_true')
    parser.add_argument('--pre_merging', dest='pre_merging', type=int, default=0)
    parser.add_argument('--merging', dest='merging', action='store_true')
    parser.add_argument('--cliques', dest='cliques', action='store_true')
    parser.add_argument('--cliques_first_it', dest='cliques_first_it', action='store_true')
    parser.add_argument('--contig_merging', dest='contig_merge', action='store_true')
    parser.add_argument('--min_overlap_perc', dest='min_overlap_perc', type=int)
    parser.add_argument('--edge_threshold', dest='edge_threshold', type=float)
    parser.add_argument('--merging_threshold', dest='merging_threshold', type=float, default=0.99)
    parser.add_argument('--min_clique_size', dest='min_clique_size', type=int, default=5)
    parser.add_argument('--fno', dest='fno', type=int, default=2)
    parser.add_argument('--continue', dest='cntinue', action='store_true')
    parser.add_argument('--remove_branches', dest='remove_branches', action='store_true')
    parser.add_argument('--transitive_edges', dest='transitive_edges', type=int)
#    parser.add_argument('--output', dest='outputdir', type=str, default='.')
    parser.add_argument('--fastq', dest='fastq', type=str, required=True)
    parser.add_argument('--overlaps', dest='overlaps', type=str, required=True)
    args = parser.parse_args()

    if not (args.merging or args.cliques or args.error_correction or args.cliques_first_it):
        print "Nothing to do; specify at least one option."
        parser.print_help()
        
    if not (args.transitive_edges in [0, 1, 2, 3]):
        print "Transitive edges needs to be 0 (default, keep all edges), 1 (remove transitive edges), 2 (remov double transitive edges), or 3 (remove triple transitive edges)."
        parser.print_help()
    
    global iteration, max_read_lengths, max_coverages_comp_times, overlap_counts, ORIGINAL_READCOUNT, read_counts, outputdir, transitive_edges
#    ORIGINAL_READCOUNT = get_original_readcount(args.fastq)
    read_counts = [ORIGINAL_READCOUNT]
    original_overlaps = analyze_overlaps(args.overlaps)
    overlap_counts = [original_overlaps]
#   outputdir = args.outputdir
    transitive_edges = args.transitive_edges
    if not args.merging and args.cliques:
        no_inclusions = "true"
    else:
        no_inclusions = "false"

    # create a global log file; after every iteration the log file is appended to this global log file
    subprocess.call(["rm", "pipeline.log"])
    subprocess.call(["touch", "pipeline.log"])
    # remove existing stats file
    subprocess.call(["rm", "stats.txt"])
    subprocess.call(["touch", "stats.txt"])

    # first iteration    
    if not args.cntinue:
        if args.pre_merging or (args.merging and not (args.cliques or args.error_correction or args.cliques_first_it)):
            run_first_it_merge(args.fastq, args.overlaps, args.merging_threshold, 100)
            pm_count = 1
            while pm_count < args.pre_merging:
#                while overlap_counts[-1] > 0 and read_counts[-1] != read_counts[-2]:
                run_merging_it(args.merging_threshold, 100, 1, "false", "false")
                pm_count += 1
        elif args.error_correction:
            run_first_it_cliques(args.fastq, args.overlaps, args.edge_threshold, args.min_overlap_perc, args.min_clique_size, args.fno, "true", no_inclusions)
        elif args.cliques_first_it or args.cliques:
            run_first_it_cliques(args.fastq, args.overlaps, args.edge_threshold, args.min_overlap_perc, args.min_clique_size, args.fno, "false", no_inclusions)
        else:
            print "No first iteration... exiting"
            return -1
    
    # error correction if specified and not performed at first iteration
    if (args.pre_merging or args.cntinue) and args.error_correction:
        run_clique_it(args.edge_threshold, args.min_overlap_perc, args.min_clique_size, args.fno, "true", no_inclusions, 1000)
    elif (args.pre_merging or args.cntinue) and args.cliques_first_it:
        run_clique_it(args.edge_threshold, args.min_overlap_perc, args.min_clique_size, args.fno, "false", no_inclusions, 1000)
    
    # remaining iterations (until nothing left to be done)
    if not (args.cliques or args.merging):
        return 0
     
    while overlap_counts[-1] > 0 and read_counts[-1] != read_counts[-2]:
        if args.merging:
            while overlap_counts[-1] > 0 and read_counts[-1] != read_counts[-2]:
                if args.cliques:
                    run_merging_it(args.edge_threshold, 100, 3, "false", "false")
                else:
                    run_merging_it(args.merging_threshold, 0, 3, "false", "false")
            if args.cliques:
                run_clique_it(args.edge_threshold, args.min_overlap_perc, 2, args.fno, "false", no_inclusions, 0)
        elif args.cliques:
            run_clique_it(args.edge_threshold, args.min_overlap_perc, args.min_clique_size, args.fno, "false", no_inclusions, 0)

#    merge_contigs()
    if args.contig_merge:
        while overlap_counts[-1] > 0 and read_counts[-1] != read_counts[-2]:
            merge_contigs()

    print "Algorithm done; in total %d iterations" %iteration
    print "Maximum read length per iteration: \t", max_read_lengths
    print "Maximum # subreads per iteration: \t", max_coverages
    print "Number of input reads per iteration: \t", read_counts
    print "Number of overlaps found per iteration: \t", overlap_counts

def merge_contigs():
    global iteration, max_read_lengths, max_coverages_comp_times, read_counts, overlap_counts, outputdir, transitive_edges
    iteration += 1
    print "\n*************************************"
    print "**** Iteration %d = merge contigs ****" %iteration
    print "*************************************"
    min_clique_size = 2
    error_correction = "false"
    FNO = 3
    subprocess.check_call([viralquasispecies,
        "--singles", "singles.fastq",
        "--paired1", "paired1.fastq",
        "--paired2", "paired2.fastq",
        "--overlaps", "overlaps.txt",
        "--threads=%d" %THREADS,
        "--first_it=false",
        "--cliques=true",
        "--min_clique_size=%d" %min_clique_size,
        "--FNO=%d" %FNO,
        "--original_readcount=%d" %ORIGINAL_READCOUNT,
        "--error_correction=%s" %error_correction,
        "--remove_trans=%d" %transitive_edges,
#        "--output=%s" %outputdir
        "--optimize=false",
        "--base_path=%s" % selfpath
    ])
    copy_files(iteration, True) # copy pipeline log only
    [readcount, n_overlaps] = analyze_results()
    read_counts.append(readcount)
    overlap_counts.append(n_overlaps)
    print "***"


def run_first_it_merge(fastq, overlaps, edge_threshold, perc):
    global iteration, max_read_lengths, max_coverages_comp_times, read_counts, overlap_counts, outputdir, transitive_edges
    iteration += 1    
    print "\n**************************************"
    print "**** Iteration %d = first_it_merge ****" %iteration
    print "**************************************"
    min_clique_size = 2
    min_overlap_perc = perc
    error_correction = "false"
    FNO = 1
    subprocess.check_call([viralquasispecies, 
        "--singles", "%s/singles.fastq" %fastq, 
        "--paired1", "%s/paired1.fastq" %fastq, 
        "--paired2", "%s/paired2.fastq" %fastq, 
        "--overlaps=%s" %overlaps, 
        "--IDs", "../../pear_id_correspondance.csv", 
        "--threads=%d" %THREADS, 
        "--edge_threshold=%f" %edge_threshold, 
        "--first_it=true", 
        "--min_clique_size=%d" %min_clique_size, 
        "--min_overlap_perc=%d" %min_overlap_perc,
        "--FNO=%d" %FNO,
        "--original_readcount=%d" %ORIGINAL_READCOUNT,
        "--error_correction=%s" %error_correction,
        "--remove_trans=%d" %transitive_edges,
#        "--output=%s" %outputdir
        "--optimize=false",
        "--base_path=%s" % selfpath
    ])
    copy_files(iteration, True)
    [readcount, n_overlaps] = analyze_results()
    read_counts.append(readcount) 
    overlap_counts.append(n_overlaps)
    print "***"


def run_first_it_cliques(fastq, overlaps, edge_threshold, min_overlap_perc, min_clique_size, FNO, error_correction, no_inclusions):
    global iteration, max_read_lengths, max_coverages_comp_times, read_counts, overlap_counts, outputdir, transitive_edges
    iteration += 1
    print "\n****************************************"
    print "**** Iteration %d = first_it_cliques ****" %iteration
    print "****************************************"
    keep_singletons = 1000
    if overlap_counts[-1]/read_counts[-1] < 5:
        perc = min_overlap_perc
    else:
        perc = min_overlap_perc
    subprocess.check_call([viralquasispecies, 
        "--singles", "%s/singles.fastq" %fastq, 
        "--paired1", "%s/paired1.fastq" %fastq, 
        "--paired2", "%s/paired2.fastq" %fastq, 
        "--overlaps=%s" %overlaps, 
#        "--IDs=/home/jasmijn/5VM_experiments/id_correspondance_pear_subset_10_percent.csv",
#        "--IDs", "../../../pear_id_correspondance.csv",
#        "--IDs", "../pear_id_correspondance.csv", 
        "--threads=%d" %THREADS, 
        "--edge_threshold=%f" %edge_threshold, 
        "--first_it=true", 
        "--min_clique_size=%d" %min_clique_size, 
        "--min_overlap_perc=%d" %perc,
        "--FNO=%d" %FNO,
        "--original_readcount=%d" %ORIGINAL_READCOUNT,
        "--error_correction=%s" %error_correction,
        "--keep_singletons=%d" %keep_singletons,
        "--cliques=true",
        "--remove_trans=%d" %transitive_edges,
#        "--output=%s" %outputdir
        "--no_inclusion_overlaps=%s" %no_inclusions,
        "--optimize=true",
        "--base_path=%s" % selfpath
    ])
    copy_files(iteration, True)
    [readcount, n_overlaps] = analyze_results()
    read_counts.append(readcount) 
    overlap_counts.append(n_overlaps)
    print "***"


def run_merging_it(edge_threshold, min_overlap_perc, FNO, optimize, cliques):
    global iteration, max_read_lengths, max_coverages_comp_times, read_counts, overlap_counts, outputdir, transitive_edges
    iteration += 1
    print "\n*******************************"
    print "**** Iteration %d = merging ****" %iteration
    print "*******************************"
    min_clique_size = 2
    error_correction = "false"
    subprocess.check_call([viralquasispecies, 
        "--singles=%s" %"singles.fastq", 
        "--paired1", "paired1.fastq", 
        "--paired2", "paired2.fastq", 
        "--overlaps", "overlaps.txt", 
        "--threads=%d" %THREADS, 
        "--edge_threshold=%f" %edge_threshold, 
        "--first_it=false", 
        "--min_clique_size=%d" %min_clique_size, 
        "--min_overlap_perc=%d" %min_overlap_perc,
        "--FNO=%d" %FNO,
        "--original_readcount=%d" %ORIGINAL_READCOUNT,
        "--error_correction=%s" %error_correction,
#        "--output=%s" %outputdir
        "--optimize=%s" %optimize,
        "--cliques=%s" %cliques,
        "--remove_trans=%d" %transitive_edges,
        "--remove_multi_occ=true",
        "--base_path=%s" % selfpath
    ])
    copy_files(iteration, True) # copy pipeline log only
    [readcount, n_overlaps] = analyze_results()
    read_counts.append(readcount) 
    overlap_counts.append(n_overlaps)
    print "***"


def run_clique_it(edge_threshold, min_overlap_perc, min_clique_size, FNO, error_correction, no_inclusions, keep_singletons):
    global iteration, max_read_lengths, max_coverages_comp_times, read_counts, overlap_counts, outputdir, transitive_edges
    iteration += 1
    print "\n*******************************"
    print "**** Iteration %d = cliques ****" %iteration
    print "*******************************"
    if overlap_counts[-1]/read_counts[-1] < 5:
        perc = 0
    else:
        perc = min_overlap_perc
    subprocess.check_call([viralquasispecies, 
        "--singles=%s" %"singles.fastq", 
        "--paired1", "paired1.fastq", 
        "--paired2", "paired2.fastq", 
        "--overlaps", "overlaps.txt", 
        "--threads=%d" %THREADS, 
        "--edge_threshold=%f" %edge_threshold, 
        "--first_it=false", 
        "--min_clique_size=%d" %min_clique_size, 
        "--min_overlap_perc=%d" %perc,
        "--FNO=%d" %FNO,
        "--original_readcount=%d" %ORIGINAL_READCOUNT,
        "--error_correction=%s" %error_correction,
        "--cliques=true",
#        "--output=%s" %outputdir
        "--no_inclusion_overlaps=%s" %no_inclusions,
        "--remove_trans=%d" %transitive_edges,
        "--keep_singletons=%d" %keep_singletons,
        "--optimize=false",
        "--base_path=%s" % selfpath
    ])
    copy_files(iteration, True)
    [readcount, n_overlaps] = analyze_results()
    read_counts.append(readcount) 
    overlap_counts.append(n_overlaps)
    print "***"
    

def copy_files(it, pipeline_only=False):
    if not pipeline_only:
        subprocess.call(["cp", "singles.fastq", "it%d_singles.fastq" %it])
        subprocess.call(["cp", "paired1.fastq", "it%d_paired1.fastq" %it])
        subprocess.call(["cp", "paired2.fastq", "it%d_paired2.fastq" %it])
        subprocess.call(["cp", "overlaps.txt", "it%d_overlaps.txt" %it])
        subprocess.call(["cp", "subreads.txt", "it%d_subreads.txt" %it])
        subprocess.call(["cp", "graph.gfa", "it%d_graph.gfa" %it])
    subprocess.call("cat viralquasispecies.log >> pipeline.log", shell=True)


def analyze_results(cliques=False):
    max_cov = analyze_coverage()
    max_coverages.append(max_cov)
    [readcount, max_len] = analyze_fastq()
    max_read_lengths.append(max_len)
    if cliques:
        analyze_cliques()
    n_overlaps = analyze_overlaps('overlaps.txt')  
    return [readcount, n_overlaps]
    

def analyze_coverage():
    cov_counts = [0 for i in xrange(10000)]
    max_cov = 0
          
    infile = "subreads.txt"
    with open(infile, 'r') as f:
        for line in f:
            reads = line.split()
            cov = len(reads)-1
            cov_counts[cov-1] += 1
            if cov > max_cov:
                max_cov = cov
                
    print "cov\tcount (top 5)"
    for i in range(max_cov-5, max_cov):
        if cov_counts[i] != 0:
            print "%d\t%d" % (i+1, cov_counts[i])
    print "\n"    
    return max_cov
    

def analyze_fastq():
    len_counts = [0 for i in xrange(10000)]
    max_len = 0
          
    infile = "singles.fastq"
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
    print "longest read: ", max_len
    print "\n"    
    return [int(c/4), max_len]
    

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
    
    print "Overlaps:"                
    print "[-+, +-, ++, --]"
    print "p-p: ", pp_count
    print "p-s: ", ps_count
    print "s-p: ", sp_count
    print "s-s: ", ss_count

    total = sum(pp_count) + sum(ps_count) + sum(sp_count) + sum(ss_count)
    print total
    print "# lines: ", c
    print "\n"
    return total


if __name__ == '__main__':
    sys.exit(main())
