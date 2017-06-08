#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import os
import sys
import random
import subprocess
import shutil
from time import clock
from sets import Set

# ------------------------------

__author__ = "Jasmijn Baaijens"
__license__ = "GPL"

version = "0.3.0"
releasedate = "March 6, 2017"

usage = """
Program: SAVAGE - Strain Aware VirAl GEnome assembly
Version: %s
Release date: %s
Contact: Jasmijn Baaijens - j.a.baaijens@cwi.nl

SAVAGE assembles individual (viral) haplotypes from NGS data. It expects as
input single- and/or paired-end Illumina sequencing reads. Please note that the
paired-end reads are expected to be in forward-forward format, as output by
PEAR.

Run savage -h for a complete description of required and optional arguments.

For more information, please visit https://bitbucket.org/jbaaijens/savage
""" % (version, releasedate)

# ------------------------------

def main():
    parser = ArgumentParser(description=usage, formatter_class=RawTextHelpFormatter)
    basic = parser.add_argument_group('basic arguments')
    basic.add_argument('-s', dest='input_s', type=str, help='path to input fastq containing single-end reads')
    basic.add_argument('-p1', dest='input_p1', type=str, help='path to input fastq containing paired-end reads (/1)')
    basic.add_argument('-p2', dest='input_p2', type=str, help='path to input fastq containing paired-end reads (/2)')
    basic.add_argument('-m', '--min_overlap_len', dest='min_overlap_len', type=int, help='minimum overlap length required between reads')
    basic.add_argument('-t', '--num_threads', dest='threads', type=int, default=1, help='allowed number of cores')
    basic.add_argument('--split', dest='split_num', type=int, required=True, help='split the data set into patches s.t. 500 < coverage/split_num < 1000')
    basic.add_argument('--revcomp', dest='revcomp', action='store_true', help='use this option when paired-end input reads are in forward-reverse orientation;\nthis option will take reverse complements of /2 reads (specified with -p2)\nplease see the SAVAGE manual for more information about input read orientations')
#    basic.add_argument('--config', dest='config', type=str, help='path to config file containing parameter settings; for an example, \nplease see the SAVAGE repository on https://bitbucket.org/jbaaijens/savage')
    ref_guided = parser.add_argument_group('reference-guided mode')
    ref_guided.add_argument('--ref', dest='reference', type=str, help='reference genome in fasta format')
#    ref_guided.add_argument('--singles', dest='singles', type=str, help='single-end read alignments in SAM format')
#    ref_guided.add_argument('--paired', dest='paired', type=str, help='paired-end read alignments in SAM format')
    advanced = parser.add_argument_group('advanced arguments')
    advanced.add_argument('--no_stage_a', dest='stage_a', action='store_false', help='skip Stage a (initial contig formation)')
    advanced.add_argument('--no_stage_b', dest='stage_b', action='store_false', help='skip Stage b (extending initial contigs)')
    advanced.add_argument('--no_stage_c', dest='stage_c', action='store_false', help='skip Stage c (merging maximized contigs into master strains)')
    advanced.add_argument('--no_overlaps', dest='compute_overlaps', action='store_false', help='skip overlap computations (use existing overlaps file instead)')
    advanced.add_argument('--no_preprocessing', dest='preprocessing', action='store_false', help='skip preprocessing procedure (i.e. creating data patches)')
#    advanced.add_argument('--overlaps', dest='overlaps', type=str, help='skip overlap computations by using given overlaps file; please make sure \nto enter the full path!')
#    advanced.add_argument('--contigs', dest='contigs', type=str, help='contigs fastq file resulting from Stage a; \n--> use this option together with --no_stage_a')
    advanced.add_argument('--ignore_subreads', dest='use_subreads', action='store_false', help='ignore subread info from previous stage')
    advanced.add_argument('--merge_contigs', dest='merge_contigs', type=float, default=0.0, help='specify maximal distance between contigs for merging into master strains (stage c)')
    advanced.add_argument('--min_clique_size', dest='min_clique_size', type=int, default=4, help='minimum clique size used during error correction')
    advanced.add_argument('--overlap_len_stage_c', dest='overlap_stage_c', type=int, default=100, help='min_overlap_len used in stage c')
    advanced.add_argument('--contig_len_stage_c', dest='contig_len_stage_c', type=int, default=100, help='minimum contig length required for stage c input contigs')
    advanced.add_argument('--keep_branches', dest='remove_branches', action='store_false', help='disable merging along branches by removing them from the graph (stage b & c)')
    advanced.add_argument('--sfo_mm', dest='sfo_mm', type=int, default=50, help='input parameter -e=SFO_MM for sfo: maximal mismatch rate 1/SFO_MM')
    advanced.add_argument('--diploid', dest='diploid', action='store_true', help='use this option for diploid genome assembly')
    advanced.add_argument('--diploid_contig_len', dest='diploid_contig_len', type=int, default=200, help='minimum contig length required for diploid step contigs')
    advanced.add_argument('--diploid_overlap_len', dest='diploid_overlap_len', type=int, default=30, help='min_overlap_len used in diploid assembly step')
    advanced.add_argument('--average_read_len', dest='average_read_len', type=float, help='average length of the input reads; will be computed from the input if not specified')

    # store the path to the SAVAGE root directory
    base_path = os.path.dirname(os.path.abspath(__file__))
    # test if the SAVAGE binary functions properly
    binary_help = base_path + '/bin/ViralQuasispecies --help'
    try:
        subprocess.check_output(binary_help, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.output.decode() + '\n')
        sys.stderr.flush()
        sys.exit(1)

    if len(sys.argv[1:])==0:
        print usage
#        parser.print_help()
#        parser.print_usage()
        parser.exit()
    args = parser.parse_args()

    print """
-------------------------------------------
SAVAGE - Strain Aware VirAl GEnome assembly
-------------------------------------------
Version: %s
Author: %s
    """ % (version, __author__)
    print "Command used:"
    print ' '.join(sys.argv)
    print
    print "Parameter values:"
    for arg in vars(args):
        print arg, "=", getattr(args, arg)
    print

    FNULL = open(os.devnull, 'w')
    remove_branches = 'true' if args.remove_branches else 'false'
#    remove_branches = 'true'

    if not (args.stage_a or args.stage_b or args.stage_c or args.preprocessing or args.compute_overlaps):
        sys.stderr.write("Nothing to be done; please specify at least one task to perform.\n")
        sys.stderr.flush()
        sys.exit(1)

    if args.stage_a and args.stage_c and not args.stage_b:
        sys.stderr.write("""Options specified suggest running stages a and c, but skipping stage b.
                 If you really want to do this, then run stage a and c separately.\n""")
        sys.stderr.flush()
        sys.exit(1)

    if not args.preprocessing:
        preprocessing = False
    elif (args.stage_b or args.stage_c) and not args.stage_a:
        preprocessing = False
    else:
        preprocessing = True

#    if preprocessing:
    if not (args.input_s or (args.input_p1 and args.input_p2)):
        sys.stderr.write("""Please enter input fastq file(s) with -s and/or -p1,-p2.\n""")
        sys.stderr.flush()
        sys.exit(1)
    elif (args.input_p1 and not args.input_p2) or (args.input_p2 and not args.input_p1):
        sys.stderr.write("""For paired-end reads, please enter the fastq file(s) separately using -p1 and -p2.\n""")
        sys.stderr.flush()
        sys.exit(1)

    if args.reference:
        denovo = False
        if not os.path.exists(args.reference):
            sys.stderr.write("""ERROR: Reference fasta not found: %s \nPlease enter full path to file.\n""" % args.reference)
            sys.stderr.flush()
            sys.exit(1)
        subprocess.check_call("bwa index %s 1>/dev/null 2>&1" % args.reference, shell=True)
    else:
        denovo = True

    # analyze single-end input reads
    if args.input_s:
        [s_seq_count, s_total_len, s_longest_seq] = analyze_fastq(args.input_s)
    else:
        s_seq_count = 0
        s_total_len = 0
        s_longest_seq = 0
    # analyze paired-end input reads
    if args.input_p1:
        [p1_seq_count, p1_total_len, p1_longest_seq] = analyze_fastq(args.input_p1)
        [p2_seq_count, p2_total_len, p2_longest_seq] = analyze_fastq(args.input_p2)
        p_total_len = p1_total_len + p2_total_len
        p_seq_count = p1_seq_count + p2_seq_count
        if p1_seq_count != p2_seq_count:
            print "ERROR: Unequal number of /1 and /2 reads. Exiting."
            sys.exit(1)
    else:
        p_seq_count = 0
        p_total_len = 0
        p_longest_seq = 0
    total_seq_len = s_total_len + p_total_len
    total_seq_count = s_seq_count + p_seq_count
    if not total_seq_len > 0:
        print "ERROR: Total input length is zero. Exiting."
        sys.exit(1)

    # compute average read length
    if args.average_read_len > 0:
        average_read_len = args.average_read_len
    else:
        average_read_len = total_seq_len/total_seq_count

    # print stats
    print "Input fastq stats:"
    print "Number of single-end reads =", s_seq_count
    print "Number of paired-end reads =", p_seq_count
    print "Total number of bases =", total_seq_len
    print "Average read length =", average_read_len

    max_tip_len = int(round(average_read_len)) # average input read length

    if not args.min_overlap_len:
        m = 0.6*average_read_len # 60% of average input read length
        min_overlap_len = int(round(m))
        assert min_overlap_len > 0
    else:
        min_overlap_len = args.min_overlap_len

    if min_overlap_len < 100:
        print "----------------------------------------------------------------"
        print "WARNING: min_overlap_len = %s" % min_overlap_len
        print "For more accurate error correction, increase the minimal overlap length using --min_overlap_len"
        print "----------------------------------------------------------------"

    # Preprocessing: split data into patches
    if preprocessing:
        print "*******************"
        print "Preprocessing input"
        overwrite_dir('stage_a')
        # split fastq into patches
        if args.input_s:
            subprocess.check_call("%s/scripts/random_split_fastq.py --input %s --output stage_a/singles --split_num %s" % (base_path, args.input_s, args.split_num), shell=True)
        if args.input_p1 and args.input_p2:
            subprocess.check_call("%s/scripts/random_split_fastq.py --input %s --input2 %s --output stage_a/paired --split_num %s" % (base_path, args.input_p1, args.input_p2, args.split_num), shell=True)
        for patch_num in range(args.split_num):
            print "\rpatch %d" % patch_num,
            sys.stdout.flush()
            # create a separate directory for every patch
            overwrite_dir('stage_a/patch%d/input_fas' % patch_num)
            # rename or create single-end reads file
            if args.input_s:
                singles_count = int(file_len('stage_a/singles.%d.fastq' % patch_num)/4)
                subprocess.check_call("%s/scripts/rename_fas.py --in stage_a/singles.%d.fastq --out stage_a/patch%d/input_fas/singles.fastq" % (base_path, patch_num, patch_num), shell=True)
            else:
                singles_count = 0
                subprocess.check_call("touch stage_a/patch%d/input_fas/singles.fastq" % patch_num, shell=True)
            # rename or create paired-end reads files
            if args.input_p1 and args.input_p2:
                subprocess.check_call("%s/scripts/rename_fas.py --in stage_a/paired1.%d.fastq --out stage_a/patch%d/input_fas/paired1.fastq --id_start %d" % (base_path, patch_num, patch_num, singles_count), shell=True)
                if args.revcomp:
                    subprocess.check_call("%s/scripts/rename_fas.py --revcomp --in stage_a/paired2.%d.fastq --out stage_a/patch%d/input_fas/paired2.fastq --id_start %d" % (base_path, patch_num, patch_num, singles_count), shell=True)
                else:
                    subprocess.check_call("%s/scripts/rename_fas.py --in stage_a/paired2.%d.fastq --out stage_a/patch%d/input_fas/paired2.fastq --id_start %d" % (base_path, patch_num, patch_num, singles_count), shell=True)
            else:
                subprocess.check_call("touch stage_a/patch%d/input_fas/paired1.fastq stage_a/patch%d/input_fas/paired2.fastq" % (patch_num, patch_num), shell=True)
            if not denovo:
                # Check for reference fasta
                if not os.path.exists(args.reference):
                    sys.stderr.write("""ERROR: Reference fasta not found: %s \nPlease enter full path to file.\n""" % args.reference)
                    sys.stderr.flush()
                    sys.exit(1)
                # Run BWA to get alignments
                fastq_path = "stage_a/patch%d/input_fas" % patch_num
                subprocess.check_call("bwa mem %s %s/singles.fastq 1> %s/singles.sam 2> /dev/null" % (args.reference, fastq_path, fastq_path), shell=True)
                subprocess.check_call("bwa mem %s %s/paired1.fastq %s/paired2.fastq 1> %s/paired.sam 2> /dev/null" % (args.reference, fastq_path, fastq_path, fastq_path), shell=True)
        # move original reads to separate directory
        overwrite_dir('stage_a/original_reads')
        subprocess.check_call("mv stage_a/*.*.fastq stage_a/original_reads/", shell=True)
        print "\rDone!" + ' ' * 40
        sys.stdout.flush()

    # For every patch, compute all suffix-prefix overlaps
    overlaps = "../original_overlaps.txt"
    if args.compute_overlaps:
        print "********************"
        print "Overlap computations"
        for patch_num in range(args.split_num):
            print "\r" + " " * 60,
            print "\rpatch %d" % patch_num,
            sys.stdout.flush()
            os.chdir('stage_a/patch%d' % patch_num)
            # find all overlaps
            if denovo:
                preprocessing_denovo(args.min_overlap_len, args.sfo_mm, args.threads, base_path)
            else:
                paired = args.input_p1 and args.input_p2
                preprocessing_ref(args.min_overlap_len, args.reference, base_path, paired)
            os.chdir('../..')
        print "\rDone!" + " " * 60
        sys.stdout.flush()
    elif args.stage_a:
        # check if overlaps file already exists for every patch
        for patch_num in range(args.split_num):
            if not os.path.exists('stage_a/patch%d/original_overlaps.txt' % patch_num):
                sys.stderr.write("""Assuming existing overlaps file 'stage_a/patch%d/original_overlaps.txt'.
                         Please make sure this file exists, or run overlap computations first\n""" % patch_num)
                sys.stderr.flush()
                sys.exit(1)

    final_contig_file = ""

    # Run SAVAGE Stage a: error correction and initial contig formation
    if args.stage_a:
        print "**************"
        print "SAVAGE Stage a"
        sys.stdout.flush()
        os.chdir('stage_a')
        # process every patch separately
        for patch_num in range(args.split_num):
            os.chdir('patch%d' % patch_num)
            # create stage_a directory FOR THIS PATCH
            overwrite_dir('stage_a')
            os.chdir('stage_a')
            edge_threshold = 0.97
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --stage a --fastq ../input_fas --overlaps %s --min_overlap_len %d --num_threads %d --remove_branches %s --max_tip_len %s --edge_threshold %s --clique_size_EC %s" %(base_path, overlaps, args.min_overlap_len, args.threads, remove_branches, max_tip_len, edge_threshold, args.min_clique_size), shell=True)
            os.chdir('../..')
        os.chdir('..')
        # combine contigs from all patches
        remove_file('stage_a/combined_singles.fastq')
        remove_file('stage_a/subreads.txt')
        subprocess.check_call("%s/scripts/combine_contigs.py --split %s --paired_to_single" % (base_path, args.split_num), shell=True)
        # now rename the merged fastq and convert to fasta
        subprocess.check_call("%s/scripts/rename_fas.py --in stage_a/combined_singles.fastq --out stage_a/singles.fastq" % base_path, shell=True)
        subprocess.check_call("%s/scripts/fastq2fasta.py stage_a/singles.fastq contigs_stage_a.fasta" % base_path, shell=True)
        print "Done!"
        final_contig_file = "contigs_stage_a.fasta"
    # else:
    #     print "Stage a skipped"

    # Run SAVAGE Stage b: build maximized contigs
    if args.stage_b:
        print "**************"
        print "SAVAGE Stage b"
        # prepare input files
        overwrite_dir('stage_b')
        if not (os.path.exists('stage_a/singles.fastq') and os.path.exists('contigs_stage_a.fasta')):
            print """Contigs file from Stage a not found. Please make sure that both
                     'stage_a/singles.fastq' and 'contigs_stage_a.fasta' exist. If
                     absent, please rerun Stage a."""
            sys.exit(1)
        elif os.stat("contigs_stage_a.fasta").st_size == 0:
            print "Empty set of contigs from Stage a (contigs_stage_a.fasta) --> Exiting SAVAGE."
            sys.exit(1)
        subprocess.call(['cp', 'stage_a/singles.fastq', 'stage_b/singles.fastq'], stdout=FNULL, stderr=FNULL)
        pident = 98
        overlaps = run_blast('a', pident, base_path, args.min_overlap_len)
        sys.stdout.flush()
        # run SAVAGE
        os.chdir('stage_b')
        if args.use_subreads:
            subprocess.check_call("cp ../stage_a/subreads.txt subreads.txt", shell=True)
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --stage b --fastq ../stage_b --overlaps %s --use_subreads --min_overlap_len %d --num_threads %d --remove_branches %s --max_tip_len %s" % (base_path, overlaps, args.min_overlap_len, args.threads, remove_branches, max_tip_len), shell=True)
        else:
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --stage b --fastq ../stage_b --overlaps %s --min_overlap_len %d --num_threads %d --remove_branches %s --max_tip_len %s" % (base_path, overlaps, args.min_overlap_len, args.threads, remove_branches, max_tip_len), shell=True) # note: not using stage a subreads
        os.chdir('..')
        subprocess.check_call("%s/scripts/fastq2fasta.py stage_b/singles.fastq contigs_stage_b.fasta" % base_path, shell=True)
        print "Done!"
        final_contig_file = "contigs_stage_b.fasta"
    # else:
    #     print "Stage b skipped"

    # Run SAVAGE Stage c: build master strains
    if args.stage_c:
        print "**************"
        print "SAVAGE Stage c"
        # parameters
        if args.overlap_stage_c:
            min_overlap_len = args.overlap_stage_c
#        else:
#            min_overlap_len = args.min_overlap_len

        if args.contig_len_stage_c:
            min_contig_len = args.contig_len_stage_c
#        else:
#            min_contig_len = int(round(average_read_len)) # average input read length
        # prepare input files
        overwrite_dir('stage_c')
        if not (os.path.exists('stage_b/singles.fastq') and os.path.exists('contigs_stage_b.fasta')):
            print """Contigs file from Stage b not found. Please make sure that both
                     'stage_b/singles.fastq' and 'contigs_stage_b.fasta' exist. If
                     absent, please rerun Stage b."""
            sys.exit(1)
        elif os.stat("contigs_stage_b.fasta").st_size == 0:
            print "Empty set of contigs from Stage b (contigs_stage_b.fasta) --> Exiting SAVAGE"
            sys.exit(1)
        subprocess.call(['cp', 'stage_b/singles.fastq', 'stage_c/singles.fastq'], stdout=FNULL, stderr=FNULL)
        pident = 100*(0.99-args.merge_contigs)
        overlaps = run_blast('b', pident, base_path, min_overlap_len)
        sys.stdout.flush()
        # run SAVAGE
        os.chdir('stage_c')
        if args.use_subreads:
            subprocess.check_call("cp ../stage_b/subreads.txt subreads.txt", shell=True)
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --fastq ../stage_c --overlaps %s --merge_contigs %f --stage c --min_overlap_len %d --use_subreads --num_threads %d --remove_branches %s --min_read_len %d --max_tip_len %s" % (base_path, overlaps, args.merge_contigs, min_overlap_len, args.threads, remove_branches, min_contig_len, max_tip_len), shell=True)
        else:
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --fastq ../stage_c --overlaps %s --merge_contigs %f --stage c --min_overlap_len %d --num_threads %d --remove_branches %s --min_read_len %d --max_tip_len %s" % (base_path, overlaps, args.merge_contigs, min_overlap_len, args.threads, remove_branches, min_contig_len, max_tip_len), shell=True)
        os.chdir('..')
        subprocess.check_call("%s/scripts/fastq2fasta.py stage_c/singles.fastq contigs_stage_c.fasta" % base_path, shell=True)
        subprocess.call("rm blastout* contigs_db*", shell=True)
        print "Done!"
        final_contig_file = "contigs_stage_c.fasta"
    # else:
    #     print "Stage c skipped"
    if args.diploid:
        # final diploid contig merging
        print "**************"
        print "SAVAGE diploid"
        # parameters
        min_overlap_len = args.diploid_overlap_len
        min_contig_len = args.diploid_contig_len
        # prepare input files
        overwrite_dir('diploid')
        if not (os.path.exists('stage_c/singles.fastq') and os.path.exists('contigs_stage_c.fasta')):
            print """Contigs file from Stage c not found. Please make sure that both
                     'stage_c/singles.fastq' and 'contigs_stage_c.fasta' exist. If
                     absent, please rerun Stage c."""
            sys.exit(1)
        elif os.stat("contigs_stage_c.fasta").st_size == 0:
            print "Empty set of contigs from Stage c (contigs_stage_c.fasta) --> Exiting SAVAGE."
            sys.exit(1)
        subprocess.call(['cp', 'stage_c/singles.fastq', 'diploid/singles.fastq'], stdout=FNULL, stderr=FNULL)
        pident = 98
        overlaps = run_blast('c', pident, base_path, min_overlap_len)
        sys.stdout.flush()
        # run SAVAGE
        os.chdir('diploid')
        if args.use_subreads:
            subprocess.check_call("cp ../stage_c/subreads.txt subreads.txt", shell=True)
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --fastq ../stage_c --overlaps %s --merge_contigs %f --stage c --min_overlap_len %d --use_subreads --num_threads %d --remove_branches %s --min_read_len %d --diploid --max_tip_len %s" % (base_path, overlaps, args.merge_contigs, min_overlap_len, args.threads, remove_branches, min_contig_len, max_tip_len), shell=True)
        else:
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --fastq ../stage_c --overlaps %s --merge_contigs %f --stage c --min_overlap_len %d --num_threads %d --remove_branches %s --min_read_len %d --diploid --max_tip_len %s" % (base_path, overlaps, args.merge_contigs, min_overlap_len, args.threads, remove_branches, min_contig_len, max_tip_len), shell=True)
        os.chdir('..')
        subprocess.check_call("%s/scripts/fastq2fasta.py diploid/singles.fastq diploid_contigs.fasta" % base_path, shell=True)
        subprocess.call("rm blastout* contigs_db*", shell=True)
        print "Done!"
        final_contig_file = "diploid_contigs.fasta"

    print """**************
SAVAGE assembly has been completed, the final contig set was written to:

        %s

Optionally, you can now apply frequency estimation using freq_est.py. Please see
the manual page for more information. If you have a reference genome for your
sample, it can also be helpful to run the strain count estimator on the resulting
contigs (estimate_strain_count.py). More information about this can also be
found in the manual.

Thank you for using SAVAGE!
    """ % final_contig_file

    FNULL.close()
    return


# ------------------------------

def file_len(fname):
    with open(fname) as f:
        i = 0
        for i, l in enumerate(f):
            pass
    if i > 0:
        linecount = i + 1
    else:
        linecount = 0
    assert linecount % 4 == 0
    return linecount

def overwrite_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
    else:
        shutil.rmtree(dir) #removes all the subdirectories!
        os.makedirs(dir)
    return

def remove_file(filename):
    if os.path.exists(filename):
        os.remove(filename)
    return

def analyze_fastq(filename):
    total_len = 0
    longest_seq = 0
    i = 0
    with open(filename) as f:
        for line in f:
            if i%4 == 1: # sequence line in fastq
                l = len(line.strip('\n'))
                total_len += l
                if l > longest_seq:
                    longest_seq = l
            i += 1
        assert i % 4 == 0 # fastq
    seq_count = int(i/4)
    return [seq_count, total_len, longest_seq]

def preprocessing_denovo(min_overlap_len, sfo_mm, threads, base_path):
    print "- De novo overlap computations",
    sys.stdout.flush()
    # prepare fasta
    subprocess.check_call("cat input_fas/singles.fastq input_fas/paired1.fastq input_fas/paired2.fastq > tmp.fastq", shell=True)
    subprocess.check_call("%s/scripts/rename_fas.py --in tmp.fastq --out s_p1_p2.fastq" % base_path, shell=True)
    subprocess.check_call("%s/scripts/fastq2fasta.py s_p1_p2.fastq s_p1_p2.fasta" % base_path, shell=True)
    singles_count = int(file_len('input_fas/singles.fastq')/4.0)
    paired_count = int(file_len('input_fas/paired1.fastq')/4.0)
    assert paired_count == int(file_len('input_fas/paired2.fastq')/4.0)
    # run SFO
    print "- Running SFO",
    sys.stdout.flush()
    subprocess.check_call("%s/sfo_2011_5/builder s_p1_p2.fasta" % base_path, shell=True)
    if paired_count > 0:
        sfo_len = int(round(min_overlap_len / 2))
    else:
        sfo_len = min_overlap_len
    # run sfoverlap
    try:
        subprocess.check_call("%s/sfo_2011_5/sfoverlap --parallel %d --indels -e %d -t %d s_p1_p2.fasta > tmp_overlaps.out 2>/dev/null" % (base_path, threads, sfo_mm, sfo_len), shell=True)
        subprocess.check_call("%s/sfo_2011_5/maxoverlaps < tmp_overlaps.out > sfoverlaps.out" % (base_path), shell=True)
        subprocess.check_call("rm tmp_overlaps.out", shell=True)
    except subprocess.CalledProcessError as e:
        print "-> sfoverlap failed, running blast instead"
        subprocess.check_call("makeblastdb -in s_p1_p2.fasta -dbtype nucl -out s_p1_p2.db 1>/dev/null 2>&1", shell=True)
        subprocess.check_call("blastn -db s_p1_p2.db -query s_p1_p2.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen' -out blastout.tsv -perc_identity 98", shell=True)
        subprocess.check_call("%s/scripts/blast2sfo.py --in blastout.tsv --out sfoverlaps.out -m %s" % (base_path, sfo_len), shell=True)
        subprocess.check_call("rm blastout.tsv", shell=True)
    # run postprocessing scripts
    print "\b" * 12 + "Processing output",
    sys.stdout.flush()
    subprocess.check_call("%s/scripts/sfo2overlaps.py --in sfoverlaps.out --out original_overlaps.txt --num_singles %d --num_pairs %d 1> /dev/null" % (base_path, singles_count, paired_count), shell=True)
#    print "Overlaps are ready!"
    subprocess.check_call("rm tmp.fastq s_p1_p2.* sfoverlaps.out", shell=True)
    return

def preprocessing_ref(min_overlap_len, reference, base_path, paired):
    print "Reference-guided overlap computations",
    sys.stdout.flush()
    # Check for reference fasta
    if not os.path.exists(reference):
        sys.stderr.write("""ERROR: Reference fasta not found: %s \nPlease enter full path to file.\n""" % reference)
        sys.stderr.flush()
        sys.exit(1)
    # Induce overlaps from alignment
    if paired:
        paired_overlap_len = int(round(min_overlap_len / 2))
        subprocess.check_call("%s/scripts/sam2overlaps.py --sam_p input_fas/paired.sam --sam_s input_fas/singles.sam --ref %s --min_overlap_len %d --out original_overlaps.txt" %(base_path, reference, paired_overlap_len), shell=True)
    else:
        subprocess.check_call("%s/scripts/sam2overlaps.py --sam_s input_fas/singles.sam --ref %s --min_overlap_len %d --out original_overlaps.txt" %(base_path, reference, min_overlap_len), shell=True)
    # if singles and paired:
    #     subprocess.check_call("%s/scripts/sam2overlaps.py --sam_p %s --sam_s %s --ref %s --min_overlap_len %d --out original_overlaps.txt" %(base_path, paired, singles, reference, min_overlap_len/2), shell=True)
    # elif singles:
    #     subprocess.check_call("%s/scripts/sam2overlaps.py --sam_s %s --ref %s --min_overlap_len %d --out original_overlaps.txt" %(base_path, singles, reference, min_overlap_len), shell=True)
    # elif paired:
    #     subprocess.check_call("%s/scripts/sam2overlaps.py --sam_p %s --ref %s --min_overlap_len %d --out original_overlaps.txt" %(base_path, paired, reference, min_overlap_len/2), shell=True)
    return

def run_blast(previous_stage, pident, base_path, min_overlap_len):
    overlaps_file = "contig_overlaps.txt"
    subprocess.check_call("makeblastdb -in contigs_stage_%s.fasta -dbtype nucl -out contigs_db 1>/dev/null 2>&1" % (previous_stage), shell=True)
    subprocess.check_call("blastn -db contigs_db -query contigs_stage_%s.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen' -out blastout_contigs.tsv -perc_identity %s" % (previous_stage, pident), shell=True)
    subprocess.check_call("%s/scripts/blast2overlaps.py --in blastout_contigs.tsv --out %s --min_overlap_len %d" % (base_path, overlaps_file, min_overlap_len), shell=True)
    overlaps_path = "../" + overlaps_file
    return overlaps_path


if __name__ == '__main__':
    sys.exit(main())
