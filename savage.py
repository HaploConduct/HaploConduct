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
from collections import namedtuple

# ------------------------------

__author__ = "Jasmijn Baaijens"
__license__ = "GPL"

version = "0.4.0"
releasedate = "July 17, 2017"

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

For the complete manual, please visit https://bitbucket.org/jbaaijens/savage
""" % (version, releasedate)

# ------------------------------

InputStruct = namedtuple("InputStruct", "input_s input_p1 input_p2 read_len fragmentsize stddev")

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
    basic.add_argument('-o', '--outdir', dest='outdir', type=str, help='specify output directory')
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
    advanced.add_argument('--no_assembly', dest='no_assembly', action='store_true', help='skip all assembly steps; only use this option when using --count_strains separate from assembly (e.g. on a denovo assembly)')
    advanced.add_argument('--count_strains', dest='count_strains', action='store_true', help='compute a lower bound on the number of strains in this sample; note: this requires a reference genome.')
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
    advanced.add_argument('--no_filtering', dest='filtering', action='store_false', help='disable kallisto-based filtering of contigs')
    advanced.add_argument('--max_tip_len', dest='max_tip_len', type=int, help='maximum extension length for a sequence to be called a tip')

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

    # store start dir
    start_dir = os.getcwd()

    FNULL = open(os.devnull, 'w')
    remove_branches = 'true' if args.remove_branches else 'false'
#    remove_branches = 'true'

    if args.reference:
        denovo = False
        if not os.path.exists(args.reference):
            sys.stderr.write("""\nERROR: Reference fasta not found: %s \nPlease enter full path to file.\n""" % args.reference)
            sys.stderr.flush()
            sys.exit(1)
        if not os.path.exists(args.reference + ".bwt"):
            print "Building index for reference genome...",
            try:
                message = subprocess.check_output("bwa index %s 2>&1" % args.reference, shell=True)
            except subprocess.CalledProcessError as e:
                print "ERROR running bwa index:"
                print e.output
                sys.exit(1)
            print "done!\n"
    else:
        denovo = True

    if args.no_assembly:
        print "Skipping assembly because --no_assembly flag was used.\n"

    if (args.no_assembly or not (args.stage_a or args.stage_b or args.stage_c or args.preprocessing or args.compute_overlaps or args.diploid)):
        if args.count_strains:
            only_strain_count = True
        else:
            sys.stderr.write("Nothing to be done; please specify at least one task to perform.\n")
            sys.stderr.write("Make sure to use --no_assembly ONLY in combination with --count_strains.\n")
            sys.stderr.flush()
            sys.exit(1)
    else:
        only_strain_count = False

    if only_strain_count:
        if os.path.exists('diploid_contigs.fasta'):
            final_contig_file = 'diploid_contigs.fasta'
        elif os.path.exists('contigs_stage_c.fasta'):
            final_contig_file = 'contigs_stage_c.fasta'
        elif os.path.exists('contigs_stage_b.fasta'):
            final_contig_file = 'contigs_stage_b.fasta'
        elif os.path.exists('contigs_stage_a.fasta'):
            final_contig_file = 'contigs_stage_a.fasta'
        else:
            sys.stderr.write("No contigs found for estimating strain count. Please run savage assembly first.\n")
            sys.stderr.flush()
            sys.exit(1)

        if args.reference:
            # run estimate_strain_count.py
            run_strain_count(args.reference, final_contig_file, base_path)
            FNULL.close()
            return # No assembly, so we are done!
        else:
            sys.stderr.write("Please specify a reference genome when using --count_strains.\n")
            sys.stderr.flush()
            sys.exit(1)
    else:
        final_contig_file = ""

    if args.stage_a and args.stage_c and not args.stage_b:
        sys.stderr.write("""Options specified suggest running stages a and c, but skipping stage b.
                 If you really want to do this, then run stage a and c separately.\n""")
        sys.stderr.flush()
        sys.exit(1)

    if args.compute_overlaps and (args.stage_b or args.stage_c) and not args.stage_a:
        sys.stderr.write("""Options specified suggest computing overlaps and running stages b and/or c, but skipping stage a.
                 Please add --no_overlaps flag to skip overlap computations as well.\n""")
        sys.stderr.flush()
        sys.exit(1)

    if not args.preprocessing:
        preprocessing = False
    elif (args.stage_a or args.stage_b or args.stage_c or args.diploid) and not args.compute_overlaps:
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


    print "Parameter values:"
    for arg in vars(args):
        print arg, "=", getattr(args, arg)
    print

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
            sys.stderr.write("""ERROR: Unequal number of /1 and /2 reads. Exiting.\n""")
            sys.stderr.flush()
            sys.exit(1)
    else:
        p_seq_count = 0
        p_total_len = 0
        p_longest_seq = 0
    total_seq_len = s_total_len + p_total_len
    total_seq_count = s_seq_count + p_seq_count
    if not total_seq_len > 0:
        sys.stderr.write("""ERROR: Total input length is zero. Exiting.\n""")
        sys.stderr.flush()
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
    print "Average sequence length = %.1f" % average_read_len
    print
#    fragmentsize = args.fragmentsize if args.fragmentsize else average_read_len
#    stddev = args.stddev if args.stddev else 20
    fragmentsize = average_read_len
    stddev = 20
    input_info = InputStruct(args.input_s, args.input_p1, args.input_p2, average_read_len, fragmentsize, stddev)

    if args.max_tip_len is None:
        max_tip_len = int(round(average_read_len)) # average input read length
        print "Using max_tip_len =", max_tip_len
        assert max_tip_len > 0
    elif args.max_tip_len >= 0:
        max_tip_len = args.max_tip_len
    else:
        sys.stderr.write("""\nERROR: invalid --max_tip_len %s. Please make sure this value is non-negative. Exiting.\n\n""" % args.max_tip_len)
        sys.stderr.flush()
        sys.exit(1)

    if args.min_overlap_len is None:
        m = 0.6*average_read_len # 60% of average input read length
        min_overlap_len = int(round(m))
        print "Using min_overlap_len =", min_overlap_len
        print
        assert min_overlap_len > 0
    else:
        min_overlap_len = args.min_overlap_len

    if min_overlap_len < 80:
        print "----------------------------------------------------------------"
        print "WARNING: min_overlap_len = %s" % min_overlap_len
        print "For more accurate error correction, increase the minimal overlap length using --min_overlap_len"
        print "----------------------------------------------------------------"

    # Preprocessing: split data into patches
    if preprocessing:
        print "*******************"
        print "Preprocessing input"
        overwrite_dir(args.outdir + '/stage_a')
        # split fastq into patches
        if args.input_s:
            subprocess.check_call("%s/scripts/random_split_fastq.py --input %s --output %s/stage_a/singles --split_num %s" % (base_path, args.input_s, args.outdir, args.split_num), shell=True)
        if args.input_p1 and args.input_p2:
            subprocess.check_call("%s/scripts/random_split_fastq.py --input %s --input2 %s --output %s/stage_a/paired --split_num %s" % (base_path, args.input_p1, args.input_p2, args.outdir, args.split_num), shell=True)
        # move to output directory
        if args.outdir not in [None, '.', './']:
            os.chdir(args.outdir)
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
                    sys.stderr.write("""\nERROR: Reference fasta not found: %s \nPlease enter full path to file.\n""" % args.reference)
                    sys.stderr.flush()
                    sys.exit(1)
                # Run BWA to get alignments
                fastq_path = "stage_a/patch%d/input_fas" % patch_num
                try:
                    message1 = subprocess.check_call("bwa mem %s %s/singles.fastq 1> %s/singles.sam 2> /dev/null" % (args.reference, fastq_path, fastq_path), shell=True)
                    message2 = subprocess.check_call("bwa mem %s %s/paired1.fastq %s/paired2.fastq 1> %s/paired.sam 2> /dev/null" % (args.reference, fastq_path, fastq_path, fastq_path), shell=True)
                except subprocess.CalledProcessError as e:
                    print "ERROR running bwa mem:"
                    print e.output
                    sys.exit(1)
                    # subprocess.check_call("bwa index %s 1>/dev/null 2>&1" % args.reference, shell=True)
                    # subprocess.check_call("bwa mem %s %s/singles.fastq 1> %s/singles.sam 2> /dev/null" % (args.reference, fastq_path, fastq_path), shell=True)
                    # subprocess.check_call("bwa mem %s %s/paired1.fastq %s/paired2.fastq 1> %s/paired.sam 2> /dev/null" % (args.reference, fastq_path, fastq_path, fastq_path), shell=True)

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
                preprocessing_denovo(min_overlap_len, args.sfo_mm, args.threads, base_path)
            else:
                paired = args.input_p1 and args.input_p2
                preprocessing_ref(min_overlap_len, args.reference, base_path, paired)
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
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --stage a --fastq ../input_fas --overlaps %s --min_overlap_len %d --num_threads %d --remove_branches %s --max_tip_len %s --edge_threshold %s --clique_size_EC %s" %(base_path, overlaps, min_overlap_len, args.threads, remove_branches, max_tip_len, edge_threshold, args.min_clique_size), shell=True)
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
    elif not (args.stage_b or args.stage_c or args.diploid or args.count_strains):
        if os.path.exists('contigs_stage_a.fasta'):
            os.remove('contigs_stage_a.fasta')

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
        pident = 99
        # find contig overlaps
        try:
            sfo_mm = 200 # 0.5% errors
            paired_count = 0
            singles_count = int(file_len('contigs_stage_a.fasta')/2)
            overlaps = run_sfo('a', sfo_mm, base_path, min_overlap_len, args.threads, singles_count, paired_count)
        except subprocess.CalledProcessError as e:
            print "\nRUNNING BLAST....\n"
            overlaps = run_blast('a', pident, base_path, min_overlap_len)
            subprocess.call("rm blastout* contigs_db*", shell=True)
        sys.stdout.flush()
        # run SAVAGE
        os.chdir('stage_b')
        subprocess.check_call(["cp", overlaps, "original_overlaps.txt"])
        if args.use_subreads:
            subprocess.check_call("cp ../stage_a/subreads.txt subreads.txt", shell=True)
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --stage b --fastq ../stage_b --overlaps %s --use_subreads --min_overlap_len %d --num_threads %d --remove_branches %s --max_tip_len %s" % (base_path, overlaps, min_overlap_len, args.threads, remove_branches, max_tip_len), shell=True)
        else:
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --stage b --fastq ../stage_b --overlaps %s --min_overlap_len %d --num_threads %d --remove_branches %s --max_tip_len %s" % (base_path, overlaps, min_overlap_len, args.threads, remove_branches, max_tip_len), shell=True) # note: not using stage a subreads
        os.remove(overlaps)
        os.chdir('..')
        subprocess.check_call("%s/scripts/fastq2fasta.py stage_b/singles.fastq contigs_stage_b.fasta" % base_path, shell=True)
        print "Done!"
        final_contig_file = "contigs_stage_b.fasta"
        # apply frequency-based filtering
        if args.filtering:
            try:
                freq_filtering("contigs_stage_b.fasta", "stage_b/singles.fastq", 0, input_info)
            except subprocess.CalledProcessError as e:
                print "\nKallisto not found - skipping this filtering step.\n"
    elif final_contig_file != "":
        if os.path.exists('stage_b'):
            shutil.rmtree('stage_b') #removes all the subdirectories!
        if os.path.exists('contigs_stage_b.fasta'):
            os.remove('contigs_stage_b.fasta')

    # Run SAVAGE Stage c: build master strains
    if args.stage_c:
        print "**************"
        print "SAVAGE Stage c"
        # parameters
        if args.overlap_stage_c:
            min_overlap_len = args.overlap_stage_c
#        else:
#            min_overlap_len = args.min_overlap_len

        min_contig_len = args.contig_len_stage_c
#        if args.contig_len_stage_c:
#            min_contig_len = args.contig_len_stage_c
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
        # find contig overlaps
        try:
            sfo_mm = 1 + (0.99 - args.merge_contigs) / (args.merge_contigs + 0.01)
            paired_count = 0
            singles_count = int(file_len('contigs_stage_b.fasta')/2)
            overlaps = run_sfo('b', sfo_mm, base_path, min_overlap_len, args.threads, singles_count, paired_count)
        except subprocess.CalledProcessError as e:
            print "\nRUNNING BLAST....\n"
            overlaps = run_blast('b', pident, base_path, min_overlap_len)
            subprocess.call("rm blastout* contigs_db*", shell=True)
        sys.stdout.flush()
        # run SAVAGE
        os.chdir('stage_c')
        subprocess.check_call(["cp", overlaps, "original_overlaps.txt"])
        if args.use_subreads:
            subprocess.check_call("cp ../stage_b/subreads.txt subreads.txt", shell=True)
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --fastq ../stage_c --overlaps %s --merge_contigs %f --stage c --min_overlap_len %d --use_subreads --num_threads %d --remove_branches %s --min_read_len %d --max_tip_len %s" % (base_path, overlaps, args.merge_contigs, min_overlap_len, args.threads, remove_branches, min_contig_len, max_tip_len), shell=True)
        else:
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --fastq ../stage_c --overlaps %s --merge_contigs %f --stage c --min_overlap_len %d --num_threads %d --remove_branches %s --min_read_len %d --max_tip_len %s" % (base_path, overlaps, args.merge_contigs, min_overlap_len, args.threads, remove_branches, min_contig_len, max_tip_len), shell=True)
        os.remove(overlaps)
        os.chdir('..')
        subprocess.check_call("%s/scripts/fastq2fasta.py stage_c/singles.fastq contigs_stage_c.fasta" % base_path, shell=True)
        print "Done!"
        final_contig_file = "contigs_stage_c.fasta"
        # apply frequency-based filtering
        if args.filtering:
            try:
                freq_filtering("contigs_stage_c.fasta", "stage_c/singles.fastq", 0, input_info)
            except subprocess.CalledProcessError as e:
                print "\nKallisto not found - skipping this filtering step.\n"
    elif final_contig_file != "":
        if os.path.exists('stage_c'):
            shutil.rmtree('stage_c') #removes all the subdirectories!
        if os.path.exists('contigs_stage_c.fasta'):
            os.remove('contigs_stage_c.fasta')
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
        # find contig overlaps
        try:
            sfo_mm = 1 + (0.99 - args.merge_contigs) / (args.merge_contigs + 0.01)
            paired_count = 0
            singles_count = int(file_len('contigs_stage_c.fasta')/2)
            overlaps = run_sfo('c', sfo_mm, base_path, min_overlap_len, args.threads, singles_count, paired_count)
        except subprocess.CalledProcessError as e:
            print "\nRUNNING BLAST....\n"
            overlaps = run_blast('c', pident, base_path, min_overlap_len)
            subprocess.call("rm blastout* contigs_db*", shell=True)
        sys.stdout.flush()
        # run SAVAGE
        os.chdir('diploid')
        subprocess.check_call(["cp", overlaps, "original_overlaps.txt"])
        if args.use_subreads:
            subprocess.check_call("cp ../stage_c/subreads.txt subreads.txt", shell=True)
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --fastq ../stage_c --overlaps %s --merge_contigs %f --stage c --min_overlap_len %d --use_subreads --num_threads %d --remove_branches %s --min_read_len %d --diploid --max_tip_len %s" % (base_path, overlaps, args.merge_contigs, min_overlap_len, args.threads, remove_branches, min_contig_len, max_tip_len), shell=True)
        else:
            subprocess.check_call("%s/scripts/pipeline_per_stage.py --fastq ../stage_c --overlaps %s --merge_contigs %f --stage c --min_overlap_len %d --num_threads %d --remove_branches %s --min_read_len %d --diploid --max_tip_len %s" % (base_path, overlaps, args.merge_contigs, min_overlap_len, args.threads, remove_branches, min_contig_len, max_tip_len), shell=True)
        os.remove(overlaps)
        os.chdir('..')
        subprocess.check_call("%s/scripts/fastq2fasta.py diploid/singles.fastq diploid_contigs.fasta" % base_path, shell=True)
        print "Done!"
        final_contig_file = "diploid_contigs.fasta"
    elif final_contig_file != "":
        if os.path.exists('diploid'):
            shutil.rmtree('diploid') #removes all the subdirectories!
        if os.path.exists('diploid_contigs.fasta'):
            os.remove('diploid_contigs.fasta')

    # cleanup
    if os.path.exists('frequencies'):
        shutil.rmtree('frequencies') #removes all the subdirectories!

    if args.count_strains:
        if args.reference:
            # run strain count on final contig set and reference genome
            run_strain_count(args.reference, final_contig_file, base_path)
        else:
            print """
For estimating a lower bound on the number of strains a reference genome is
required; please specify a reference using --ref. If you want to keep your
current contigs, make sure to add the --no_assembly flag.
"""

    print """**************
SAVAGE assembly has been completed, the final contig set was written to:

        %s

Optionally, you can now apply frequency estimation using freq-est.py. Please see
the manual page for more information: https://bitbucket.org/jbaaijens/savage.

Thank you for using SAVAGE!
    """ % final_contig_file

    # return to start directory
    os.chdir(start_dir)

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
    print "- Running rust-overlaps",
    sys.stdout.flush()
    sfo_err = 1 / sfo_mm
    if paired_count > 0:
        sfo_len = int(round(min_overlap_len / 2))
    else:
        sfo_len = min_overlap_len
    # run sfoverlap
    try:
#        print "sfo_err:", sfo_err
        subprocess.check_call("rust-overlaps -i -r -w %d s_p1_p2.fasta sfoverlaps.out %f %d" % (threads, sfo_err, sfo_len), shell=True)
    except subprocess.CalledProcessError as e:
        print "-> sfoverlap failed, running blast instead"
        subprocess.check_call("makeblastdb -in s_p1_p2.fasta -dbtype nucl -out s_p1_p2.db 1>/dev/null 2>&1", shell=True)
        subprocess.check_call("blastn -db s_p1_p2.db -query s_p1_p2.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen' -out blastout.tsv -perc_identity 98", shell=True)
        subprocess.check_call("%s/scripts/blast2sfo.py --in blastout.tsv --out sfoverlaps.out -m %s" % (base_path, sfo_len), shell=True)
        subprocess.check_call("rm blastout.tsv", shell=True)
    # run postprocessing scripts
    print "\b" * 15 + "Processing output",
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
        sys.stderr.write("""\nERROR: Reference fasta not found: %s \nPlease enter full path to file.\n""" % reference)
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

def run_sfo(previous_stage, sfo_mm, base_path, min_overlap_len, threads, singles_count, paired_count):
    overlaps_file = "contig_overlaps.txt"
    sfo_err = 1 / sfo_mm
#    print "sfo_err: ", sfo_err
    subprocess.check_call("rust-overlaps -w %d -i -r contigs_stage_%s.fasta sfoverlaps.out %f %d" % (threads, previous_stage, sfo_err, min_overlap_len), shell=True)
    subprocess.check_call("%s/scripts/sfo2overlaps.py --in sfoverlaps.out --out %s --num_singles %d --num_pairs %d 1> /dev/null" % (base_path, overlaps_file, singles_count, paired_count), shell=True)
    subprocess.check_call("rm sfoverlaps.out", shell=True)
    overlaps_path = "../" + overlaps_file
    return overlaps_path

def freq_filtering(contig_fasta, contig_fastq, min_TPM, input_info): # fragmentsize, stddev, forward, reverse=""):
    print "\n--- Filtering contigs ---"
    # run kallisto
    abundance_file = run_kallisto(contig_fasta, input_info)
    # read TPMs from file
    TPM_dict = {}
    with open(abundance_file, 'r') as f:
        c = 0
        for line in f:
            c += 1
            if c == 1:
                continue
            [target_id, length, eff_length, est_counts, tpm] = line.split('\t')
            TPM_dict[target_id] = float(tpm)

    input_count = 0
    output_count = 0

    if min(TPM_dict.itervalues()) > min_TPM:
        print "Nothing filtered out.\n"
    else:
        # rename old contig files
        contigs_name, extension = os.path.splitext(contig_fasta)
        renamed_fasta = contigs_name + ".unfiltered.fasta"
        renamed_fastq = contigs_name + ".unfiltered.fastq"
        subprocess.check_call(['mv', contig_fasta, renamed_fasta])
        subprocess.check_call(['mv', contig_fastq, renamed_fastq])
        # filter contig set and write to new files
        output_fasta = open(contig_fasta, 'w')
        output_fastq = open(contig_fastq, 'w')
        with open(renamed_fastq, 'r') as f:
            c = 0
            for line in f:
                c += 1
                if (c % 4) == 1:
                    # ID line
                    id_line = line
                    cur_id = id_line.lstrip('@').rstrip('\n')
                elif (c % 4) == 2:
                    # seq line
                    seq_line = line
                elif (c % 4) == 3:
                    # + line
                    continue
                else:
                    # qual line
                    assert (c % 4) == 0
                    qual_line = line
                    input_count += 1
                    if TPM_dict[cur_id] > min_TPM:
                        # write to new files
                        output_fasta.write('>' + cur_id + '\n' + seq_line)
                        output_fastq.write(id_line + seq_line + '+\n' + qual_line)
                        output_count += 1
        output_fasta.close()
        output_fastq.close()
        os.remove(renamed_fasta)
        os.remove(renamed_fastq)
        print "Filtered %s down to %s contigs.\n" % (contig_fasta, output_count)
    return

def run_kallisto(contigs, input_info): #fragmentsize, stddev, forward, reverse=""):
    kallisto = "kallisto" # kallisto executable
    FNULL = open(os.devnull, 'w')
    # create output directory
    subprocess.call(['mkdir', '-p', 'frequencies'])
    # index construction
    contigs_name, extension = os.path.splitext(contigs)
    index_file = 'frequencies/' + contigs_name + '.idx'
    print "Kallisto index construction... "
    subprocess.check_call([kallisto, 'index', '-i', index_file, contigs], stdout=FNULL, stderr=FNULL)
    # estimate abundances
    print "Kallisto abundance quantification... "
    kallisto_out = 'frequencies/' + contigs_name
    if input_info.input_s:
        fragmentsize = str(input_info.fragmentsize)
        stddev = str(input_info.stddev)
        if input_info.input_p1 and input_info.input_p2:
            subprocess.check_call([kallisto, 'quant', '-i', index_file, '-o', kallisto_out, '-b', '100', '-l', fragmentsize, '-s', stddev, '--single', input_info.input_s, input_info.input_p1, input_info.input_p2], stdout=FNULL, stderr=FNULL)
        else:
            subprocess.check_call([kallisto, 'quant', '-i', index_file, '-o', kallisto_out, '-b', '100', '-l', fragmentsize, '-s', stddev, '--single', input_info.input_s], stdout=FNULL, stderr=FNULL)
    else:
        subprocess.check_call([kallisto, 'quant', '-i', index_file, '-o', kallisto_out, '-b', '100', input_info.input_p1, input_info.input_p2], stdout=FNULL, stderr=FNULL)
    # if reverse:
    #     subprocess.check_call([kallisto, 'quant', '-i', index_file, '-o', kallisto_out, '-b', '100', '-l', fragmentsize, '-s', stddev, forward, reverse])
    # else:
    #     subprocess.check_call([kallisto, 'quant', '-i', index_file, '-o', kallisto_out, '-b', '100', '-l', fragmentsize, '-s', stddev, '--single', forward])
    abundance_file = kallisto_out + '/abundance.tsv'
    FNULL.close()
    return abundance_file

def run_strain_count(reference, contigs_fastq, base_path):
    print "Estimating strain count on %s" % contigs_fastq
    # align contigs to ref
    contigs_sam = contigs_fastq.rstrip('.fastq') + ".sam"
    subprocess.check_call("bwa mem %s %s 1> %s 2> /dev/null" % (reference, contigs_fastq, contigs_sam), shell=True)
    subprocess.check_call("%s/estimate_strain_count.py --sam %s --ref %s" % (base_path, contigs_sam, reference), shell=True)
    return


if __name__ == '__main__':
    sys.exit(main())
