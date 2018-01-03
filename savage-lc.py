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
from math import floor

# ------------------------------

__author__ = "Jasmijn Baaijens"
__license__ = "GPL"

version = "0.1.0"
releasedate = ""

usage = """
Program: SAVAGE-LC - Strain Aware Assembly on Low Coverage Data
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
original_fastq = ""
denovo = True

def main():
    parser = ArgumentParser(description=usage, formatter_class=RawTextHelpFormatter)
    basic = parser.add_argument_group('basic arguments')
    basic.add_argument('-s', dest='input_s', type=str, help='path to input fastq containing single-end reads')
    basic.add_argument('-p1', dest='input_p1', type=str, help='path to input fastq containing paired-end reads (/1)')
    basic.add_argument('-p2', dest='input_p2', type=str, help='path to input fastq containing paired-end reads (/2)')
    basic.add_argument('-m', '--min_overlap_len', dest='min_overlap_len', type=int, default=50, help='minimum overlap length required between reads')
    basic.add_argument('-m_EC', '--min_overlap_len_EC', dest='min_overlap_len_EC', type=int, help='minimum overlap length required between reads during error correction')
    basic.add_argument('-t', '--num_threads', dest='threads', type=int, default=1, help='allowed number of cores')
#    basic.add_argument('--revcomp', dest='revcomp', action='store_true', help='use this option when paired-end input reads are in forward-reverse orientation;\nthis option will take reverse complements of /2 reads (specified with -p2)\nplease see the SAVAGE manual for more information about input read orientations')
    basic.add_argument('--hap_cov', dest='hap_cov', type=float, required=True, help='average coverage per haplotype')
    basic.add_argument('--insert_size', dest='insert_size', type=float, required=True, help='mean insert size for paired-end input')
    basic.add_argument('--stddev', dest='stddev', type=float, required=True, help='standard deviation of insert size for paired-end input')
    ref_guided = parser.add_argument_group('reference-guided mode')
    ref_guided.add_argument('--ref', dest='reference', type=str, help='reference genome in fasta format')
    ref_guided.add_argument('--ref_guided_mode', dest='ref_guided_mode', action='store_true', help='perform reference-guided assembly')
    advanced = parser.add_argument_group('advanced arguments')
    advanced.add_argument('--no_EC', dest='error_correction', action='store_false', help='skip error correction in initial iteration (i.e. no cliques)')
    advanced.add_argument('--no_overlaps', dest='compute_overlaps', action='store_false', help='skip overlap computations (use existing overlaps file instead)')
    advanced.add_argument('--no_preprocessing', dest='preprocessing', action='store_false', help='skip preprocessing procedure')
    advanced.add_argument('--no_assembly', dest='assembly', action='store_false', help='skip all assembly steps; only use this option when using --count_strains separate from assembly (e.g. on a denovo assembly)')
    advanced.add_argument('--count_strains', dest='count_strains', action='store_true', help='compute a lower bound on the number of strains in this sample; note: this requires a reference genome.')
    advanced.add_argument('--mismatch_rate', dest='merge_contigs', type=float, default=0.0, help='specify maximal distance between contigs for merging into master strains (stage c)')
    advanced.add_argument('--min_clique_size', dest='min_clique_size', type=int, default=3, help='minimum clique size used during error correction')
    advanced.add_argument('--sfo_err', dest='sfo_err', type=float, default=0.02, help='input parameter for sfo: maximal mismatch rate')
    advanced.add_argument('--diploid', dest='diploid', action='store_true', help='use this option for diploid genome assembly')
    advanced.add_argument('--diploid_contig_len', dest='diploid_contig_len', type=int, default=0, help='minimum contig length required for diploid step contigs')
    advanced.add_argument('--diploid_overlap_len', dest='diploid_overlap_len', type=int, help='min_overlap_len used in diploid assembly step; by default equal to assembly min_overlap_len.')
    advanced.add_argument('--average_read_len', dest='average_read_len', type=float, help='average length of the input reads; will be computed from the input if not specified')
#    advanced.add_argument('--no_filtering', dest='filtering', action='store_false', help='disable kallisto-based filtering of contigs')
    advanced.add_argument('--max_tip_len', dest='max_tip_len', type=int, help='maximum extension length for a sequence to be called a tip')
    #advanced.add_argument('--min_evidence', dest='min_evidence', type=int, required=True, help='minimum number of uniquely matching reads to resolve branches')
    split_mode = parser.add_argument_group('split mode')
    split_mode.add_argument('--original_SE_count', dest='original_SE_count', type=int, default=-1, help='number of single-end sequences in input before splitting')
    split_mode.add_argument('--original_PE_count', dest='original_PE_count', type=int, default=-1, help='number of paired-end sequences in input before splitting')
    split_mode.add_argument('--original_fastq', dest='original_fastq', type=str, default='', help='original fastq file before splitting')

    # store current working directory
    cwd = os.getcwd()
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
------------------------------------------------------
SAVAGE-LC - Strain Aware Assembly on Low Coverage Data
------------------------------------------------------
Version: %s
Author: %s
    """ % (version, __author__)
    print "Command used:"
    print ' '.join(sys.argv)
    print

    FNULL = open(os.devnull, 'w')

    # always perform de novo assembly; reference is only used for strain counting
    global denovo;
    denovo = False if args.ref_guided_mode else True

    if args.reference:
        if not os.path.exists(args.reference):
            sys.stderr.write("""\nERROR: Reference fasta not found: %s \nPlease enter full path to file.\n""" % args.reference)
            sys.stderr.flush()
            sys.exit(1)
        if not os.path.exists(args.reference + ".bwt"):
            print "Building index for reference genome...",
            subprocess.check_call("bwa index %s 1>/dev/null 2>&1" % args.reference, shell=True)
            print "done!\n"

    if not args.assembly:
        print "Skipping assembly because --no_assembly flag was used.\n"

    if not (args.preprocessing or args.compute_overlaps or args.assembly or args.diploid):
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
        if os.path.exists('contigs_diploid.fasta'):
            final_contig_file = 'contigs_diploid.fasta'
        elif os.path.exists('contigs.fasta'):
            final_contig_file = 'contigs.fasta'
        else:
            sys.stderr.write("No contigs found for estimating strain count. Please run savage assembly first.\n")
            sys.stderr.flush()
            sys.exit(1)

        if args.reference:
            # run estimate_strain_count
            run_strain_count(args.reference, final_contig_file, base_path)
            FNULL.close()
            return # No assembly, so we are done!
        else:
            sys.stderr.write("Please specify a reference genome when using --count_strains.\n")
            sys.stderr.flush()
            sys.exit(1)
    else:
        final_contig_file = ""

    if not args.preprocessing:
        preprocessing = False
    elif (args.assembly or args.diploid) and not args.compute_overlaps:
        preprocessing = False
    elif args.diploid and not args.assembly:
        preprocessing = False
    else:
        preprocessing = True

    if not args.compute_overlaps:
        compute_overlaps = False
    elif args.diploid and not args.assembly:
        compute_overlaps = False
    else:
        compute_overlaps = True

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
    print "Number of paired-end reads = %d" % int(p_seq_count/2)
    print "Total number of bases =", total_seq_len
    print "Average sequence length = %.1f" % average_read_len
    print

    global original_fastq;
    if args.original_fastq != '':
        original_fastq = args.original_fastq
    else:
        original_fastq = cwd + "/assembly/s_p1_p2.fastq"

    if args.max_tip_len is None:
        max_tip_len = int(round(average_read_len)) # average input sequence length
        # if p_seq_count > s_seq_count:
        #     max_tip_len = 2*int(round(average_read_len)) # average input fragment length
        # else:
        #     max_tip_len = int(round(average_read_len)) # average input sequence length
        print "Using max_tip_len =", max_tip_len
        assert max_tip_len > 0
    elif args.max_tip_len >= 0:
        max_tip_len = args.max_tip_len
    else:
        sys.stderr.write("""\nERROR: invalid --max_tip_len %s. Please make sure this value is non-negative. Exiting.\n\n""" % args.max_tip_len)
        sys.stderr.flush()
        sys.exit(1)

    if args.min_overlap_len_EC is None:
        #m = 0.6*average_read_len # 60% of average input read length
        m = 2+0.5*average_read_len
        min_overlap_len_EC = int(round(m))
        print "Using min_overlap_len_EC =", min_overlap_len_EC
        print
        assert min_overlap_len_EC > 0
    else:
        min_overlap_len_EC = args.min_overlap_len_EC

    if min_overlap_len_EC < 80:
        print "----------------------------------------------------------------"
        print "WARNING: min_overlap_len_EC = %s" % min_overlap_len_EC
        print "For more accurate error correction, increase the minimal overlap length used during error correction using --min_overlap_len_EC"
        print "----------------------------------------------------------------"

    # Preprocessing: rename and reorganize reads
    if preprocessing:
        print "*******************"
        print "Preprocessing input"
        overwrite_dir('assembly')
        sys.stdout.flush()
        # combine all reads into one big fastq file, concatenating singles-paired1-paired2
        if args.input_s:
            if args.input_p1 and args.input_p2:
                subprocess.check_call("cat %s %s %s > tmp.fastq" % (args.input_s, args.input_p1, args.input_p2), shell=True)
            else:
                subprocess.check_call("cat %s > tmp.fastq" % (args.input_s), shell=True)
        else:
            subprocess.check_call("cat %s %s > tmp.fastq" % (args.input_p1, args.input_p2), shell=True)
        subprocess.check_call("%s/scripts/rename_fas.py --in tmp.fastq --out assembly/s_p1_p2.fastq" % (base_path), shell=True)
        subprocess.check_call("rm tmp.fastq", shell=True)
        # align reads to reference if reference-guided
        if not denovo:
            # Check for reference fasta
            if not os.path.exists(args.reference):
                sys.stderr.write("""\nERROR: Reference fasta not found: %s \nPlease enter full path to file.\n""" % args.reference)
                sys.stderr.flush()
                sys.exit(1)
            # Run BWA to get alignments
            fastq_path = "assembly"
            try: # just try if reference has already been indexed
                subprocess.check_call("bwa mem %s %s/s_p1_p2.fastq 1> %s/s_p1_p2.sam 2> /dev/null" % (args.reference, fastq_path, fastq_path), shell=True)
            except subprocess.CalledProcessError as e:
                subprocess.check_call("bwa index %s 1>/dev/null 2>&1" % args.reference, shell=True)
                subprocess.check_call("bwa mem %s %s/s_p1_p2.fastq 1> %s/s_p1_p2.sam 2> /dev/null" % (args.reference, fastq_path, fastq_path), shell=True)
        print "\rDone!" + ' ' * 40
        sys.stdout.flush()

    # Compute all suffix-prefix overlaps
    overlaps = "original_overlaps.txt"
    if compute_overlaps:
        print "********************"
        print "Overlap computations"
        sys.stdout.flush()
        os.chdir('assembly')
        # find all overlaps
        if denovo:
            preprocessing_denovo(min_overlap_len_EC, args.sfo_err, args.threads, base_path)
        else:
            preprocessing_ref(min_overlap_len_EC, args.reference, base_path)
        os.chdir('..')
        print "- Done!"
        sys.stdout.flush()
    elif args.assembly:
        # check if overlaps file already exists
        if not os.path.exists('assembly/original_overlaps.txt'):
            sys.stderr.write("""Assuming existing overlaps file 'assembly/original_overlaps.txt'.
                     Please make sure this file exists, or run overlap computations first\n""")
            sys.stderr.flush()
            sys.exit(1)

    # settings required for savage-lc-split
    if args.original_fastq != '':
        SE_count = args.original_SE_count
        PE_count = args.original_PE_count
    else:
        SE_count = s_seq_count
        PE_count = p_seq_count
    original_readcount = SE_count + PE_count

    assert args.insert_size > 0
    readlen = average_read_len
    intseg = args.insert_size - 2*readlen # internal segment size, can be negative
    stddev = args.stddev
    hcov = args.hap_cov

    # Run SAVAGE Stage a: error correction and initial contig formation
    if args.assembly:
        print "***************"
        print "SAVAGE assembly"
        sys.stdout.flush()
        os.chdir('assembly')
        diploid = "false"
        verbose = "false"
        error_rate = 0
        min_read_len = 0
        # min_evidence = args.min_evidence
        branch_reduction = [args.hap_cov, SE_count, int(PE_count/2)]
        EC = "true" if args.error_correction else "false"
        build_threshold_table(readlen, intseg, stddev, hcov, base_path)
        run_savage_assembly(EC, args.min_overlap_len, min_overlap_len_EC, args.min_clique_size, min_read_len, max_tip_len, branch_reduction, error_rate, args.threads, original_readcount, diploid, base_path, verbose)
        os.chdir('..')
        subprocess.check_call("%s/scripts/fastq2fasta.py assembly/singles.fastq contigs.fasta" % base_path, shell=True)
        print "Done!"
        final_contig_file = "contigs.fasta"
    elif not (args.diploid or args.count_strains):
        if os.path.exists('contigs.fasta'):
            os.remove('contigs.fasta')

    if args.diploid:
        # final diploid contig merging
        print "**************"
        print "SAVAGE diploid"
        # parameters
        if args.diploid_overlap_len:
            min_overlap_len = args.diploid_overlap_len
        else:
            min_overlap_len = args.min_overlap_len
        min_contig_len = args.diploid_contig_len
        # prepare input files
        overwrite_dir('diploid')
        if not (os.path.exists('assembly/singles.fastq') and os.path.exists('contigs.fasta')):
            print """Contigs file from assembly stage not found. Please make sure that both
                     'assembly/singles.fastq' and 'contigs.fasta' exist. If absent, please
                     rerun assembly."""
            sys.exit(1)
        elif os.stat("contigs.fasta").st_size == 0:
            print "Empty set of contigs from assembly stage (contigs.fasta) --> Exiting SAVAGE."
            sys.exit(1)
        subprocess.call(['cp', 'assembly/singles.fastq', 'diploid/s_p1_p2.fastq'], stdout=FNULL, stderr=FNULL)
        pident = 98
        # find contig overlaps
        paired_count = 0
        contig_count = int(file_len('contigs.fasta')/2)
        os.chdir('diploid')
        reversals = True
        run_sfo('../contigs.fasta', args.merge_contigs, base_path, min_overlap_len, args.threads, contig_count, paired_count, "original_overlaps.txt", reversals)
        sys.stdout.flush()
        # run SAVAGE
        subprocess.check_call("cp ../assembly/subreads.txt subreads.txt", shell=True)
        diploid = "true"
        verbose = "false"
        EC = "false"
        error_rate = args.merge_contigs
        min_read_len = 0
        # min_evidence = max(1, args.min_evidence-1)
        min_clique_size = 2
        branch_reduction = [args.hap_cov, SE_count, int(PE_count/2)]
        build_threshold_table(readlen, intseg, stddev, hcov, base_path)
        run_savage_assembly(EC, min_overlap_len, min_overlap_len, min_clique_size, min_read_len, max_tip_len, branch_reduction, error_rate, args.threads, original_readcount, diploid, base_path, verbose)
        os.chdir('..')
        if os.path.exists('diploid/singles.fastq'):
            subprocess.check_call("%s/scripts/fastq2fasta.py diploid/singles.fastq contigs_diploid.fasta" % base_path, shell=True)
            final_contig_file = "contigs_diploid.fasta"
            print "Done!"
        else:
            print "No contigs produced in diploid stage."
    elif final_contig_file != "":
        if os.path.exists('diploid'):
            shutil.rmtree('diploid') #removes all the subdirectories!
        if os.path.exists('contigs_diploid.fasta'):
            os.remove('contigs_diploid.fasta')

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

    FNULL.close()
    return


# ------------------------------

def file_len(fname):
    i = 0
    if os.path.exists(fname):
        with open(fname, 'r') as f:
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
#        subprocess.check_call("rm -rf %s" %dir, shell=True)
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
    if not os.path.isfile(filename):
        return [0, 0, 0]
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

def preprocessing_denovo(min_overlap_len, sfo_err, threads, base_path):
    sys.stdout.flush()
    # prepare fasta
    subprocess.check_call("%s/scripts/fastq2fasta.py s_p1_p2.fastq s_p1_p2.fasta" % base_path, shell=True)
    read_count = int(file_len('s_p1_p2.fastq')/4.0)
    # run SFO
    print "Running rust-overlaps",
    sys.stdout.flush()
    sfo_len = min_overlap_len
    # run sfoverlap
    try:
        subprocess.check_call("rust-overlaps -i -r -w %d s_p1_p2.fasta sfoverlaps.out %f %d 1> /dev/null" % (threads, sfo_err, sfo_len), shell=True)
    except subprocess.CalledProcessError as e:
        print "ERROR: rust-overlaps not installed"
        sys.exit(1)
    # run postprocessing scripts
    print "- Processing output",
    sys.stdout.flush()
    subprocess.check_call("%s/scripts/sfo2overlaps.py --in sfoverlaps.out --out original_overlaps.txt --num_singles %d --num_pairs 0 1> /dev/null" % (base_path, read_count), shell=True)
#    print "Overlaps are ready!"
    subprocess.check_call("rm sfoverlaps.out", shell=True)
    return

def preprocessing_ref(min_overlap_len, reference, base_path):
    print "Reference-guided overlap computations",
    sys.stdout.flush()
    # Check for reference fasta
    if not os.path.exists(reference):
        sys.stderr.write("""\nERROR: Reference fasta not found: %s \nPlease enter full path to file.\n""" % reference)
        sys.stderr.flush()
        sys.exit(1)
    # Induce overlaps from alignment
    subprocess.check_call("%s/scripts/sam2overlaps.py --sam_s s_p1_p2.sam --ref %s --min_overlap_len %d --out original_overlaps.txt" %(base_path, reference, min_overlap_len), shell=True)
    return

def run_sfo(fasta, sfo_err, base_path, min_overlap_len, threads, s_count, p_count, overlaps_file, reversals):
#    print "sfo_err: ", sfo_err
    if reversals:
        subprocess.check_call("rust-overlaps -w %d -i -r %s sfoverlaps.out %f %d 1> /dev/null" % (threads, fasta, sfo_err, min_overlap_len), shell=True)
    else:
        subprocess.check_call("rust-overlaps -w %d -i %s sfoverlaps.out %f %d 1> /dev/null" % (threads, fasta, sfo_err, min_overlap_len), shell=True)
    subprocess.check_call("%s/scripts/sfo2overlaps.py --in sfoverlaps.out --out %s --num_singles %d --num_pairs %d 1> /dev/null" % (base_path, overlaps_file, s_count, p_count), shell=True)
    subprocess.check_call("rm sfoverlaps.out", shell=True)
    return

def run_strain_count(reference, contigs_fastq, base_path):
    print "Estimating strain count on %s" % contigs_fastq
    # align contigs to ref
    contigs_sam = contigs_fastq.rstrip('.fastq') + ".sam"
    try:
        subprocess.check_call("bwa mem %s %s 1> %s 2> /dev/null" % (reference, contigs_fastq, contigs_sam), shell=True)
    except subprocess.CalledProcessError as e:
        subprocess.check_call("bwa index %s 1>/dev/null 2>&1" % reference, shell=True)
        subprocess.check_call("bwa mem %s %s 1> %s 2> /dev/null" % (reference, contigs_fastq, contigs_sam), shell=True)
    subprocess.check_call("%s/estimate_strain_count.py --sam %s --ref %s" % (base_path, contigs_sam, reference), shell=True)
    return

def run_savage_assembly(EC, min_overlap_len, min_overlap_len_EC, min_clique_size, min_read_len, max_tip_len, branch_reduction, error_rate, threads, original_readcount, diploid, base_path, verbose=False):
    # create a global log file; after every iteration the log file is appended to this global log file
    FNULL = open(os.devnull, 'w')
    subprocess.call(["rm", "pipeline.log"], stdout=FNULL, stderr=FNULL)
    subprocess.call(["touch", "pipeline.log"])
    # remove existing stats file
    subprocess.call(["rm", "stats.txt"], stdout=FNULL, stderr=FNULL)
    subprocess.call(["touch", "stats.txt"])
    # remove existing tips file
    subprocess.call(["rm", "removed_tip_sequences.fastq"], stdout=FNULL, stderr=FNULL)
    subprocess.call(["touch", "removed_tip_sequences.fastq"])
    # keep stats to find out when the algo has converged
    original_overlaps = int(file_len('original_overlaps.txt'))
    const_read_its = 0
    iteration = 1
    read_counts = [original_readcount]
    overlap_counts = [original_overlaps]
    edge_counts = []
    max_read_lengths = []
    stats = [iteration, read_counts, overlap_counts, edge_counts, max_read_lengths]
    branch_red = [0, 0, 0]
    final_it = False
    remove_tips = "false"
    # run first iteration (error correction)
    if EC == "true":
        edge_threshold = 0.95 #0.97
        first_it = "true"
        cliques = "true"
        stats = run_viralquasispecies(stats, "s_p1_p2.fastq", "original_overlaps.txt", min_overlap_len_EC, min_overlap_len, min_clique_size, edge_threshold, min_read_len, max_tip_len, first_it, cliques, EC, branch_red, error_rate, threads, original_readcount, diploid, verbose, final_it, remove_tips)
    elif diploid == "true":
        edge_threshold = 1 # select edges based on mismatch rate rather than overlap score
        first_it = "false"
        cliques = "true"
        branch_red = branch_reduction
        stats = run_viralquasispecies(stats, "s_p1_p2.fastq", "original_overlaps.txt", min_overlap_len, min_overlap_len, min_clique_size, edge_threshold, min_read_len, max_tip_len, first_it, cliques, EC, branch_red, error_rate, threads, original_readcount, diploid, verbose, final_it, remove_tips)
    else:
        edge_threshold = 0.95 #0.97
        first_it = "true"
        cliques = "true"
        stats = run_viralquasispecies(stats, "s_p1_p2.fastq", "original_overlaps.txt", min_overlap_len_EC, min_overlap_len, min_clique_size, edge_threshold, min_read_len, max_tip_len, first_it, cliques, EC, branch_red, error_rate, threads, original_readcount, diploid, verbose, final_it, remove_tips)

    # now iterate until convergence
    [iteration, read_counts, overlap_counts, edge_counts, max_read_lengths] = stats # update stats
    first_it = "false"
    EC = "false"
    min_clique_size = 2
    edge_threshold = 1
    while read_counts[-1] > 0 and overlap_counts[-1] > 0 and edge_counts[-1] > 0 and const_read_its < 2:
        while read_counts[-1] > 0 and overlap_counts[-1] > 0 and edge_counts[-1] > 0 and const_read_its < 2:
            # merge simple paths
            cliques = "false"
            branch_red = [0, 0, 0]
            stats = run_viralquasispecies(stats, "singles.fastq", "overlaps.txt", min_overlap_len, min_overlap_len, min_clique_size, edge_threshold, min_read_len, max_tip_len, first_it, cliques, EC, branch_red, error_rate, threads, original_readcount, diploid, verbose, final_it, remove_tips)
            [iteration, read_counts, overlap_counts, edge_counts, max_read_lengths] = stats # update stats
            if read_counts[-1] == read_counts[-2]:
                const_read_its += 1
            else:
                const_read_its = 0
        # apply read-based branch reduction and merge along remaining branches
        cliques = "true"
        branch_red = branch_reduction
        print "iteration %d -> BranchReduction" % iteration
        stats = run_viralquasispecies(stats, "singles.fastq", "overlaps.txt", min_overlap_len, min_overlap_len, min_clique_size, edge_threshold, min_read_len, max_tip_len, first_it, cliques, EC, branch_red, error_rate, threads, original_readcount, diploid, verbose, final_it, remove_tips)
        [iteration, read_counts, overlap_counts, edge_counts, max_read_lengths] = stats # update stats
        if read_counts[-1] == read_counts[-2]:
            const_read_its += 1
        else:
            const_read_its = 0

#     if diploid=="true":
#         # diploid final round with a relaxed evidence threshold after removing tips
#         const_read_its = 0
#         branch_reduction[0] = 1 #int(floor(branch_reduction[0]/2))
# #        remove_tips = "true"
#         while read_counts[-1] > 0 and overlap_counts[-1] > 0 and const_read_its < 2:
#             # apply read-based branch reduction and merge along remaining branches
#             cliques = "true"
#             branch_red = branch_reduction
#             print "iteration %d -> BranchReduction" % iteration
#             stats = run_viralquasispecies(stats, "singles.fastq", "overlaps.txt", min_overlap_len, min_overlap_len, min_clique_size, edge_threshold, min_read_len, max_tip_len, first_it, cliques, EC, branch_red, error_rate, threads, original_readcount, diploid, verbose, final_it, remove_tips)
#             [iteration, read_counts, overlap_counts, edge_counts, max_read_lengths] = stats # update stats
#             if read_counts[-1] == read_counts[-2]:
#                 const_read_its += 1
#             else:
#                 const_read_its = 0

    # one final iteration to remove singletons and tips (?) from contig file
    final_it = True
    cliques = "false"
    remove_tips = "false" #"true"
    branch_red = [0, 0, 0]
    if read_counts[-1] > 0:
        stats = run_viralquasispecies(stats, "singles.fastq", "overlaps.txt", min_overlap_len, min_overlap_len, min_clique_size, edge_threshold, min_read_len, max_tip_len, first_it, cliques, EC, branch_red, error_rate, threads, original_readcount, diploid, verbose, final_it, remove_tips)

    print "Assembly done in %d iterations" % iteration
    print "Maximum read length per iteration: \t", max_read_lengths
    print "Number of contigs per iteration: \t", read_counts[1:]
    print "Number of overlaps per iteration: \t", overlap_counts
    FNULL.close()
    return

def run_viralquasispecies(stats, fastq, overlaps, min_overlap_len, next_min_overlap, min_clique_size, edge_threshold, min_read_len, max_tip_len, first_it, cliques, EC, branch_reduction, error_rate, threads, original_readcount, diploid, verbose, final_it, remove_tips):
    [iteration, read_counts, overlap_counts, edge_counts, max_read_lengths] = stats
    selfpath = sys.path[0]
    viralquasispecies = selfpath + "/bin/ViralQuasispecies"
    COPYFILES = False
    if EC=="true":
        keep_singletons = 1000
    elif diploid=="true" and final_it:
        #max_tip_len *= 2
        keep_singletons = max_tip_len
        #keep_singletons = 0
    else:
        keep_singletons = 0
#    fno = 3 if cliques=="true" else 1
    fno = 1
    remove_trans = 2 if EC=="true" else 1
    separate_tips = "true" if final_it else "false"
    remove_inclusions = "true" if (final_it and diploid=="true") else "false"
#    remove_tips = "true" if final_it and diploid=="true" else "false"
    [hap_cov, branch_SE_c, branch_PE_c] = branch_reduction
    if (cliques=="false" or (EC=="false" and hap_cov==0)):
        remove_branches = "true"
    else:
        remove_branches = "false"
    if verbose == 'true':
        print "\n*********************"
        print "**** Iteration %d ****" %iteration
        print "*********************"
    shell_command = [viralquasispecies,
        "--singles=%s" %fastq,
        "--overlaps=%s" %overlaps,
        "--threads=%d" %threads,
        "--edge_threshold=%f" %edge_threshold,
        "--first_it=%s" %first_it,
        "--cliques=%s" %cliques,
        "--error_correction=%s" %EC,
        "--keep_singletons=%d" %keep_singletons,
        "--min_clique_size=%d" %min_clique_size,
        "--remove_branches=%s" %remove_branches,
        "--remove_tips=%s" %remove_tips,
        "--min_overlap_len=%d" %min_overlap_len,
        "--merge_contigs=%f" %error_rate,
        "--FNO=%d" %fno,
        "--original_readcount=%d" %original_readcount,
        "--remove_trans=%d" %remove_trans,
        "--optimize=false",
        "--verbose=%s" %verbose,
        "--base_path=%s" % selfpath,
        "--min_read_len=%s" % min_read_len,
        "--max_tip_len=%s" % max_tip_len,
        "--separate_tips=%s" % separate_tips,
        "--ignore_inclusions=%s" % remove_inclusions,
        "--diploid=%s" % diploid,
        "--min_qual=0" # never insert N's
    ]
    if hap_cov > 0:
        shell_command.append("--branch_reduction=true")
        shell_command.append("--original_fastq=%s" % original_fastq)
        shell_command.append("--branch_SE_c=%s" % branch_SE_c)
        shell_command.append("--branch_PE_c=%s" % branch_PE_c)
    subprocess.check_call(shell_command)
    if COPYFILES:
        copy_files(stats[0])
    copy_log()
    # recompute overlaps on contigs
    if denovo:
        s_count = int(file_len('singles.fastq')/4)
        if s_count > 0:
            subprocess.check_call("%s/scripts/fastq2fasta.py singles.fastq singles.fasta" % selfpath, shell=True)
            sfo_err = 0
            reversals = True #if first_it=="true" else False # allow reversals only in the first two iterations
            run_sfo("singles.fasta", sfo_err, selfpath, next_min_overlap, threads, s_count, 0, "overlaps.txt", reversals)
    # update pipeline statistics
    stats = analyze_results(stats)
    if verbose == 'true':
        print "***"
    return stats

def copy_files(it):
    subprocess.call(["cp", "singles.fastq", "it%d_singles.fastq" %it])
    subprocess.call(["cp", "overlaps.txt", "it%d_overlaps.txt" %it])
    subprocess.call(["cp", "subreads.txt", "it%d_subreads.txt" %it])
    subprocess.call(["cp", "graph.gfa", "it%d_graph.gfa" %it])
    subprocess.call(["cp", "graph_trimmed.gfa", "it%d_graph_trimmed.gfa" %it])

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

def analyze_results(stats):
    [iteration, read_counts, overlap_counts, edge_counts, max_read_lengths] = stats
    iteration += 1
    [contig_count, total_len, max_len] = analyze_fastq("singles.fastq")
    max_read_lengths.append(max_len)
    read_counts.append(contig_count)
    if os.path.isfile('overlaps.txt'):
        n_overlaps = int(file_len('overlaps.txt'))
    else:
        n_overlaps = 0
    overlap_counts.append(n_overlaps)
    n_edges = get_edge_count()
    edge_counts.append(n_edges)
    return [iteration, read_counts, overlap_counts, edge_counts, max_read_lengths]

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
        return max_cov
    else:
        return 0

def build_threshold_table(readlen, intseg, stddev, hcov, base_path):
    filename = "evidence_threshold_table.tsv"
    if os.path.isfile(filename):
        # threshold table file exists, check if it has the right parameters
        params = {}
        with open(filename, 'r') as f:
            for line in f:
                if line[0] != '#':
                    continue
                line = line.rstrip('\n').split()
                if len(line) != 3:
                    continue
                params[line[1]] = float(line[2])
        if (readlen == params["readlen"] and intseg == params["intseg"] and
                stddev == params["stddev"] and hcov == params["hcov"]):
            # correct file exists
            return
    print "Building the evidence threshold table..."
    # build a new table
    shell_command = "{}/scripts/min_ev_table.py".format(base_path)
    shell_command += " -l {}".format(readlen)
    shell_command += " -i {}".format(intseg)
    shell_command += " --stddev {}".format(stddev)
    shell_command += " --hcov {}".format(hcov)
    shell_command += " -o {}".format(filename)
    subprocess.check_call(shell_command, shell=True)
    return



if __name__ == '__main__':
    sys.exit(main())
