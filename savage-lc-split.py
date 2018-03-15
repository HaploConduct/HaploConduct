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
from multiprocessing import Pool
from functools import partial

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
    ref_guided_split = parser.add_argument_group('reference-guided vertical data splitting')
    ref_guided_split.add_argument('--ref', dest='reference', type=str, required=True, help='reference genome in fasta format')
    ref_guided_split.add_argument('--split_size', dest='split_size', type=int, default=10000, help='size of regions over which the reads are divided')
    ref_guided_split.add_argument('--split_overlap', dest='split_overlap', type=int, default=1000, help='size of overlap between regions over which the reads are divided')
    ref_guided_split.add_argument('--split_only', dest='split_only', action='store_true', help="don't do assembly, only data splitting")
    advanced = parser.add_argument_group('advanced arguments')
#    advanced.add_argument('--no_EC', dest='error_correction', action='store_false', help='skip error correction in initial iteration (i.e. no cliques)')
#    advanced.add_argument('--no_overlaps', dest='compute_overlaps', action='store_false', help='skip overlap computations (use existing overlaps file instead)')
#    advanced.add_argument('--no_preprocessing', dest='preprocessing', action='store_false', help='skip preprocessing procedure')
    advanced.add_argument('--no_assembly', dest='assembly', action='store_false', help='skip all assembly steps; only use this option when using --count_strains separate from assembly (e.g. on a denovo assembly)')
    advanced.add_argument('--count_strains', dest='count_strains', action='store_true', help='compute a lower bound on the number of strains in this sample; note: this requires a reference genome.')
    advanced.add_argument('--mismatch_rate', dest='merge_contigs', type=float, default=0.0, help='specify maximal distance between contigs for merging into master strains (stage c)')
    advanced.add_argument('--min_clique_size', dest='min_clique_size', type=int, default=3, help='minimum clique size used during error correction')
    advanced.add_argument('--sfo_err', dest='sfo_err', type=float, default=0.02, help='input parameter for sfo: maximal mismatch rate')
    advanced.add_argument('--diploid', dest='diploid', action='store_true', help='use this option for diploid genome assembly')
    advanced.add_argument('--diploid_contig_len', dest='diploid_contig_len', type=int, default=0, help='minimum contig length required for diploid step contigs')
    advanced.add_argument('--diploid_overlap_len', dest='diploid_overlap_len', type=int, default=30, help='min_overlap_len used in diploid assembly step')
    advanced.add_argument('--average_read_len', dest='average_read_len', type=float, help='average length of the input reads; will be computed from the input if not specified')
#    advanced.add_argument('--no_filtering', dest='filtering', action='store_false', help='disable kallisto-based filtering of contigs')
    advanced.add_argument('--max_tip_len', dest='max_tip_len', type=int, help='maximum extension length for a sequence to be called a tip')
#    advanced.add_argument('--min_evidence', dest='min_evidence', type=int, required=True, help='minimum number of uniquely matching reads to resolve branches')
    pool = parser.add_argument_group('thread pool settings')
    pool.add_argument('--pool_size', dest='pool_size', type=int, default=1, help='number of regions to be processed in parallel')
    pool.add_argument('--sfo_threads', dest='sfo_threads', type=int, help='number of threads used per region')

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
    denovo = True
    # analyze reference genome
    if not os.path.exists(args.reference):
        sys.stderr.write("""\nERROR: Reference fasta not found: %s \nPlease enter full path to file.\n""" % args.reference)
        sys.stderr.flush()
        sys.exit(1)
    if not os.path.exists(args.reference + ".bwt"):
        print "Building index for reference genome...",
        subprocess.check_call("bwa index %s 1>/dev/null 2>&1" % args.reference, shell=True)
        print "done!\n"
    ref_pieces = {} # store tuples (ID, length)
    chrom2regions = {}
    with open(args.reference) as f:
        l = 0
        for line in f:
            if line[0] == '>':
                if l > 0:
                    ref_pieces[ID] = l
                    chrom2regions[ID] = []
                ID = line.lstrip('>').rstrip('\n').split()[0]
                if any(c in ID for c in '|/\()}}{{[]'):
                    print "ERROR: chromosome id contains |/\()}}{{[]; Exiting."
                    sys.exit(1)
                l = 0
            else:
                l += len(line.rstrip('\n'))
        ref_pieces[ID] = l # add final reference line
        chrom2regions[ID] = []

    if not args.assembly:
        print "Skipping assembly because --no_assembly flag was used.\n"

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

    original_readcount = s_seq_count + p_seq_count

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
        m = 2+0.5*average_read_len # 50% of average input read length
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

    if args.assembly:
        # Preprocessing: rename and reorganize reads
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

        # Align to reference and output all possible alignments (bwa mem -a)
        print "Align reads to reference using bwa mem -a"
        try:
            subprocess.check_call("bwa mem -a -t %s %s assembly/s_p1_p2.fastq 2> /dev/null | samtools sort > assembly/s_p1_p2.bam" % (args.threads, args.reference), shell=True)
        except subprocess.CalledProcessError as e:
            subprocess.check_call("bwa index %s 1>/dev/null 2>&1" % args.reference, shell=True)
            subprocess.check_call("bwa mem -a -t %s %s assembly/s_p1_p2.fastq 2> /dev/null | samtools sort > assembly/s_p1_p2.bam" % (args.threads, args.reference), shell=True)
        subprocess.check_call('samtools index assembly/s_p1_p2.bam', shell=True)
        subprocess.check_call('samtools depth assembly/s_p1_p2.bam > assembly/s_p1_p2.depth', shell=True)

        # TODO: Analyze sequencing depth and target size

        # find alignment start and end points globally
        print "Analyze regions of alignments"
        with open('assembly/s_p1_p2.depth', 'r') as f:
            chr_name = ""
            pos = "0"
            for line in f:
                oldpos = int(pos)
                [chrom, pos, depth] = line.rstrip('\n').split('\t')
                if chr_name != chrom: # new chromosome start
                    start = 0
                    end = 0
                    maxdepth = 0
                chr_name = chrom
                if start == 0:
                    start = int(pos)
                    maxdepth = max(maxdepth, int(depth))
                elif int(pos)-oldpos > 1:
                    end = oldpos
                    if maxdepth > 2 and end-start > average_read_len:
                        chrom2regions[chrom].append([start, end])
                        # print chrom, start, end
                    start = int(pos)
                    end = 0
                    maxdepth = 0
                else:
                    maxdepth = max(maxdepth, int(depth))
                    continue
            end = int(pos)
            if start > 0 and maxdepth > 2 and end-start > average_read_len:
                chrom2regions[chrom].append([start, end])
                #print chrom, start, end

        # Split into regions of 10kb, overlapping by 1000 bp on either side, and
        # divide reads from fastq into separate input files
        print "Split reads into regions of {}bp".format(args.split_size)
        chrom2finalsplit = {}
        for chrom, regions in chrom2regions.iteritems():
            final_split = []
            length = ref_pieces[chrom]
            idx = 0
            pos = args.split_size + regions[idx][0]
            while idx < len(regions):
                if regions[idx][0] >= pos:
                    pos = regions[idx][0] + args.split_size
                region_lb = max(0, pos - (args.split_size + args.split_overlap))
                region_ub = min(length, pos)
                # select reads and write to file
                dirname = 'assembly/%s_%s_%s' % (chrom, region_lb, region_ub)
                overwrite_dir(dirname)
                overwrite_dir('%s/assembly' % dirname)
                subprocess.check_call('samtools view -F4 -bh assembly/s_p1_p2.bam %s:%s-%s | samtools fastq - 1> %s/assembly/s_p1_p2.fastq 2>/dev/null' % (chrom, region_lb, region_ub, dirname), shell=True)
                if file_len('%s/assembly/s_p1_p2.fastq' % dirname) >= 400:
                    final_split.append([region_lb, region_ub])
                else:
                    subprocess.check_call('rm -rf %s' % dirname, shell=True)
                while idx < len(regions) and regions[idx][1] < pos:
                    idx += 1
                pos += args.split_size
            chrom2finalsplit[chrom] = final_split

        # TODO: If a region has less than 1kb covered, merge with previous/next region
        # depending on alignments

        print "\rDone!" + ' ' * 40
        sys.stdout.flush()

        if args.split_only:
            # split_only flag means we only do preprocessing, no assembly
            sys.exit()

        # Run savage-lc on each region
        print "Run savage-lc per split region"
        original_fastq = cwd + "/assembly/s_p1_p2.fastq"
        settings = [args, base_path, s_seq_count, p_seq_count, original_fastq,
                        min_overlap_len_EC, average_read_len, max_tip_len]
        if args.pool_size == 1:
            for chrom, final_split in chrom2finalsplit.iteritems():
                for region in final_split:
                    run_savage_lc(settings, chrom, region)
        else:
            pool = Pool(args.pool_size)
            for chrom, final_split in chrom2finalsplit.iteritems():
                pool.map(partial(run_savage_lc, settings, chrom), final_split)
            pool.close()
            pool.join()

            # for [region_lb, region_ub] in final_split:
                # dirname = 'assembly/%s_%s_%s' % (chrom, region_lb, region_ub)
                # if file_len('%s/assembly/s_p1_p2.fastq' % dirname) < 100:
                #     continue
                # os.chdir(dirname)
                # # run savage-lc
                # savage_command = "%s/savage-lc.py --no_preprocessing" % base_path
                # savage_command += " -s assembly/s_p1_p2.fastq"
                # savage_command += " --hap_cov %s" % args.hap_cov
                # savage_command += " --insert_size %s" % args.insert_size
                # savage_command += " --stddev %s" % args.stddev
                # savage_command += " --original_SE_count %s" % s_seq_count
                # savage_command += " --original_PE_count %s" % p_seq_count
                # savage_command += " --original_fastq %s" % original_fastq
                # savage_command += " -m %s" % args.min_overlap_len
                # savage_command += " -m_EC %s" % min_overlap_len_EC
                # savage_command += " -t %s" % args.threads
                # savage_command += " --mismatch_rate %s" % args.merge_contigs
                # savage_command += " --min_clique_size %s" % args.min_clique_size
                # savage_command += " --sfo_err %s" % args.sfo_err
                # savage_command += " --average_read_len %s" % average_read_len
                # savage_command += " --max_tip_len %s" % max_tip_len
                # #savage_command += " --min_evidence %s" % args.min_evidence
                # if args.diploid:
                #     if args.diploid_overlap_len:
                #         diploid_overlap_len = args.diploid_overlap_len
                #     else:
                #         diploid_overlap_len = args.min_overlap_len
                #     savage_command += " --diploid"
                #     savage_command += " --diploid_contig_len %s" % args.diploid_contig_len
                #     savage_command += " --diploid_overlap_len %s" % diploid_overlap_len
                # if args.count_strains:
                #     savage_command += " --count_strains --ref %s" % args.reference
                # savage_command += " > savage.log 2>&1"
                # subprocess.check_call(savage_command, shell=True)
                # os.chdir('../..')

        # Combine contigs, rename, and update subreads
        new_subreads = open('assembly/subreads.txt', 'w')
        combined_count = 0
        for chrom, final_split in chrom2finalsplit.iteritems():
            for [region_lb, region_ub] in final_split:
                dirname = 'assembly/%s_%s_%s' % (chrom, region_lb, region_ub)
                # add contigs to combined contig file
                if args.diploid and os.path.exists('%s/diploid/singles.fastq' % dirname):
                    # print "diploid contigs from", dirname
                    subprocess.check_call("cat %s/diploid/singles.fastq >> assembly/tmp_contigs.fastq" % dirname, shell=True)
                    subreads = read_subreads("%s/diploid/subreads.txt" % dirname)
                elif os.path.exists('%s/assembly/singles.fastq' % dirname):
                    # print "contigs from", dirname
                    subprocess.check_call("cat %s/assembly/singles.fastq >> assembly/tmp_contigs.fastq" % dirname, shell=True)
                    subreads = read_subreads("%s/assembly/subreads.txt" % dirname)
                else:
                    continue
                # process subreads
                for line in subreads:
                    splitline = line.split('\t')
                    contig_id = int(splitline[0])
                    new_contig_id = contig_id + combined_count
                    splitline[0] = str(new_contig_id)
                    new_subreads.write('\t'.join(splitline))
                combined_count += len(subreads)
        new_subreads.close()
        # rename combined contigs
        subprocess.check_call("%s/scripts/rename_fas.py --in assembly/tmp_contigs.fastq --out assembly/combined_contigs.fastq" % (base_path), shell=True)
        subprocess.check_call("rm assembly/tmp_contigs.fastq", shell=True)
        print "combined contig count:", combined_count
        subprocess.check_call("bwa mem -a -t %s %s assembly/combined_contigs.fastq > assembly/combined_contigs.sam" % (args.threads, args.reference), shell=True)
        #subprocess.check_call("bwa mem -a -t %s %s assembly/combined_contigs.fastq 2> /dev/null | samtools sort > assembly/combined_contigs.bam" % (args.threads, args.reference), shell=True)
        #subprocess.check_call('samtools index assembly/combined_contigs.bam', shell=True)

    # continue assembly on combined contig file
    # OPTION1: compute exact overlaps for neighboring regions

    # OPTION2: run reference-guided savage
    combined_dir = "combined"
    overwrite_dir(combined_dir)
    os.chdir(combined_dir)
    os.mkdir("assembly")
    shutil.copyfile("../assembly/combined_contigs.fastq", "assembly/s_p1_p2.fastq")
    savage_command = "%s/savage-lc.py" % base_path
    savage_command += " -s ../assembly/combined_contigs.fastq"
    savage_command += " --ref %s" % args.reference
    savage_command += " --ref_guided_mode"
    savage_command += " --hap_cov=0"
    savage_command += " --insert_size %s" % args.insert_size
    savage_command += " --stddev %s" % args.stddev
    savage_command += " -m %s" % args.min_overlap_len
    savage_command += " --no_EC"
    savage_command += " --min_clique_size=2"
    savage_command += " -t %s" % args.threads
    savage_command += " --mismatch_rate %s" % args.merge_contigs
    savage_command += " --average_read_len %s" % average_read_len
    savage_command += " --max_tip_len %s" % max_tip_len
    if args.diploid:
        if args.diploid_overlap_len:
            diploid_overlap_len = args.diploid_overlap_len
        else:
            diploid_overlap_len = args.min_overlap_len
        savage_command += " --diploid"
        savage_command += " --diploid_contig_len %s" % args.diploid_contig_len
        savage_command += " --diploid_overlap_len %s" % diploid_overlap_len
        final_contig_file = combined_dir + "/contigs_diploid.fasta"
    else:
        final_contig_file = combined_dir + "/contigs.fasta"
    if args.count_strains:
        savage_command += " --count_strains"
    savage_command += " > savage.log 2>&1"
    subprocess.check_call(savage_command, shell=True)
    os.chdir("..")

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

def read_subreads(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    return lines

def run_savage_lc(settings, chrom, region):
    # print "run_savage_lc({}, {})".format(chrom, region)
    [region_lb, region_ub] = region
    [args, base_path, s_seq_count, p_seq_count, original_fastq,
                min_overlap_len_EC, average_read_len, max_tip_len] = settings
    if args.sfo_threads:
        threads = args.sfo_threads
    else:
        threads = args.threads
    dirname = 'assembly/%s_%s_%s' % (chrom, region_lb, region_ub)
    os.chdir(dirname)
    # run savage-lc
    savage_command = "%s/savage-lc.py --no_preprocessing" % base_path
    savage_command += " -s assembly/s_p1_p2.fastq"
    savage_command += " --hap_cov %s" % args.hap_cov
    savage_command += " --insert_size %s" % args.insert_size
    savage_command += " --stddev %s" % args.stddev
    savage_command += " --original_SE_count %s" % s_seq_count
    savage_command += " --original_PE_count %s" % p_seq_count
    savage_command += " --original_fastq %s" % original_fastq
    savage_command += " -m %s" % args.min_overlap_len
    savage_command += " -m_EC %s" % min_overlap_len_EC
    savage_command += " -t %s" % threads
    savage_command += " --mismatch_rate %s" % args.merge_contigs
    savage_command += " --min_clique_size %s" % args.min_clique_size
    savage_command += " --sfo_err %s" % args.sfo_err
    savage_command += " --average_read_len %s" % average_read_len
    savage_command += " --max_tip_len %s" % max_tip_len
    #savage_command += " --min_evidence %s" % args.min_evidence
    if args.diploid:
        if args.diploid_overlap_len:
            diploid_overlap_len = args.diploid_overlap_len
        else:
            diploid_overlap_len = args.min_overlap_len
        savage_command += " --diploid"
        savage_command += " --diploid_contig_len %s" % args.diploid_contig_len
        savage_command += " --diploid_overlap_len %s" % diploid_overlap_len
    if args.count_strains:
        savage_command += " --count_strains --ref %s" % args.reference
    savage_command += " > savage.log 2>&1"
    try:
        subprocess.check_call(savage_command, shell=True)
    except subprocess.CalledProcessError as e:
        print "%s failed" % dirname
    os.chdir('../..')
    return


if __name__ == '__main__':
    sys.exit(main())
