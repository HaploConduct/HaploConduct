#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import os
import sys
import random
import subprocess
from time import clock
from sets import Set


__author__ = "Jasmijn Baaijens"
__license__ = "GPL"

usage = """

SAVAGE: Strain Aware VirAl GEnome assembly

REQUIRED DIRECTORY STRUCTURE:
input_fas directory in current path, containing singles.fastq, paired1.fastq,
paired2.fastq

RESULTS:
Directories: stage_a, stage_b, stage_c
Files in current directory: stage_a_contigs.fasta, stage_b_contigs.fasta,
stage_c_contigs.fasta

Run savage.py -h for a complete description of optional arguments.

"""

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


def main():
    parser = ArgumentParser(description=usage, formatter_class=RawTextHelpFormatter)
    parser.add_argument('--ref', dest='reference', type=str, help='reference genome in fasta format')
    parser.add_argument('--singles', dest='singles', type=str, help='single-end read alignments in SAM format')
    parser.add_argument('--paired', dest='paired', type=str, help='paired-end read alignments in SAM format')
#    parser.add_argument('--threshold', dest='threshold', type=float, default=0.97, help='edge threshold for Stage a')
    parser.add_argument('--overlaps', dest='overlaps', type=str, help='skip overlap computations by using given overlaps file; make sure to enter the full path!!')
    parser.add_argument('--contigs', dest='contigs', type=str, help='contigs fastq file resulting from Stage a; \n--> use this option together with --no_stage_a')
    parser.add_argument('--no_use_subreads', dest='use_subreads', action='store_false', help='use subread info from previous stage; default=false')
    parser.add_argument('--no_stage_a', dest='stage_a', action='store_false', help='skip Stage a (initial contig formation)')
    parser.add_argument('--no_stage_b', dest='stage_b', action='store_false', help='skip Stage b (extending initial contigs)')
    parser.add_argument('--no_stage_c', dest='stage_c', action='store_false', help='skip Stage c (merging maximized contigs into master strains)')
    parser.add_argument('--preprocessing', dest='preprocessing', action='store_true', help='run preprocessing procedure (i.e. overlaps algorithm)')
    parser.add_argument('--sfo_mm', dest='sfo_mm', type=int, default=50, help='input parameter -e=M for sfo: maximal mismatch rate 1/M')
    parser.add_argument('--min_overlap_len', dest='min_overlap_len', type=int, default=150)
    parser.add_argument('--merge_contigs', dest='merge_contigs', type=float, default=0.01, help='specify maximal distance between contigs for merging into master strains (stage c)')
    parser.add_argument('--keep_branches', dest='remove_branches', action='store_false', help='disable merging along branches by removing them from the graph (stage b & c)')
    parser.add_argument('--num_threads', dest='threads', type=int, default=1)
    args = parser.parse_args()

    base_path = sys.path[0]
    print "base path: %s" % base_path
    FNULL = open(os.devnull, 'w')
    remove_branches = 'true' if args.remove_branches else 'false'

    if args.overlaps and (args.reference or args.singles or args.paired):
        print "Overlaps file given, so reference / singles sam / paired sam are ignored"

    if args.stage_a and args.stage_c and not args.stage_b:
        print """Options specified suggest running stages a and c, but skipping stage b.
                 If you really want to do this, then run stage a and c separately."""
        print
        return -1

    if args.overlaps or (not args.stage_a and not args.preprocessing):
        preprocessing = False
    elif args.reference and (args.singles or args.paired):
        preprocessing = True
        denovo = False
    elif args.singles or args.paired:
        print usage
        print "Provide a reference fasta file"
        print
        return -1
    elif args.reference:
        print usage
        print "Provide samfile(s)"
        print
        return -1
    else:
        preprocessing = True
        denovo = True

    if preprocessing:
        overlaps = "../original_overlaps.txt"
        if denovo:
            print "Preprocessing: get de novo overlaps using SFO"
            # prepare fasta
            subprocess.check_call("cat input_fas/singles.fastq input_fas/paired1.fastq input_fas/paired2.fastq > s_p1_p2.fastq", shell=True)
            subprocess.check_call("python %s/scripts/fastq2fasta.py s_p1_p2.fastq s_p1_p2.fasta" % base_path, shell=True)
            singles_count = int(file_len('input_fas/singles.fastq')/4.0)
            paired_count = int(file_len('input_fas/paired1.fastq')/4.0)
            assert paired_count == int(file_len('input_fas/paired2.fastq')/4.0)
            # run SFO
            print "Running SFO..."
            subprocess.check_call("%s/sfo_2011_5/builder s_p1_p2.fasta" % base_path, shell=True)
            if paired_count > 0:
                sfo_len = round(args.min_overlap_len / 2)
            else:
                sfo_len = args.min_overlap_len
            subprocess.check_call("%s/sfo_2011_5/sfoverlap --parallel %d --indels -e %d -t %d s_p1_p2.fasta | %s/sfo_2011_5/maxoverlaps > sfoverlaps.out" % (base_path, args.threads, args.sfo_mm, sfo_len, base_path), shell=True)
            # run postprocessing scripts
            print "Processing SFO output..."
            subprocess.check_call("python %s/scripts/sfo2overlaps.py --in sfoverlaps.out --out original_overlaps.txt --num_singles %d --num_pairs %d" % (base_path, singles_count, paired_count), shell=True)
            print "Overlaps are ready!"
        else:
            print "Preprocessing: get induced overlaps from alignments"
            # Induce overlaps from alignment
            if args.singles and args.paired:
                print
                subprocess.check_call("python %s/scripts/sam2overlaps.py --sam_p %s --sam_s %s --ref %s --min_overlap_len %d --out original_overlaps.txt" %(base_path, args.paired, args.singles, args.reference, args.min_overlap_len/2), shell=True)
            elif args.singles:
                subprocess.check_call("python %s/scripts/sam2overlaps.py --sam_s %s --ref %s --min_overlap_len %d --out original_overlaps.txt" %(base_path, args.singles, args.reference, args.min_overlap_len), shell=True)
            elif args.paired:
                subprocess.check_call("python %s/scripts/sam2overlaps.py --sam_p %s --ref %s --min_overlap_len %d --out original_overlaps.txt" %(base_path, args.paired, args.reference, args.min_overlap_len/2), shell=True)
    elif args.stage_a:
        print "Using overlaps from %s" % args.overlaps
        overlaps = args.overlaps

    # create directories
    subprocess.call(['mkdir', 'stage_a'], stdout=FNULL, stderr=FNULL)
    subprocess.call(['mkdir', 'stage_b'], stdout=FNULL, stderr=FNULL)
    subprocess.call(['mkdir', 'stage_c'], stdout=FNULL, stderr=FNULL)


    # Run SAVAGE Stage a: error correction and initial contig formation
    if args.stage_a:
        print "**************"
        print "SAVAGE Stage a"
        sys.stdout.flush()
        os.chdir('stage_a')
        subprocess.check_call("python %s/scripts/pipeline_per_stage.py --stage a --fastq ../input_fas --overlaps %s --min_overlap_len %d --num_threads %d --remove_branches %s" %(base_path, overlaps, args.min_overlap_len, args.threads, remove_branches), shell=True)
        os.chdir('..')
        subprocess.check_call("python %s/scripts/fastq2fasta.py stage_a/singles.fastq contigs_stage_a.fasta" % base_path, shell=True)
        print "Done!"
    else:
        print "Stage a skipped"

    # Run SAVAGE Stage b: build maximized contigs
    if args.stage_b:
        print "**************"
        print "SAVAGE Stage b"
        # prepare input files
        overlaps = ""
        if not args.stage_a:
            if args.contigs:
                subprocess.call(['cp', args.contigs, 'stage_a/singles.fastq'], stdout=FNULL, stderr=FNULL)
                subprocess.call(['touch', 'stage_a/paired1.fastq'], stdout=FNULL, stderr=FNULL)
                subprocess.call(['touch', 'stage_a/paired2.fastq'], stdout=FNULL, stderr=FNULL)
                subprocess.check_call("python %s/scripts/fastq2fasta.py stage_a/singles.fastq contigs_stage_a.fasta" % base_path, shell=True)
            if args.overlaps:
                print "Using overlaps from %s" % args.overlaps
                overlaps = args.overlaps
#        # consider paired-end reads as singles
#        subprocess.call('cat stage_a/singles.fastq stage_a/paired1.fastq stage_a/paired2.fastq > stage_a/combined_reads.fastq', shell=True)
#        subprocess.check_call( # rename the reads in the combined fastq
#            'python rename_fas.py --in stage_a/combined_reads.fastq --out stage_b/singles.fastq', shell=True)
#        subprocess.check_call("python %s/scripts/fastq2fasta.py stage_b/singles.fastq contigs_stage_a_combined.fasta" % base_path, shell=True)
        subprocess.call(['cp', 'stage_a/singles.fastq', 'stage_b/singles.fastq'], stdout=FNULL, stderr=FNULL)
        if overlaps == "":
            subprocess.check_call("makeblastdb -in contigs_stage_a.fasta -dbtype nucl -out contigs_db", shell=True, stdout=FNULL, stderr=FNULL)
            subprocess.check_call("blastn -db contigs_db -query contigs_stage_a.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen' -out blastout_contigs.tsv -perc_identity 98", shell=True)
            subprocess.check_call("python %s/scripts/blast2overlaps.py --in blastout_contigs.tsv --out contig_overlaps.txt --min_overlap_len %d" % (base_path, args.min_overlap_len), shell=True)
            overlaps = "../contig_overlaps.txt"
        #
        sys.stdout.flush()
        # run SAVAGE
        os.chdir('stage_b')
        if args.use_subreads:
            subprocess.check_call("cp ../stage_a/subreads.txt subreads.txt", shell=True)
            subprocess.check_call("python %s/scripts/pipeline_per_stage.py --stage b --fastq ../stage_b --overlaps %s --use_subreads --min_overlap_len %d --num_threads %d --remove_branches %s" % (base_path, overlaps, args.min_overlap_len, args.threads, remove_branches), shell=True)
        else:
            subprocess.check_call("python %s/scripts/pipeline_per_stage.py --stage b --fastq ../stage_b --overlaps %s --min_overlap_len %d --num_threads %d --remove_branches %s" % (base_path, overlaps, args.min_overlap_len, args.threads, remove_branches), shell=True) # note: not using stage a subreads
        os.chdir('..')
        subprocess.check_call("python %s/scripts/fastq2fasta.py stage_b/singles.fastq contigs_stage_b.fasta" % base_path, shell=True)
        print "Done!"
#        if not args.stage_a:
#            print "Note: stage a was skipped, so stage b did not use the subread information from stage a. For frequency estimation of stage b or stage c contigs, make sure to provide the subreads files from stage a (see the README)."
    else:
        print "Stage b skipped"

    # Run SAVAGE Stage c: build master strains
    if args.stage_c:
        print "**************"
        print "SAVAGE Stage c"
        # prepare input files
        overlaps = ""
        if not args.stage_b:
            if args.contigs:
                subprocess.call(['cp', args.contigs, 'stage_b/singles.fastq'], stdout=FNULL, stderr=FNULL)
                subprocess.check_call("python %s/scripts/fastq2fasta.py stage_b/singles.fastq contigs_stage_b.fasta" % base_path, shell=True)
            if args.overlaps:
                print "Using overlaps from %s" % args.overlaps
                overlaps = args.overlaps
        subprocess.call(['cp', 'stage_b/singles.fastq', 'stage_c/singles.fastq'], stdout=FNULL, stderr=FNULL)
        if overlaps == "":
            pident = 100*(1-args.merge_contigs)
            subprocess.check_call("makeblastdb -in contigs_stage_b.fasta -dbtype nucl -out contigs_db", shell=True, stdout=FNULL, stderr=FNULL)
            subprocess.check_call("blastn -db contigs_db -query contigs_stage_b.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen' -out blastout_contigs.tsv -perc_identity %s" % pident, shell=True)
            subprocess.check_call("python %s/scripts/blast2overlaps.py --in blastout_contigs.tsv --out contig_overlaps.txt --min_overlap_len %d" % (base_path, args.min_overlap_len), shell=True)
            overlaps = "../contig_overlaps.txt"
        #
        sys.stdout.flush()
        # run SAVAGE
        os.chdir('stage_c')
        if args.use_subreads:
            subprocess.check_call("cp ../stage_b/subreads.txt subreads.txt", shell=True)
            subprocess.check_call("python %s/scripts/pipeline_per_stage.py --fastq ../stage_c --overlaps %s --merge_contigs %f --stage c --min_overlap_len %d --use_subreads --num_threads %d --remove_branches %s" % (base_path, overlaps, args.merge_contigs, args.min_overlap_len, args.threads, remove_branches), shell=True)
        else:
            subprocess.check_call("python %s/scripts/pipeline_per_stage.py --fastq ../stage_c --overlaps %s --merge_contigs %f --stage c --min_overlap_len %d --num_threads %d --remove_branches %s" % (base_path, overlaps, args.merge_contigs, args.min_overlap_len, args.threads, remove_branches), shell=True)
        os.chdir('..')
        subprocess.check_call("python %s/scripts/fastq2fasta.py stage_c/singles.fastq contigs_stage_c.fasta" % base_path, shell=True)
        subprocess.call("rm blastout* contigs_db* contig_overlaps.txt", shell=True)
        print "Done!"
    else:
        print "Stage c skipped"

if __name__ == '__main__':
    sys.exit(main())
