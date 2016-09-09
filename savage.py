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
pear_reads directory in current path, containing singles.fastq, paired1.fastq, paired2.fastq 
 
RESULTS:
Directories: stage_a, stage_b, stage_c
Files in current directory: stage_a_contigs.fasta, stage_b_contigs.fasta, stage_c_contigs.fasta

Run savage.py -h for a complete description of optional arguments.

"""

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
    

def main():
    parser = ArgumentParser(description=usage, formatter_class=RawTextHelpFormatter)
    parser.add_argument('--ref', dest='reference', type=str, help='reference genome in fasta format')
    parser.add_argument('--singles', dest='singles', type=str, help='single-end read alignments in SAM format')
    parser.add_argument('--paired', dest='paired', type=str, help='paired-end read alignments in SAM format')
    parser.add_argument('--threshold', dest='threshold', type=float, default=0.97, help='edge threshold for Stage a')
    parser.add_argument('--overlaps', dest='overlaps', type=str, help='skip overlap computations by using given overlaps file; make sure to enter the full path!!')
    parser.add_argument('--no_stage_a', dest='stage_a', action='store_false', help='skip Stage a (initial contig formation); \n--> use this option together with --contigs')
    parser.add_argument('--no_stage_b', dest='stage_b', action='store_false', help='skip Stage b (extending initial contigs); \n--> this automatically skips Stage c as well')
    parser.add_argument('--no_stage_c', dest='stage_c', action='store_false', help='skip Stage c (merging maximized contigs into master strains)')
    parser.add_argument('--contigs', dest='contigs', type=str, help='contigs fastq file resulting from Stage a; \n--> use this option together with --no_stage_a')
    args = parser.parse_args()
    
    base_path = sys.path[0]
    print "base path: %s" % base_path
    FNULL = open(os.devnull, 'w')
    
    if args.overlaps and (args.reference or args.singles or args.paired):
        print "Overlaps file given, so reference / singles sam / paired sam are ignored"
        
    if args.overlaps:
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
    
    print "Edge threshold is set to %.3f" % args.threshold
    
    if preprocessing:
        overlaps = "../original_overlaps.txt"
        if denovo:
            print "Preprocessing: get de novo overlaps using SFO"
            # prepare fasta
            subprocess.check_call("cat pear_reads/singles.fastq pear_reads/paired1.fastq pear_reads/paired2.fastq > s_p1_p2.fastq", shell=True)
            subprocess.check_call("python %s/scripts/fastq2fasta.py s_p1_p2.fastq s_p1_p2.fasta" % base_path, shell=True)
            singles_count = int(file_len('pear_reads/singles.fastq')/4.0)
            paired_count = int(file_len('pear_reads/paired1.fastq')/4.0)
            assert paired_count == int(file_len('pear_reads/paired2.fastq')/4.0)
            # run SFO
            print "Running SFO..."
            subprocess.check_call("%s/sfo_2011_5/builder s_p1_p2.fasta" % base_path, shell=True)
            subprocess.check_call("%s/sfo_2011_5/sfoverlap --parallel 0 --indels -e 50 -t 167 s_p1_p2.fasta | %s/sfo_2011_5/maxoverlaps > sfoverlaps.out" % (base_path, base_path), shell=True)
            # run postprocessing scripts
            print "Processing SFO output..."
            subprocess.check_call("python %s/scripts/sfo2overlaps_script1.py %d %d sfoverlaps.out" % (base_path, singles_count, paired_count), stdout=FNULL, shell=True)
            subprocess.check_call("python %s/scripts/sfo2overlaps_script2.py %d %d" % (base_path, singles_count, paired_count), stdout=FNULL, shell=True)
            # combine output and remove intermediate files
            print "Finishing and cleaning up..."
            subprocess.check_call("cat *formatted9 > original_overlaps.txt", shell=True)
            subprocess.check_call("rm *_sorted", shell=True)
            subprocess.check_call("rm *formatted9", shell=True)
            subprocess.check_call("rm s_p1_p2.fast*", shell=True)
            print "Overlaps are ready!"
        else:
            print "Preprocessing: get induced overlaps from alignments"
            # Induce overlaps from alignment
            if args.singles and args.paired:
                print
                subprocess.check_call("python %s/scripts/sam2overlaps.py --sam_p %s --sam_s %s --ref %s --min_overlap_len 150 --out original_overlaps.txt" %(base_path, args.paired, args.singles, args.reference), shell=True)
            elif args.singles:
                subprocess.check_call("python %s/scripts/sam2overlaps.py --sam_s %s --ref %s --min_overlap_len 150 --out original_overlaps.txt" %(base_path, args.singles, args.reference), shell=True)
            elif args.paired:
                subprocess.check_call("python %s/scripts/sam2overlaps.py --sam_p %s --ref %s --min_overlap_len 150 --out original_overlaps.txt" %(base_path, args.paired, args.reference), shell=True)
    else:
        print "Using overlaps from %s" % args.overlaps
        overlaps = args.overlaps

    if args.stage_a:
        print "**************"
        print "SAVAGE Stage a"
        # Run SAVAGE Stage a: error correction and initial contig formation
        subprocess.call(['mkdir', 'stage_a'], stdout=FNULL, stderr=FNULL)
        os.chdir('stage_a')
        subprocess.check_call("python %s/scripts/pipeline_stage_a.py --fastq ../pear_reads --overlaps %s --merging --cliques --error_correction --edge_threshold %f --min_overlap_perc 0 --fno 3 --min_clique_size 4 --transitive_edges 2" %(base_path, overlaps, args.threshold), stdout=FNULL, shell=True)
        os.chdir('..')
        subprocess.check_call("python %s/scripts/fastq2fasta.py stage_a/singles.fastq contigs_stage_a.fasta" % base_path, shell=True)
        print "Done!"  
    else:
        if args.contigs:
            subprocess.call(['mkdir', 'stage_a'], stdout=FNULL, stderr=FNULL)
            subprocess.call(['cp', args.contigs, 'stage_a/singles.fastq'], stdout=FNULL, stderr=FNULL)
            subprocess.check_call("python %s/scripts/fastq2fasta.py stage_a/singles.fastq contigs_stage_a.fasta" % base_path, shell=True)
        print "Stage a skipped"

    if args.stage_b and args.stage_a:
        print "**************"
        print "SAVAGE Stage b"
        # Run SAVAGE Stage b: build maximized contigs
        subprocess.check_call("makeblastdb -in contigs_stage_a.fasta -dbtype nucl -out contigs_db", shell=True)
        subprocess.check_call("blastn -db contigs_db -query contigs_stage_a.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen' -out blastout_contigs.tsv -perc_identity 98", shell=True)
        subprocess.check_call("python %s/scripts/blast2overlaps.py --in blastout_contigs.tsv --out contig_overlaps.txt --min_overlap_len 50" % base_path, shell=True)
        subprocess.call(['mkdir', 'stage_b'], stdout=FNULL, stderr=FNULL)
        subprocess.check_call("cp stage_a/subreads.txt stage_b/subreads.txt", shell=True)
        os.chdir('stage_b')
        subprocess.check_call("python %s/scripts/pipeline_stages_b_c.py --fastq ../stage_a/singles.fastq --overlaps ../contig_overlaps.txt --transitive_edges 1 --merge_contigs 0 --min_qual 0.999 --use_subreads" % base_path, stdout=FNULL, shell=True)
        os.chdir('..')
        subprocess.check_call("python %s/scripts/fastq2fasta.py stage_b/singles.fastq contigs_stage_b.fasta" % base_path, shell=True)
        print "Done!"
    elif args.stage_b: # (but not stage a)
        print "**************"
        print "SAVAGE Stage b"
        # Run SAVAGE Stage b: build maximized contigs
        subprocess.check_call("makeblastdb -in contigs_stage_a.fasta -dbtype nucl -out contigs_db", shell=True)
        subprocess.check_call("blastn -db contigs_db -query contigs_stage_a.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen' -out blastout_contigs.tsv -perc_identity 98", shell=True)
        subprocess.check_call("python %s/scripts/blast2overlaps.py --in blastout_contigs.tsv --out contig_overlaps.txt --min_overlap_len 50" % base_path, shell=True)
        subprocess.call(['mkdir', 'stage_b'], stdout=FNULL, stderr=FNULL)
        subprocess.check_call("cp stage_a/subreads.txt stage_b/subreads.txt", shell=True)
        os.chdir('stage_b')
        subprocess.check_call("python %s/scripts/pipeline_stages_b_c.py --fastq ../stage_a/singles.fastq --overlaps ../contig_overlaps.txt --transitive_edges 1 --merge_contigs 0 --min_qual 0.999" % base_path, stdout=FNULL, shell=True) # note: not using stage a subreads
        os.chdir('..')
        subprocess.check_call("python %s/scripts/fastq2fasta.py stage_b/singles.fastq contigs_stage_b.fasta" % base_path, shell=True)
        print "Done!"
        print "Note: stage a was skipped, so stage b did not use the subread information from stage a. For frequency estimation of stage b or stage c contigs, make sure to provide the subreads files from stage a (see the README)."
    else:
        print "Stage b skipped"


    if args.stage_c:
        print "**************"
        print "SAVAGE Stage c"
        # Run SAVAGE Stage c: build master strains
        subprocess.check_call("makeblastdb -in contigs_stage_b.fasta -dbtype nucl -out contigs_db", shell=True)
        subprocess.check_call("blastn -db contigs_db -query contigs_stage_b.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen' -out blastout_contigs.tsv -perc_identity 98", shell=True)
        subprocess.check_call("python %s/scripts/blast2overlaps.py --in blastout_contigs.tsv --out contig_overlaps.txt --min_overlap_len 50" % base_path, shell=True)
        subprocess.call(['mkdir', 'stage_c'], stdout=FNULL, stderr=FNULL)
        subprocess.check_call("cp stage_b/subreads.txt stage_c/subreads.txt", shell=True)
        os.chdir('stage_c')
        subprocess.check_call("python %s/scripts/pipeline_stages_b_c.py --fastq ../stage_b/singles.fastq --overlaps ../contig_overlaps.txt --transitive_edges 1 --merge_contigs 0.01 --min_qual 0.999 --use_subreads" % base_path, stdout=FNULL, shell=True)
        os.chdir('..')
        subprocess.check_call("python %s/scripts/fastq2fasta.py stage_c/singles.fastq contigs_stage_c.fasta" % base_path, shell=True)
        subprocess.call("rm blastout* contigs_db* contig_overlaps.txt", shell=True)
        print "Done!"  
    else:
        print "Stage c skipped"
      
if __name__ == '__main__':
    sys.exit(main())
