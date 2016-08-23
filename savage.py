#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random
import subprocess
from time import clock
from sets import Set


__author__ = "Jasmijn Baaijens"

usage = """%prog [options]

SAVAGE: Strain Aware VirAl GEnome assembly

# REQUIRED DIRECTORY STRUCTURE:
# pear_reads directory in current path, containing singles.fastq, paired1.fastq, paired2.fastq 
# 
# Resulting directory structure:
# Directories: stage_a, stage_b, stage_c
# Files in current directory: stage_a_contigs.fasta, stage_b_contigs.fasta, stage_c_contigs.fasta

Run savage.py -h for all arguments.

"""

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument('--ref', dest='reference', type=str)
    parser.add_argument('--singles', dest='singles', type=str)
    parser.add_argument('--paired', dest='paired', type=str)
    parser.add_argument('--no_stage_a', dest='stage_a', action='store_false')
    parser.add_argument('--no_stage_b', dest='stage_b', action='store_false')
    parser.add_argument('--no_stage_c', dest='stage_c', action='store_false')
    parser.add_argument('--edge_threshold', dest='threshold', default=0.97)
    args = parser.parse_args()
    
    base_path = sys.path[0]
    print "base path: %s" % base_path
    FNULL = open(os.devnull, 'w')
    
    if args.overlaps and (args.reference or args.singles or args.paired):
        print "Overlaps file given, so reference / singles sam / paired sam are ignored"
        
    if args.reference and (args.singles or args.paired):
        preprocessing = True
    elif args.singles or args.paired:
        print usage
        print "Provide a reference fasta file"
        print 
        return -1
    else:
        print usage
        print "Provide samfile(s)"
        print 
        return -1
    
    print "Edge threshold is set to %f" % args.threshold
    
    if preprocessing:
        print "Preprocessing: get induced overlaps"
        # Induce overlaps from alignment
        if args.singles and args.paired:
            print
            subprocess.check_call("python %s/scripts/sam2overlaps.py --sam_p %s --sam_s %s --ref %s --min_overlap_len 150 --out original_overlaps.txt" %(base_path, args.paired, args.singles, args.reference), shell=True)
        elif args.singles:
            subprocess.check_call("python %s/scripts/sam2overlaps.py --sam_s %s --ref %s --min_overlap_len 150 --out original_overlaps.txt" %(base_path, args.singles, args.reference), shell=True)
        elif args.paired:
            subprocess.check_call("python %s/scripts/sam2overlaps.py --sam_p %s --ref %s --min_overlap_len 150 --out original_overlaps.txt" %(base_path, args.paired, args.reference), shell=True)
    else:
        print "Preprocessing: get de novo overlaps"
        print "This option is not yet supported but will appear soon!"
        return -1

    if args.stage_a:
        print "**************"
        print "SAVAGE Stage a"
        # Run SAVAGE Stage a: error correction and initial contig formation
        subprocess.call(['mkdir', 'stage_a'], stdout=FNULL, stderr=FNULL)
        os.chdir('stage_a')
        subprocess.check_call("python %s/scripts/pipeline_stage_a.py --fastq ../pear_reads --overlaps ../original_overlaps.txt --merging --cliques --error_correction --edge_threshold %f --min_overlap_perc 0 --fno 3 --min_clique_size 4 --transitive_edges 2" %(base_path, args.threshold), stdout=FNULL, shell=True)
        os.chdir('..')
        subprocess.check_call("python %s/scripts/fastq2fasta.py stage_a/singles.fastq contigs_stage_a.fasta" % base_path, shell=True)
        print "Done!"  
    else:
        print "Stage a skipped"

    if args.stage_b:
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
