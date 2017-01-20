#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random


__author__ = "Jasmijn Baaijens"

usage = """

Rename the read identifiers in a fasta/fastq file with positive integers, 
starting from a specified int value (id_start). 

"""

def revcomp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    s = ''.join(bases)
    return s[::-1]

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument('--in', dest='infile', type=str, required=True, 
                        help="input file: overlap results from SFO")
    parser.add_argument('--out', dest='outfile', type=str, required=True, 
                        help="output file: overlap results for SAVAGE")
    parser.add_argument('--id_start', dest='id_start', type=int, default=0, 
                        help="new identifiers will be: id_start, id_start + 1, id_start + 2, etc.")
    parser.add_argument('--revcomp', dest='revcomp', action='store_true', 
                        help="build reverse complementary sequences")
    args = parser.parse_args()
          
    infile = args.infile
    outfile = args.outfile
    
    if infile[-1] == 'a': # fasta format 
        k = 2
    elif infile[-1] == 'q': # fastq format 
        k = 4
    else:
        print """input file must be in fasta format (.fa or .fasta) 
            or fastq format (.fq or .fastq)"""
        
    with open(infile, 'rb') as f1:
        with open(outfile, 'wb') as f2:
            c = 0
            cur_id = args.id_start
            for line in f1:
                c += 1
                if (c % k) == 1: 
                    # ID line
                    f2.write(line[0] + str(cur_id) + '\n')
                    cur_id += 1
                elif (k == 4 and (c % k) == 2) or (k == 2 and (c % k) == 0): 
                    # seq line
                    if args.revcomp:
                        seq = revcomp(line.strip())
                        f2.write(seq + '\n')
                    else:
                        f2.write(line)
                elif k == 4 and (c % k) == 3:
                    # + line
                    f2.write('+\n')
                else:
                    # qual line
                    assert k == 4 and (c % k) == 0
                    if args.revcomp:
                        rev_qual = line.strip()[::-1]
                        f2.write(rev_qual + '\n')
                    else:
                        f2.write(line)
                    
                     
if __name__ == '__main__':
        sys.exit(main())
