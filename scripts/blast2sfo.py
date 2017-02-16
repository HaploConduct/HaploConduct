#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random
import subprocess
from time import clock


__author__ = "Jasmijn Baaijens"
__license__ = "GPL"

usage = """

Convert blast overlaps to SFO format

"""

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument('--in', dest='infile', type=str)
    parser.add_argument('--out', dest='outfile', type=str)
    parser.add_argument('-m', dest='min_overlap_len', required=True, type=int)
    args = parser.parse_args()

    if not (args.infile and args.outfile):
        print "Specify input and output files."
        parser.print_help()
        sys.exit(1)

    with open(args.outfile, 'w') as f1:
        with open(args.infile, 'r') as f2:
            too_short_count = 0
            overlap_count = 0
            for line in f2:
                [qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, qlen, slen] = line.strip('\n').split('\t')
                qori = int(qstart) <= int(qend)
                sori = int(sstart) <= int(send)
                assert qori # otherwise script needs extra cases to find correct pos1
                if int(length) < args.min_overlap_len:
                    too_short_count += 1
                    continue
                if int(qseqid) <= int(sseqid):
                    idA = qseqid
                    idB = sseqid
                    ori = 'N' if sori else 'I'
                    OHA = int(qstart)-int(sstart)
                    OLA = length
                    OHB = int(slen)-int(send)-(int(qlen)-int(qend))
                    OLB = length
                else:
                    idA = sseqid
                    idB = qseqid
                    ori = 'N' if sori else 'I'
                    OHA = int(sstart)-int(qstart)
                    OLA = length
                    OHB = int(qlen)-int(qend)-(int(slen)-int(send))
                    OLB = length
                sfo_line = '\t'.join([idA, idB, ori, str(OHA), str(OHB), OLA, OLB]) + '\t0\n'
                f1.write(sfo_line)
                overlap_count += 1
            print "overlaps shorter than %d bp: %d" %(args.min_overlap_len, too_short_count)
            print "total overlaps found: ", overlap_count

if __name__ == '__main__':
    sys.exit(main())
