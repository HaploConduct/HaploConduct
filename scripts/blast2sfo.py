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
                if qseqid == sseqid:
                    # self-overlap
                    continue
                qori = int(qstart) <= int(qend)
                sori = int(sstart) <= int(send)
                assert qori # otherwise script needs extra cases to find correct pos1
                if int(length) < args.min_overlap_len:
                    too_short_count += 1
                    continue

                idA = qseqid
                idB = sseqid
                ori = 'N' if sori else 'I'
                OLA = length
                OLB = OLA
                if sori:
                    OHA = int(qstart)-int(sstart)
                    OHB = int(slen)-int(sstart)-(int(qlen)-int(qstart))
                else:
                    OHA = int(qstart)-(int(slen)-int(sstart)+1)
                    OHB = int(sstart)-(int(qlen)-int(qstart)+1)
#                if int(idA) > int(idB):
                if idA > idB:
                    # swap order such that id1 < id2
                    idA, idB = idB, idA
                    if ori == 'N':
                        OHA *= -1
                        OHB *= -1
                    else:
                        # swap orientations such that id1 sequence is forward
                        OHA, OHB = OHB, OHA
                sfo_line = '\t'.join([idA, idB, ori, str(OHA), str(OHB), OLA, OLB, mismatch]) + '\n'
                f1.write(sfo_line)
                overlap_count += 1
            print "overlaps shorter than %d bp: %d" %(args.min_overlap_len, too_short_count)
            print "total overlaps found: ", overlap_count

if __name__ == '__main__':
    sys.exit(main())
