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

Convert savage overlaps to SFO format

savage format: [id1, id2, pos1, pos2, order, ori1, ori2, perc1, perc2, len1, len2, type1, type2]
SFO format: [idA, idB, ori, OHA, OHB, OLA, OLB, K]

"""

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument('--overlaps_in', dest='infile', type=str)
    parser.add_argument('--overlaps_out', dest='outfile', type=str)
    parser.add_argument('--fasta', dest='fasta', type=str)
    parser.add_argument('-m', dest='min_overlap_len', default=0, type=int)
    args = parser.parse_args()

    if not (args.infile and args.outfile):
        print "Specify input and output files."
        parser.print_help()
        sys.exit(1)

    id2len = {}
    with open(args.fasta, 'r') as f:
        c = 0
        for line in f:
            c += 1
            if (c % 2) == 1:
                # ID line
                ID = line.lstrip('>').rstrip('\n')
            else:
                # seq line
                seqlen = len(line.strip('\n'))
                id2len[ID] = seqlen
        assert int(c/2) == len(id2len)

    with open(args.outfile, 'w') as f1:
        with open(args.infile, 'r') as f2:
            too_short_count = 0
            overlap_count = 0
            for line in f2:
                [id1, id2, pos1, pos2, order, ori1, ori2, perc1, perc2, len1, len2, type1, type2] = line.strip('\n').split('\t')
                assert type1 == type2 == 's' # only allow s-s overlaps
                length = int(len1)
                if int(length) < args.min_overlap_len:
                    too_short_count += 1
                    continue
                # translate to sfo format
                idA = id1
                idB = id2
                if [ori1, ori2] == ['+', '-'] or [ori1, ori2] == ['-', '+']:
                    ori = 'I'
                else:
                    ori = 'N'
                OHA = int(pos1)
                OLA = len1
                if id2len[id2] == int(len1):
                    OHB = int(pos1) + int(len1) - id2len[id1]
                else:
                    OHB = id2len[id2] - int(len1)
                OLB = len1

#                if int(idA) > int(idB):
                if idA > idB:
                    # swap order such that id1 < id2
                    idA, idB = idB, idA
                    ori1, ori2 = ori2, ori1
                    OHA *= -1
                    OHB *= -1

                if ori1 == "-":
                    # swap orientations such that id1 sequence is forward
                    OHA, OHB = -OHB, -OHA

                sfo_line = '\t'.join([idA, idB, ori, str(OHA), str(OHB), OLA, OLB]) + '\t.\n'
                f1.write(sfo_line)
                overlap_count += 1
            print "overlaps shorter than %d bp: %d" %(args.min_overlap_len, too_short_count)
            print "total overlaps found: ", overlap_count

if __name__ == '__main__':
    sys.exit(main())
