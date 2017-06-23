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

Create an overlaps file for viral quasispecies assembly
based on the alignments found by running blast
(single-end reads only)

"""

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument('--in', dest='infile', type=str)
    parser.add_argument('--out', dest='outfile', type=str)
    parser.add_argument('--min_overlap_len', dest='min_overlap_len', type=int, default=0)
    parser.add_argument('--no_revcomp', dest="norevcomp", action='store_true')
    args = parser.parse_args()

    if not (args.infile and args.outfile):
        print "Specify input and output files."
        parser.print_help()
        sys.exit(1)

    with open(args.outfile, 'w') as f1:
        with open(args.infile, 'r') as f2:
            rev_comp_count = 0
            too_short_count = 0
            overlap_count = 0
            for line in f2:
                [qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, qlen, slen] = line.strip('\n').split('\t')
                if qseqid == sseqid:
                    # self-overlap
                    continue
                if int(qstart) > int(qend) or int(sstart) > int(send):
                    rev_comp_count += 1
                    if args.norevcomp:
                        continue
                qori = int(qstart) <= int(qend)
                sori = int(sstart) <= int(send)
                assert qori # otherwise script needs extra cases to find correct pos1
                if not sori:
                    sstart = int(slen)-int(sstart)+1
                    send = int(slen)-int(send)+1
                len1 = length
                if int(len1) < args.min_overlap_len:
                    too_short_count += 1
                    continue
                if int(qstart)-int(sstart) >= 0:
                    id1 = qseqid
                    id2 = sseqid
                    pos1 = str(int(qstart)-int(sstart))
                    ori1 = "+" if qori else "-"
                    ori2 = "+" if sori else "-"
                    if int(pos1) >= int(qlen):
                        print "qstart > 1"
                        print pos1
                        print qlen
                    assert int(pos1) < int(qlen)
                else:
                    id1 = sseqid
                    id2 = qseqid
                    ori1 = "+" if sori else "-"
                    ori2 = "+" if qori else "-"
                    pos1 = str(int(sstart)-int(qstart))
                    if int(pos1) >= int(slen):
                        print "sstart > 1"
                        print pos1
                        print slen
                    assert int(pos1) < int(slen)
                assert pos1 >= 0
                perc = int(round(100*max(float(len1)/float(qlen), float(len1)/float(slen))))
                if perc > 100:
                    perc = 100
                perc1 = str(perc)
                len2 = "-"
                perc2 = "-"
                order = "-"
                pos2 = "-"
#                ori1 = "+"
#                ori2 = "+"
                type1 = "s"
                type2 = "s"
                overlap = '\t'.join([id1, id2, pos1, pos2, order, ori1, ori2, perc1, perc2, len1, len2, type1, type2])
#                if (ori1 == "+" or ori2 == "-"):
#                    continue
#                    print line
#                    print overlap
                f1.write(overlap + '\n')
                overlap_count += 1
            print "reverse complements found: ", rev_comp_count
            print "overlaps shorter than %d bp: %d" %(args.min_overlap_len, too_short_count)
            print "total overlaps found: ", overlap_count

if __name__ == '__main__':
    sys.exit(main())
