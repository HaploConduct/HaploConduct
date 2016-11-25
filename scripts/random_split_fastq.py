#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random
import math


__author__ = "Jasmijn Baaijens"

usage = """%prog [options]

Distribute fastq uniformly over a specified number of files.

"""

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument("--input", dest="input", type=str, required=True, help="input fastq filename")
    parser.add_argument("--input2", dest="input2", type=str, required=False, help="input fastq filename (2)")
    parser.add_argument("--output", dest="output", type=str, required=True, help="prefix for output files")
    parser.add_argument("--split_num", dest="split_num", type=int, required=True, help="number of desired output files")
    parser.add_argument('--illumina_IDs', dest='id_sep', type=bool, default=False)
    args = parser.parse_args()
    fastq = args.input

    # read input file and randomly assign output files to read IDs
    ID_to_set = {}
    with open(fastq, 'r') as f:
        i = 0
        for line in f:
            if i%4 == 0:
                illumina_ID = line.split()[0]
                if args.id_sep:
                    cur_ID = int(illumina_ID.split('.')[1])
                else:
                    cur_ID = int(illumina_ID.strip('@'))
#                if i > 0 and cur_ID <= prev_ID:
#                    break
                prev_ID = cur_ID
                rand_set = random.randint(0, args.split_num-1)
                ID_to_set[cur_ID] = rand_set
            i += 1

#    print "ID_to_set size: ", len(ID_to_set)
#    print "total lines read: ", int(i/4)

    num_files = args.split_num
    totalcount = len(ID_to_set)
    reads_per_file = math.floor(totalcount / num_files)
    remaining = totalcount - reads_per_file * num_files
    count_per_file = [reads_per_file + 1 if j < remaining else reads_per_file for j in xrange(num_files)]
    # redistribute to make sure output files are exactly of the same size?

    # run through input file again, now split lines as assigned
    split_lines = [[] for i in xrange(args.split_num)]
    with open(fastq, 'r') as f:
        i = 0
        tup = []
        for line in f:
            if i%4 == 0:
                illumina_ID = line.split()[0]
                if args.id_sep:
                    ID = int(illumina_ID.split('.')[1])
                else:
                    ID = int(illumina_ID.strip('@'))
                tup = [illumina_ID + '\n']
            else:
                tup.append(line)

            if len(tup) == 4:
                rand_set = ID_to_set[ID]
                split_lines[rand_set].append(tup)
            i += 1

    # if given, run through second input file (paired reads)
    if args.input2:
        split_lines2 = [[] for i in xrange(args.split_num)]
        with open(args.input2, 'r') as f:
            i = 0
            tup = []
            for line in f:
                if i%4 == 0:
                    illumina_ID = line.split()[0]
                    if args.id_sep:
                        ID = int(illumina_ID.split('.')[1])
                    else:
                        ID = int(illumina_ID.strip('@'))
                    tup = [illumina_ID + '\n']
                else:
                    tup.append(line)

                if len(tup) == 4:
                    rand_set = ID_to_set[ID]
                    split_lines2[rand_set].append(tup)
                i += 1

    # write each set to its own file
    for i in xrange(args.split_num):
        if args.input2:
            filename = args.output + "1.%d.fastq" % i
        else:
            filename = args.output + ".%d.fastq" % i
        with open(filename, 'w') as f:
            for tup in split_lines[i]:
                f.write(''.join(tup))

    if args.input2:
        for i in xrange(args.split_num):
            filename = args.output + "2.%d.fastq" % i
            with open(filename, 'w') as f:
                for tup in split_lines2[i]:
                    f.write(''.join(tup))


if __name__ == '__main__':
        sys.exit(main())
