#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random
import subprocess
from time import clock
from tqdm import tqdm # progress tracker


__author__ = "Jasmijn Baaijens"

usage = """%prog [options]

Randomly insert mutations in input fasta sequence(s) at a specified mutation rate.

"""

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument('--in', dest='infile', type=str)
    parser.add_argument('--out', dest='outfile', type=str)
    parser.add_argument('--sub', dest='sub_rate', type=float, default=0)
    parser.add_argument('--ins', dest='ins_rate', type=float, default=0)
    parser.add_argument('--del', dest='del_rate', type=float, default=0)
    parser.add_argument('--lb', dest='lb', type=int, default=0)
    parser.add_argument('--ub', dest='ub', type=int, default=-1)
    args = parser.parse_args()

    if args.ub >= 0:
        assert args.lb <= args.ub
        
    mut_rate = args.sub_rate + args.ins_rate + args.del_rate
    if not (args.infile and args.outfile and mut_rate > 0):
        print "Specify input and output files, as well as the desired mutation rates."
        parser.print_help()

    with open(args.infile, 'r') as f1:
        with open(args.outfile, 'w') as f2:
            i = 0
            seq = ""
            for line in f1:
                if line[0] == ">":
                    if len(seq) > 0:
                        print "Sequence length: {}".format(len(seq))
                        new_seq = seq[args.lb : args.ub]
                        new_seq = introduce_deletions(new_seq, args.del_rate)
                        new_seq = introduce_substitutions(new_seq, args.sub_rate)
                        new_seq = introduce_insertions(new_seq, args.ins_rate)
                        f2.write(seq[:args.lb] + new_seq + seq[args.ub:] + '\n')
                        i += 1
                    seq = ""
                    f2.write(line)
                else:
                    seq += line.strip('\n')
            if len(seq) > 0:
                print "Sequence length: {}".format(len(seq))
                new_seq = seq[args.lb : args.ub]
                new_seq = introduce_deletions(new_seq, args.del_rate)
                new_seq = introduce_substitutions(new_seq, args.sub_rate)
                new_seq = introduce_insertions(new_seq, args.ins_rate)
                f2.write(seq[:args.lb] + new_seq + seq[args.ub:] + '\n')
                i += 1
            print "Number of sequences processed: {}".format(i)


def introduce_substitutions(genome, mut_rate):
    nucs = ["A", "C", "T", "G", "A", "C", "T", "G"]
    mut_genome = ""
    mut_positions = []
    random.seed()
    genome_len = len(genome)
    mut_count = int(round(genome_len * mut_rate))
    # print "Sequence length: %d" % genome_len
    # print "Mutation count: %d" % mut_count
    while len(mut_positions) < mut_count:
        k = random.randint(0, genome_len-1)
        if k not in mut_positions:
            mut_positions.append(k)
    old_pos = -1
    for pos in tqdm(sorted(mut_positions)):
        mut_genome += genome[old_pos+1 : pos]
        base = genome[pos]
        if base in nucs:
            base_idx = nucs.index(genome[pos])
            k2 = random.randint(1, 3)
        else:
            base_idx = 0
            k2 = random.randint(1, 4)
        mut_genome += nucs[base_idx + k2]
        old_pos = pos
    if old_pos+1 < len(genome):
        mut_genome += genome[old_pos+1 :]
    assert len(genome) == len(mut_genome)
    print "{} substitutions done".format(mut_count)
    return mut_genome


def introduce_insertions(genome, mut_rate):
    random.seed()
    genome_len = len(genome)
    mut_count = int(round(genome_len * mut_rate))
    mut_genome = ""
    mut_positions = {}
    while len(mut_positions) < mut_count:
        k = random.randint(0, genome_len-1)
        size = random.randint(1, 10)
        if k not in mut_positions:
            mut_positions[k] = size
    old_k = 0
    for k in sorted(mut_positions.iterkeys()): # sort by position
        size = mut_positions[k]
        mut_genome += genome[old_k:k] + random_seq(size)
        old_k = k
    print "{} insertions done".format(mut_count)
    return mut_genome


def random_seq(size):
    nucs = ["A", "C", "T", "G"]
    seq = ""
    random.seed()
    while len(seq) < size:
        seq += nucs[random.randint(0, 3)]
    return seq


def introduce_deletions(genome, mut_rate):
    random.seed()
    genome_len = len(genome)
    mut_count = int(round(genome_len * mut_rate))
    mut_genome = genome
    i = 0
    while i < mut_count:
        i += 1
        k = random.randint(0, genome_len-1)
        size = random.randint(1, min(10, genome_len-k))
        mut_genome = genome[:k] + genome[k+size:]
        genome_len -= size
        # print size, len(genome), len(mut_genome)
        assert genome_len == len(mut_genome)
        genome = mut_genome
    print "{} deletions done".format(mut_count)
    return mut_genome


if __name__ == '__main__':
    sys.exit(main())
