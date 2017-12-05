#!/usr/bin/env python3

from __future__ import division
import sys
from argparse import ArgumentParser
from scipy.stats import norm
from math import floor
import operator as op
from functools import reduce

# input:
# - read length
# - internal segment size (mean, stddev)
# - haplotype coverage

# output:
# a table with two columns (distance, expected evidence, threshold) and a line
# for every possible distance (1,...,L), where L is maximum distance with an
# expected amount of evidence greater than zero.


__author__ = "Jasmijn Baaijens"

usage = """

Build a table of minimum evidence values per variation distance, depending
on read length, internal segment size and haplotype coverage.

"""

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument('-l', dest='readlen', type=float, required=True,
                        help="sequence length (per read end)")
    parser.add_argument('-i', dest='intseg', type=float, required=True,
                        help="mean internal segment size")
    parser.add_argument('-s', '--stddev', dest='stddev', type=float, required=True,
                        help="standard deviation from mean internal segment size")
    parser.add_argument('-c', '--hcov', dest='hcov', type=float, required=True,
                        help="average coverage per haplotype")
    parser.add_argument('-o', dest='outfile', type=str,
                        default="exp_ev_table.tsv", help="store table here")
    parser.add_argument('-v', dest='verbose', action='store_true',
                        help="verbose mode")
    args = parser.parse_args()

    readlen = args.readlen
    intseg = args.intseg
    stddev = args.stddev
    hcov = args.hcov

    fragsize = intseg + 2*readlen
    assert fragsize > 0

    seq_err = 0.01
    accuracy = 1e-3

    outfile = open(args.outfile, 'w')
    outfile.write("# INPUT:\n")
    outfile.write("# readlen {}\n".format(readlen))
    outfile.write("# intseg {}\n".format(intseg))
    outfile.write("# stddev {}\n".format(stddev))
    outfile.write("# hcov {}\n".format(hcov))
    outfile.write("# OUTPUT:\n")
    outfile.write("# dist\texp_ev\tmin_ev\n")

    # set initial values
    exp_ev = 0
    old_exp_ev = 0
    dist = 1
    exp_ev_list = []
    norm_cdf = {} # store computed cdf values here for reuse (major speedup)

    print("Computing expected evidence values")
    #print("Evidence boundaries:", end=" ")
    while exp_ev > 0 or dist == 1:
        exp_ev = 0
        # count evidence from single-end sequences
        exp_ev += hcov * max(0, fragsize - dist) / fragsize
        # count evidence from paired-end sequences;
        p_l = []
        for x in range(0, int(floor(readlen))):
            p, norm_cdf = prob(x, dist, readlen, intseg, stddev, norm_cdf)
            p_l.append(p)
        exp_ev += hcov * sum(p_l) / readlen
        exp_ev = int(floor(exp_ev))
        exp_ev_list.append(exp_ev)
        if exp_ev == 0:
            break
        elif exp_ev != old_exp_ev:
            #print(dist, end=" ")
            old_exp_ev = exp_ev
        dist += 1
        if dist > fragsize + 2*stddev:
            break

    print("Compute corresponding thresholds in increasing order")
    # compute min_evidence values in increasing order
    ev_values = sorted(set(exp_ev_list))
    ev_to_threshold = {}
    min_ev = 1
    for exp_ev in ev_values:
        min_ev = findMinEv(exp_ev, min_ev, seq_err, accuracy)
        ev_to_threshold[exp_ev] = min_ev

    print("Write results to {}".format(args.outfile))
    # now write output
    for i, exp_ev in enumerate(exp_ev_list):
        dist = i+1
        min_ev = ev_to_threshold[exp_ev]
        outfile.write("{}\t{}\t{}\n".format(dist, exp_ev, min_ev))
    outfile.close()
    return



def prob(x, dist, readlen, intseg, stddev, norm_cdf):
    min_insert = dist - 2*readlen + x + 1
    if min_insert in norm_cdf:
        p1 = norm_cdf[min_insert]
    else:
        p1 = norm(intseg, stddev).cdf(min_insert)
        norm_cdf[min_insert] = p1
    max_insert = dist - readlen + x
    if max_insert in norm_cdf:
        p2 = norm_cdf[max_insert]
    else:
        p2 = norm(intseg, stddev).cdf(max_insert)
        norm_cdf[max_insert] = p2
    return p2 - p1, norm_cdf

def findMinEv(c, m1, seq_err, accuracy):
    p1 = 0
    for m in range(m1, c):
        p1 += choose(c,m)*seq_err**(m)*(1-seq_err)**(c-m)

    while p1 > accuracy:
        m1 += 1
        p = 0
        for m in range(m1, c):
            p += choose(c,m)*seq_err**(m)*(1-seq_err)**(c-m)
        p1 = p
    return m1

def choose(n,k):
    k = min(k, n-k)
    if k == 0:
        return 1
    numer = reduce(op.mul, range(n, n-k, -1))
    denom = reduce(op.mul, range(1, k+1))
    return numer//denom

    return min_ev

if __name__ == '__main__':
    sys.exit(main())
