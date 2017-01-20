#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random
import math


__author__ = "Jasmijn Baaijens"

usage = """%prog [options] 

Analyze overlaps file (SAVAGE format).

"""

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument("--input", dest="input", type=str, required=True, help="input overlaps filename")
    args = parser.parse_args()
    analyze_overlaps(args.input)
    

def analyze_overlaps(filename):   
    pp_count = [0 for i in xrange(4)]
    ps_count = [0 for i in xrange(4)]
    sp_count = [0 for i in xrange(4)]
    ss_count = [0 for i in xrange(4)]
    c = 0
    with open(filename) as f:
        for line in f:
            c += 1
            line = line.strip().split('\t')
            if line[11] == 'p' and line[12] == 'p':
                if line[5] == '-' and line[6] == '+':
                    pp_count[0] += 1
                elif line[5] == '+' and line[6] == '-':
                    pp_count[1] += 1
                elif line[5] == '+' and line[6] == '+':
                    pp_count[2] += 1
                elif line[5] == '-' and line[6] == '-':
                    pp_count[3] += 1
                else:
                    print 'orientation not found...'
            elif line[11] == 'p' and line[12] == 's':
                if line[5] == '-' and line[6] == '+':
                    ps_count[0] += 1
                elif line[5] == '+' and line[6] == '-':
                    ps_count[1] += 1
                elif line[5] == '+' and line[6] == '+':
                    ps_count[2] += 1
                elif line[5] == '-' and line[6] == '-':
                    ps_count[3] += 1
                else:
                    print 'orientation not found...'
            elif line[11] == 's' and line[12] == 'p':
                if line[5] == '-' and line[6] == '+':
                    sp_count[0] += 1
                elif line[5] == '+' and line[6] == '-':
                    sp_count[1] += 1
                elif line[5] == '+' and line[6] == '+':
                    sp_count[2] += 1
                elif line[5] == '-' and line[6] == '-':
                    sp_count[3] += 1
                else:
                    print 'orientation not found...'
            elif line[11] == 's' and line[12] == 's':
                if line[5] == '-' and line[6] == '+':
                    ss_count[0] += 1
                elif line[5] == '+' and line[6] == '-':
                    ss_count[1] += 1
                elif line[5] == '+' and line[6] == '+':
                    ss_count[2] += 1
                elif line[5] == '-' and line[6] == '-':
                    ss_count[3] += 1
                else:
                    print 'orientation not found...'
            else:
                print 'read types not recognized...'
    
    print "Overlaps:"                
    print "[-+, +-, ++, --]"
    print "p-p: ", pp_count
    print "p-s: ", ps_count
    print "s-p: ", sp_count
    print "s-s: ", ss_count

    total = sum(pp_count) + sum(ps_count) + sum(sp_count) + sum(ss_count)
    print total
    print "# lines: ", c
    print "\n"
    return total
    
                        
if __name__ == '__main__':
        sys.exit(main())
