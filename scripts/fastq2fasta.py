#!/usr/bin/env python
from __future__ import division
from optparse import OptionParser
import os
import sys
import random


__author__ = "Jasmijn Baaijens"

usage = """%prog <infile> <outfile>

 Convert a fastq file to a fasta file.

    <infiles>           input fasta file
    <outfiles>          output fastq file

"""


def main():
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()
    
    if len(args) < 2:
        parser.print_help()
        return 1
          
    infile = os.path.abspath(args[0])
    outfile = os.path.abspath(args[1])
        
    with open(infile, 'rb') as f1:
        with open(outfile, 'wb') as f2:
            c = 0
            for line in f1:
                c += 1
                if (c%4) == 1:
                    f2.write('>' + line[1:])
                elif (c%4) == 2:
                    f2.write(line)
                     
if __name__ == '__main__':
        sys.exit(main())
