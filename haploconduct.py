#!/usr/bin/env python
from __future__ import division
import argparse
import subprocess
import os
import sys


__author__ = "Jasmijn Baaijens"
__license__ = "GPL"

usage = """
Program: HaploConduct
Version: 0.2
Release date: December 23, 2019
GitHub: https://github.com/HaploConduct/HaploConduct

HaploConduct is a package designed for reconstruction of individual haplotypes
from next generation sequencing data, in particular Illumina. HaploConduct
consists of two methods: SAVAGE and POLYTE.

SAVAGE is a computational tool for reconstructing individual haplotypes of
intra-host virus strains, a viral quasispecies (of unknown ploidy).

POLYTE aims to reconstruct haplotigs for diploid and polyploid genomes of known
ploidy and comes in two modes: polyte and polyte-split. The latter mode is
intended for regions larger than 100 kb. This pipeline splits the data into
regions (bins) of 10 kb, runs POLYTE on each region, and combines the result
into a final assembly.

Run haploconduct -h for a complete description.

For the complete manual, please visit https://github.com/HaploConduct/HaploConduct
"""



def main():
    parser =  argparse.ArgumentParser(prog=usage)
    parser.add_argument('subprogram',
                        choices=['savage', 'polyte', 'polyte-split'],
                        help='HaploConduct offers three modes: savage, polyte, and polyte-split. Please see the manual for details.')

    # test if the haploconduct core functions properly
    base_path = os.path.dirname(os.path.abspath(__file__))
    binary_help = base_path + '/bin/ViralQuasispecies --help'
    try:
        subprocess.check_output(binary_help, stderr=subprocess.STDOUT,
                                shell=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.output.decode() + '\n')
        sys.stderr.flush()
        sys.exit(1)

    if len(sys.argv[1:])==0:
        print usage
        parser.exit()
    elif '-h' in sys.argv or '--help' in sys.argv:
        if sys.argv[1] == 'savage':
            subprocess.check_call('{}/savage.py -h'.format(base_path),
                                  shell=True)
            parser.exit()
        elif sys.argv[1] == 'polyte':
            subprocess.check_call('{}/polyte.py -h'.format(base_path),
                                  shell=True)
            parser.exit()
        elif sys.argv[1] == 'polyte-split':
            subprocess.check_call('{}/polyte-split.py -h'.format(base_path),
                                  shell=True)
            parser.exit()

    args, params = parser.parse_known_args()

    if args.subprogram == 'savage':
        subprocess.check_call(['{}/savage.py'.format(base_path)] + params)
    elif args.subprogram == 'polyte':
        subprocess.check_call(['{}/polyte.py'.format(base_path)] + params)
    else:
        subprocess.check_call(['{}/polyte-split.py'.format(base_path)] + params)

    if params:
        print "HaploConduct completed successfully. Have a nice day:)\n"
    return


if __name__ == '__main__':
    sys.exit(main())
