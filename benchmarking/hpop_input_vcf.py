#!/usr/bin/python3
import sys, os
import argparse


"""
Convert a standard VCF file to a VCF accepted by H-PoP for polyploid haplotype
assembly: removes all homozygous calls and any format fields other than 'GT'
"""

def main():
    parser = argparse.ArgumentParser("convert standard VCF to H-PoP VCF")
    parser.add_argument('-i', dest='infile', type=str, required=True)
    parser.add_argument('-o', dest='outfile', type=str, required=True)
    args = parser.parse_args()

    with open(args.infile, 'r') as vcf_in:
        with open(args.outfile, 'w') as vcf_out:
            for line in vcf_in:
                if line[0] == '#':
                    # copy header
                    vcf_out.write(line)
                    continue
                # check genotype field
                format = line.split('\t')[8]
                gt_idx = format.split(':').index('GT')
                data = line.split('\t')[9]
                gt = data.split(':')[gt_idx]
                if len(set(gt.split('/'))) == 1:
                    # homozygous call
                    continue
                # remove additional fields
                new_line = line.split('\t')[0:8]
                new_line.append('GT')
                new_line.append(gt)
                vcf_out.write('\t'.join(new_line) + '\n')
    return

if __name__ == '__main__':
    sys.exit(main())
