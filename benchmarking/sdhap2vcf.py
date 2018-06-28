#!/usr/bin/python3
import sys, os
import argparse

"""
Transform unphased VCF (input) to phased VCF (output) using phasing results
from SDhaP or H-PoP.
"""

def main():
    parser = argparse.ArgumentParser("build phased VCF from SDhaP output")
    parser.add_argument('--phased', dest='phased', type=str, required=True)
    parser.add_argument('--vcf', dest='vcf', type=str, required=True)
    parser.add_argument('-o', dest='outfile', type=str, required=True)
    parser.add_argument('--format', dest='format', type=str, default='sdhap')
    args = parser.parse_args()

    if args.format not in ['sdhap', 'hpop']:
        print("InputError: format must be sdhap or hpop")
        sys.exit(1)

    idx2phase = {}
    with open(args.phased, 'r') as phased:
        block = 0
        for line in phased:
            if line[0] == 'B':
                block += 1
                continue
            elif line[0] == '*':
                continue
            line = line.rstrip().split('\t')
            idx = int(line[0])
            if '-' in line[1:]:
                # not fully phased
                continue
            else:
                if args.format == 'sdhap':
                    haps = [str(int(x)-1) if x != '-' else '.' for x in line[1:]]
                else:
                    haps = [x if x != '-' else '.' for x in line[1:]]
            if max([int(x) for x in haps]) >= len(haps):
                # remove unusable output
                continue
            phase = "|".join(haps)
            idx2phase[idx] = [phase, block]

    ploidy = len(haps)
    print("ploidy = {}".format(ploidy))
    if ploidy > 4:
        print("ploidy > 4, need extra output VCF ---> extend script. Exiting.")
        sys.exit(1)
    print("{} largest index".format(max(idx2phase.keys())))

    # parse through VCF
    # substitute GT field by phased haplotypes
    # add phase set tag

    if args.format == 'sdhap':
        # 0-based vcf index
        idx = -1
    else:
        idx = 0
    unphased = 0
    artifacts = 0
    hom_ref = 0
    hom_alt = 0
    heterozygous = 0
    mnv_count = 0
    if ploidy > 2:
        vcf2 = open(args.outfile + ".2", 'w')
    with open(args.vcf, 'r') as vcf_in:
        with open(args.outfile, 'w') as vcf_out:
            for line in vcf_in:
                if line[0] == '#':
                    # copy header and add PS field description
                    if line[1] != '#':
                        vcf_out.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">\n')
                        if ploidy > 2:
                            vcf2.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">\n')
                    vcf_out.write(line)
                    if ploidy > 2:
                        vcf2.write(line)
                    continue
                elif len(line.split('\t')[3]) > 1 or len(line.split('\t')[4]) > 1:
                    mnv_count += 1
                    # continue
                # check genotype field
                line = line.rstrip()
                format = line.split('\t')[8]
                gt_idx = format.split(':').index('GT')
                data = line.split('\t')[9]
                gt = data.split(':')[gt_idx]
                pos = line.split('\t')[1]
                new_line = line.split('\t')[0:8]
                new_line2 = line.split('\t')[0:8]
                gt_set = set(gt.split('/'))
                if len(gt_set) == 1:
                    if sum([int(x) for x in set(gt.split('/'))]) == 0:
                        hom_ref += 1
                        continue
                    else:
                        hom_alt += 1
                else:
                    heterozygous += 1
                idx += 1
                data_split = data.split(':')
                data_split2 = data.split(':')
                if idx in idx2phase:
                    [phase, block] = idx2phase[idx]
                    if gt_set != set(phase.split('|')):
                        print(idx, pos, gt, phase, block)
                        artifacts += 1
                    if ploidy == 2:
                        data_split[gt_idx] = phase
                    elif ploidy == 3:
                        data_split[gt_idx] = phase[:3]
                        data_split2[gt_idx] = phase[2:]
                    else:
                        data_split[gt_idx] = phase[:3]
                        data_split2[gt_idx] = phase[4:]
                    data_split.append(str(block))
                    data_split2.append(str(block))
                    format += ":PS"
                else:
                    unphased += 1
                    if ploidy == 2:
                        data_split[gt_idx] = gt
                    elif ploidy == 3:
                        data_split[gt_idx] = gt[:3]
                        data_split2[gt_idx] = gt[2:]
                    else:
                        data_split[gt_idx] = gt[:3]
                        data_split2[gt_idx] = gt[4:]
                data = ":".join(data_split)
                new_line.append(format)
                new_line.append(data)
                data2 = ":".join(data_split2)
                new_line2.append(format)
                new_line2.append(data2)
                vcf_out.write('\t'.join(new_line) + '\n')
                if ploidy > 2:
                    vcf2.write('\t'.join(new_line2) + '\n')
    if ploidy > 2:
        vcf2.close()
    print("{} variants considered".format(idx+1))
    print("{} phased".format(len(idx2phase)))
    print("{} unphased".format(unphased))
    print("{} artifacts".format(artifacts))
    print("{}-0/0\t{}-0/1\t{}-1/1\t{}-mnv".format(
                    hom_ref, heterozygous, hom_alt, mnv_count))
    return

if __name__ == '__main__':
    sys.exit(main())
