#!/usr/bin/env python3
from __future__ import division
from argparse import ArgumentParser
import sys, os
import vcf
import subprocess
import re
import shutil

usage = """Build contigs per block from a phased VCF using a reference sequence"""

def main():
    """Build contigs per block from a phased VCF using a reference sequence"""
    parser = ArgumentParser(description=usage)
    parser.add_argument('--vcf', dest='vcf', type=str, required=True, help='path to phased vcf (bgzipped and tabix-indexed)')
    parser.add_argument('-r', dest='ref', type=str, required=True, help='path to reference fasta')
    parser.add_argument('-c', dest='chrom', type=str, required=True, help='chromosome name in reference')
    parser.add_argument('--phaser', dest='phaser', action='store_true', help='vcf in phaser format')
    parser.add_argument('-o', dest='out_file', type=str, required=True, help='final output fasta')
    parser.add_argument('--region', dest='region', type=str, help='specify region of interest')
    parser.add_argument('--only-h1', dest='no_h2', action='store_true')
    parser.add_argument('--only-h2', dest='no_h1', action='store_true')
    args = parser.parse_args()

    # vcf_file = "k2_10x.phased.vcf.gz"
    # ref = "/export/scratch2/baaijens/hg38/hg38.chr6.fa"
    # chrom = "chr6"
    # phaser = False
    # out_file = vcf_file.rstrip("phased.vcf.gz") + ".contigs.fasta"

    if args.no_h1 and args.no_h2:
        print("Skipping h1 AND h2, nothing to do...")
        sys.exit(1)

    # read vcf
    vcf_reader = vcf.Reader(open(args.vcf, 'rb'))
    # check if tabix index is present, otherwise build index

    # count trailing Ns in reference
    if args.region:
        start_pos = int(args.region.split('-')[0])
        end_pos = int(args.region.split('-')[1])
    else:
        front_Ns, back_Ns, ref_len = count_ref_Ns(args.ref, args.chrom)
        start_pos = front_Ns + 1
        end_pos = ref_len - back_Ns

    # create tmp dir
    tmp_dir = "tmp"
    overwrite_dir(tmp_dir)

    if args.phaser:
        ps_old = -1
    else:
        ps_old = start_pos
    pos_old = start_pos
    block_id = 0
    for record in vcf_reader:
        fields = record.FORMAT.split(':')
        gt_data = record.samples[0]
        if args.phaser:
            GT_field_name = 'PG'
            PS_field_name = 'PI'
        else:
            GT_field_name = 'GT'
            PS_field_name = 'PS'

        if GT_field_name in fields:
            # check position
            pos = int(record.POS)
            if pos > end_pos:
                break
            elif pos < start_pos:
                continue
            # get phasing
            genotype = gt_data[GT_field_name].split('|')
            if len(genotype) <= 1: # unphased variant
                continue
            assert PS_field_name in fields
            phase_set = int(gt_data[PS_field_name])
            if phase_set != ps_old:
                # build contig
                if not args.no_h1:
                    subprocess.check_call(
                        "samtools faidx {} {}:{}-{} | bcftools consensus -H 1 {} > {}/{}.h1.fasta".format(
                            args.ref, args.chrom, pos_old, pos-1, args.vcf, tmp_dir, block_id
                        ),
                        shell=True)
                if not args.no_h2:
                    subprocess.check_call(
                        "samtools faidx {} {}:{}-{} | bcftools consensus -H 2 {} > {}/{}.h2.fasta".format(
                            args.ref, args.chrom, pos_old, pos-1, args.vcf, tmp_dir, block_id
                        ),
                        shell=True)
                ps_old = phase_set
                pos_old = pos
                block_id += 1

    # build contig for final block
    if not args.no_h1:
        subprocess.check_call(
            "samtools faidx {} {}:{}-{} | bcftools consensus -H 1 {} > {}/{}.h1.fasta".format(
                args.ref, args.chrom, pos_old, end_pos, args.vcf, tmp_dir, block_id
            ),
            shell=True)
    if not args.no_h2:
        subprocess.check_call(
            "samtools faidx {} {}:{}-{} | bcftools consensus -H 2 {} > {}/{}.h2.fasta".format(
                args.ref, args.chrom, pos_old, end_pos, args.vcf, tmp_dir, block_id
            ),
            shell=True)
    block_id += 1

    # concatenate resulting fasta files
    subprocess.check_call(
        "cat {0}/*.fasta > {0}/all.fasta".format(tmp_dir),
        shell=True
    )
    # rename reads
    subprocess.check_call(
        "~/Software/bbmap/rename.sh in={}/all.fasta out={} overwrite=t".format(
            tmp_dir, args.out_file
        ),
        shell=True
    )
    # remove intermediate files
    subprocess.check_call("rm -rf {}".format(tmp_dir), shell=True)



def overwrite_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
    else:
        shutil.rmtree(dir) #removes all the subdirectories!
        os.makedirs(dir)
    return


def count_ref_Ns(ref, chrom):
    with open(ref, 'r') as f:
        front_Ns = 0
        back_Ns = 0
        ref_len = 0
        nuc_seen = False
        process_chrom = False
        for line in f:
            if line[0] == '>':
                if line.strip().split()[0][1:] == chrom:
                    print(line.rstrip())
                    process_chrom = True
                else:
                    process_chrom = False
                continue
            elif not process_chrom:
                continue
            line = line.strip()
            ref_len += len(line)
            N_split = re.split("[^N]+", line)
            if not nuc_seen:
                # front N sequence
                front_Ns += len(N_split[0])
                if len(N_split) > 1:
                    nuc_seen = True
                    back_Ns = len(N_split[-1])
            elif len(N_split) > 1:
                back_Ns = len(N_split[-1])
            else:
                # back N sequence
                back_Ns += len(N_split[-1])
        print("front Ns: {}".format(front_Ns))
        print("back Ns: {}".format(back_Ns))
    return front_Ns, back_Ns, ref_len

if __name__ == '__main__':
    sys.exit(main())
