#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import subprocess

__author__ = "Jasmijn Baaijens"

usage = """

Combine contigs from SAVAGE Stage a for continuing the assembly algorithm.
This script produces a combined contig file as well as the corresponding
subreads file.

"""

def main():
    print "combine_contigs.py"
    parser = ArgumentParser(description=usage)
    parser.add_argument('--split', dest='split_num', type=int, required=True,
                            help='specify the number of patches')
    parser.add_argument('--paired_to_single', dest='paired_to_single', action='store_true',
                            help='reuse paired-end contigs as single-end contigs')
    parser.add_argument('--paired_to_paired', dest='paired_to_paired', action='store_true',
                            help='reuse paired-end contigs as single-end contigs')
    # parser.add_argument('--out_contigs', dest='out_contigs', type=str, required=True,
    #                         help='output file for combined contigs')
    # parser.add_argument('--out_subreads', dest='out_subreads', type=str, required=True,
    #                         help='output file for combined subreads')
    args = parser.parse_args()
    if args.split_num < 1:
        print "--split must be an integer >= 1"
        sys.exit(1)
    if args.paired_to_single and args.paired_to_paired:
        print "options --paired_to_single and --paired_to_paired are mutually exclusive"
        sys.exit(1)
    # combine contigs and update subread info
    [SE_contigs, PE_contigs] = combine_contigs(args.split_num, args.paired_to_single, args.paired_to_paired)
#    print SE_contigs
#    print PE_contigs
#    print SE_contigs + 2*PE_contigs
    return

def combine_contigs(split_num, paired_to_single, paired_to_paired):
    current_SE_contigs = 0
    current_PE_contigs = 0
    # count the total number of single-end contigs
    subprocess.check_call("cat stage_a/patch*/stage_a/singles.fastq > stage_a/tmp_singles.fastq", shell=True)
    final_SE_contigs = round(file_len('stage_a/tmp_singles.fastq')/4)
    subprocess.check_call("rm stage_a/tmp_singles.fastq", shell=True)
    if os.path.exists("stage_a/combined_singles.fastq"):
        subprocess.check_call("rm stage_a/combined_singles.fastq", shell=True)
    # now build new contig & subreads files
    ALL_TYPES = ["singles", "paired1", "paired2"]
    with open('stage_a/subreads.txt', 'w') as new_subreads_file:
        s_id = 0
        for patch_num in range(split_num):
            # build combined fastq files
            if paired_to_single:
                # concatenate all contigs to combined singles file
                subprocess.check_call("cat stage_a/patch%d/stage_a/singles.fastq >> stage_a/combined_singles.fastq" % patch_num, shell=True)
                subprocess.check_call("cat stage_a/patch%d/stage_a/paired1.fastq >> stage_a/combined_singles.fastq" % patch_num, shell=True)
                subprocess.check_call("cat stage_a/patch%d/stage_a/paired2.fastq >> stage_a/combined_singles.fastq" % patch_num, shell=True)
            elif paired_to_paired:
                # concatenate contigs per type to combined fastq files
                subprocess.check_call("cat stage_a/patch%d/stage_a/singles.fastq >> stage_a/combined_singles.fastq" % patch_num, shell=True)
                subprocess.check_call("cat stage_a/patch%d/stage_a/paired1.fastq >> stage_a/combined_paired1.fastq" % patch_num, shell=True)
                subprocess.check_call("cat stage_a/patch%d/stage_a/paired2.fastq >> stage_a/combined_paired2.fastq" % patch_num, shell=True)
            else:
                subprocess.check_call("cat stage_a/patch%d/stage_a/singles.fastq >> stage_a/combined_singles.fastq" % patch_num, shell=True)
            # count number of contigs in this patch
            singles_count = int(file_len('stage_a/patch%d/stage_a/singles.fastq' % patch_num)/4)
            paired_count = int(file_len('stage_a/patch%d/stage_a/paired1.fastq' % patch_num)/4)
            renamed2originals = {}
            i = 0
            # first convert subread IDs to original IDs (before splitting the data)
            # this avoids duplicate IDs when combining the separate patches
            original_singles = 'stage_a/original_reads/%s.%d.fastq' % (ALL_TYPES[0], patch_num)
            if os.path.isfile(original_singles):
                with open(original_singles, 'r') as f1:
                    for line in f1:
                        if i%4 == 0:
#                            old_id = line.strip('\n')[1:]
                            old_id = str(s_id) # note: old_id equals line number in original fastq to avoid problems with non-int identifiers
                            s_id += 1
                            new_id = str(int(round(i/4)))
                            renamed2originals[new_id] = old_id
                        i += 1
                    assert i%4 == 0
            original_pairs = 'stage_a/original_reads/%s.%d.fastq' % (ALL_TYPES[1], patch_num)
            if os.path.isfile(original_pairs):
                with open(original_pairs, 'r') as f2:
                    for line in f2:
                        if i%4 == 0:
#                            old_id = line.strip('\n')[1:]
                            old_id = str(s_id) # note: old_id equals line number in original fastq to avoid problems with non-int identifiers
                            s_id += 1
                            new_id = str(int(round(i/4)))
                            renamed2originals[new_id] = old_id
                        i += 1
            # keep track of contig counts for paired_to_single contig IDs
            split_start = current_SE_contigs + 2*current_PE_contigs
#            split_p1_start = split_start + singles_count
            split_p2_start = split_start + paired_count
            # convert subreads info
            with open('stage_a/patch%d/stage_a/subreads.txt' % patch_num, 'r') as f3:
                for line in f3:
                    split_line = line.strip('\n').split('\t')
                    contig_id = split_line[0]
                    if paired_to_single:
                        if int(contig_id) < singles_count: # single-end contig
                            new_contig_id = str(int(contig_id) + split_start)
                            new_line = [new_contig_id]
                            for subread_info in split_line[1:]:
                                [ID, ori, poslist, lenlist] = subread_info.split(':')
                                new_subread_id = renamed2originals[ID]
                                new_info = new_subread_id + ':' + ori + ':' + poslist + ':' + lenlist
                                new_line.append(new_info)
                            new_subreads_file.write('\t'.join(new_line) + '\n')
                        else:
                            # split paired-end contig into 2 single-end contigs
                            new_contig_id_1 = str(int(contig_id) + split_start)
                            new_contig_id_2 = str(int(contig_id) + split_p2_start)
                            new_line_1 = [new_contig_id_1]
                            new_line_2 = [new_contig_id_2]
                            for subread_info in split_line[1:]:
                                [ID, ori, poslist, lenlist] = subread_info.split(':')
                                new_subread_id = renamed2originals[ID]
                                [index1, index2] = poslist.split(',')
                                [len1, len2] = lenlist.split(',')
                                new_info_1 = new_subread_id + ':' + ori + ':' + index1 + ':' + len1
                                new_info_2 = new_subread_id + ':' + ori + ':' + index2 + ':' + len2
                                new_line_1.append(new_info_1)
                                new_line_2.append(new_info_2)
                            new_subreads_file.write('\t'.join(new_line_1) + '\n')
                            new_subreads_file.write('\t'.join(new_line_2) + '\n')
                    else: # keep subread info as it is, just update contig IDs
                        if int(contig_id) < singles_count: # single-end contig
                            new_contig_id = str(int(contig_id) + current_SE_contigs)
                        elif paired_to_paired: # paired-end contig
                            new_contig_id = str(int(contig_id) + final_SE_contigs + current_PE_contigs)
                        else:
                            continue
                        new_line = [new_contig_id]
                        for subread_info in split_line[1:]:
                            [ID, ori, poslist, lenlist] = subread_info.split(':')
                            new_subread_id = renamed2originals[ID]
                            new_info = new_subread_id + ':' + ori + ':' + poslist + ':' + lenlist
                            new_line.append(new_info)
                        new_subreads_file.write('\t'.join(new_line) + '\n')
            current_SE_contigs += singles_count
            current_PE_contigs += paired_count
    return [current_SE_contigs, current_PE_contigs]


def file_len(fname):
    with open(fname) as f:
        i = 0
        for i, l in enumerate(f):
            pass
    if i > 0:
        linecount = i + 1
    else:
        linecount = 0
    assert linecount % 4 == 0
    return linecount

if __name__ == '__main__':
    sys.exit(main())
