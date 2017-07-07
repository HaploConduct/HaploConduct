#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random
import subprocess
from time import clock
import itertools

__author__ = "Jasmijn Baaijens"

usage = """

Create an overlaps file for viral quasispecies assembly
based on the alignments to a reference genome. Assumes an
interleaved sam file as input.

"""

# globals
verbose = False

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument('--sam_s', dest='infile_s', type=str)
    parser.add_argument('--sam_p', dest='infile_p', type=str)
    parser.add_argument('--ref', dest='reference', type=str)
    parser.add_argument('--out', dest='outfile', type=str)
    parser.add_argument('--min_overlap_len', dest='min_overlap_len', type=int, default=0)
    parser.add_argument('--verbose', dest='verbose', action='store_true')
    args = parser.parse_args()

    global verbose
    verbose = args.verbose

    if not ((args.infile_s or args.infile_p) and args.outfile):
        print "Specify input and output files."
        parser.print_help()

    try:
        os.remove(args.outfile)
    except OSError:
        pass

    ref_list = []
    ref_dict = {}
    with open(args.reference, 'r') as f:
        lines = f.readlines()
        if len(lines) == 0:
            print "empty reference fasta... exiting."
            sys.exit(1)
#        elif len(lines) % 2 != 0:
#            print "invalid reference fasta... exiting."
        ref_id = ""
        ref_seq = ""
        idx = 0
        for line in lines:
#            if i%2 == 0:
            if line[0] == '>':
                if ref_seq != "" and ref_id != "":
                    ref_list.append(ref_seq)
                    ref_dict[ref_id] = idx
                    idx += 1
                    ref_seq = ""
                id_line = line.strip('\n')
                ref_id = id_line.split()[0][1:]
            else:
                ref_seq += line.strip('\n')
        if ref_seq == "":
            print "invalid fasta file... exiting"
            sys.exit(1)
        else:
            ref_list.append(ref_seq)
            ref_dict[ref_id] = idx
            idx += 1
            ref_seq = ""


    if args.infile_s:
        sam_records_s = read_sam_to_list(args.infile_s)
    else:
        sam_records_s = []
    if args.infile_p:
        sam_records_p = read_paired_sam_to_list(args.infile_p)
    else:
        sam_records_p = []

    # split sam records (single-end) per reference genome
    sam_records_s_per_ref = [[] for i in xrange(len(ref_list))]
    for record in sam_records_s:
        ref_id = record[2]
        ref_idx = ref_dict[ref_id]
        sam_records_s_per_ref[ref_idx].append(record)

    # split sam records (paired-end) per reference genome
    sam_records_p_per_ref = [[] for i in xrange(len(ref_list))]
    for record in sam_records_p:
        ref_id = record[0][2]
        ref_idx = ref_dict[ref_id]
        sam_records_p_per_ref[ref_idx].append(record)

    # find overlaps per reference genome
    for idx in xrange(len(ref_list)):
        ref_seq = ref_list[idx]
        sam_singles = sam_records_s_per_ref[idx]
        sam_paired = sam_records_p_per_ref[idx]
        process_sam(ref_seq, sam_singles, sam_paired, args.outfile, args.min_overlap_len)
    return


def power_find(n):
    try:
        int(n)
    except TypeError:
        print "power_find TypeError"
        print "n = ", n
        sys.exit(1)
    result = []
    binary = bin(n)[:1:-1]
    for x in range(len(binary)):
        if int(binary[x]):
            result.append(2**x)
    return result

def read_sam_to_list(sam):
    records = []
    header = True
    with open(sam, 'r') as f:
        unmapped = 0
        reverse = 0
        for line in f:
            if header and line[0] == '@':
                continue
            header = False
            aln = line.strip('\n').split('\t')
            [ID, FLAG, REF, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL] = aln[0:11]
            # check if read is mapped
            bits_flag = power_find(int(FLAG))
            if 4 in bits_flag:
                unmapped += 1
                continue
            elif 16 in bits_flag:
                reverse += 1
            # check for clipping
            cigar = ["".join(x) for _, x in itertools.groupby(CIGAR, key=str.isdigit)]
            softclip = CIGAR.split('S')[0]
            hardclip = CIGAR.split('H')[0]
            hardclip_end = CIGAR.split('H')[-1]
            if cigar[1] == 'S':
                # subtract soft-clipped part from alignment position
                cpos = int(POS) - int(cigar[0])
                cseq = SEQ
                cqual = QUAL
            elif cigar[1] == 'H':
                # subtract hard-clipped part from alignment position
                cpos = int(POS) - int(cigar[0])
                # add dummy nucleotides to front of sequence
                cseq = int(cigar[0])*'N' + SEQ
                cqual = int(cigar[0])*'$' + QUAL
            else:
                cpos = int(POS)
                cseq = SEQ
                cqual = QUAL
            if cigar[-1] == 'H':
                # add dummy nucleotides to back of sequence
                cseq +=	int(cigar[-2])*'N'
                cqual += int(cigar[-2])*'$'

            record = [ID, int(FLAG), REF, cpos, int(MAPQ), CIGAR, RNEXT, int(PNEXT), int(TLEN), cseq, cqual]
            records.append(record)
        if verbose:
            print "Number of singles unmapped: ", unmapped
            print "Number of singles reversed: ", reverse
    return records

def read_paired_sam_to_list(sam):
    records = []
    header = True
    with open(sam, 'r') as f:
        paired_read = []
        i = 0
        discarded = 0
        unmapped = 0
        reverse = 0
        for line in f:
            if header and line[0] == '@':
                continue
            header = False
            aln = line.strip('\n').split('\t')
            [ID, FLAG, REF, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL] = aln[0:11]
            # check if read is mapped
            if 4 in power_find(int(FLAG)):
                unmapped += 1
                continue
            # check for clipping
            cigar = ["".join(x) for _, x in itertools.groupby(CIGAR, key=str.isdigit)]
            softclip = CIGAR.split('S')[0]
            hardclip = CIGAR.split('H')[0]
            hardclip_end = CIGAR.split('H')[-1]
            if cigar[1] == 'S':
                # subtract soft-clipped part from alignment position
                cpos = int(POS) - int(cigar[0])
                cseq = SEQ
                cqual = QUAL
            elif cigar[1] == 'H':
                # subtract hard-clipped part from alignment position
                cpos = int(POS) - int(cigar[0])
                # add dummy nucleotides to front of sequence
                cseq = int(cigar[0])*'N' + SEQ
                cqual = int(cigar[0])*'$' + QUAL
            else:
                cpos = int(POS)
                cseq = SEQ
                cqual = QUAL
            if cigar[-1] == 'H':
                # add dummy nucleotides to back of sequence
                cseq += int(cigar[-2])*'N'
                cqual += int(cigar[-2])*'$'

            record = [ID, int(FLAG), REF, cpos, int(MAPQ), CIGAR, RNEXT, int(PNEXT), int(TLEN), cseq, cqual]
            paired_read.append(record)
            # store both sequences of a paired-end read together in a tuple
            assert len(paired_read) <= 2
            if i%2 == 1:
                if len(paired_read) == 2:
#                    print paired_read[0][0]
#                    print paired_read[1][0]
                    if paired_read[0][0] != paired_read[1][0]:
                        paired_read = [paired_read[1]]
                        discarded += 1
                        continue
                    elif (paired_read[0][3] >= paired_read[1][3]):
                        # reverse complement
                        if 16 in power_find(paired_read[0][1]) and 16 in power_find(paired_read[1][1]):
                            records.append([paired_read[1], paired_read[0], True])
                            reverse += 1
                        else:
                            discarded += 1
                    elif (paired_read[0][3] <= paired_read[1][3]):
                        # normal orientation
                        if (16 not in power_find(paired_read[0][1])) and (16 not in power_find(paired_read[1][1])):
                            records.append([paired_read[0], paired_read[1], False])
                        else:
                            discarded += 1
                else:
                    discarded += 1
                paired_read = []
            i += 1
        if verbose:
            print "Number of read ends discarded: ", discarded
            print "Number of read ends unmapped: ", unmapped
            print "Number of reverse complements considered: ", reverse
    return records

def compute_overlap_pos(pos1, pos2, len1, len2, CIGAR1, CIGAR2):
    # note: read #2 is in front!!!
    cigar1 = ["".join(x) for _, x in itertools.groupby(CIGAR1, key=str.isdigit)]
    cigar2 = ["".join(x) for _, x in itertools.groupby(CIGAR2, key=str.isdigit)]
    # compute length in alignment
    front_seq_len = 0
    front_ref_len = 0
    p = 0 # current position on reference
    total_back_ref_len = sum([int(cigar1[j]) if cigar1[j+1] != 'I' else 0 for j in range(0, len(cigar1), 2)])
    max_len = pos1 - pos2 + total_back_ref_len
    for i in range(0, len(cigar2), 2):
        aln_type = cigar2[i+1]
        aln_len = int(cigar2[i])
        if p < max_len:
            if aln_type != 'D':
                front_seq_len += min(aln_len, max_len - p)
            if aln_type != 'I':
                front_ref_len += min(aln_len, max_len - p)
                p += aln_len
    if front_ref_len <= (pos1 - pos2):
        # no actual overlap
        overlap_pos = -1
        overlap_len = 0
    else:
        back_ref_len = front_ref_len - (pos1 - pos2)
        back_seq_len = 0
        p = 0
        for i in range(0, len(cigar1), 2):
            aln_type = cigar1[i+1]
            aln_len = int(cigar1[i])
            assert aln_len > 0
            if p < back_ref_len:
                if aln_type != 'D':
                    back_seq_len += min(aln_len, back_ref_len - p)
                if aln_type != 'I':
                    p += aln_len
        overlap_pos = (pos1 - pos2) - ((front_ref_len - front_seq_len) - (back_ref_len - back_seq_len))
        if overlap_pos >= 0:
            overlap_len = min(len2 - overlap_pos, len1)
        else:
            # no acceptable overlap
            overlap_pos = -1
            overlap_len = 0
    return [overlap_pos, overlap_len]

def get_overlap_line(read1, read2, pos, ovlen):
#    assert pos >= 0
#    assert ovlen >= 0
    # SAM format: ID FLAG REF POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
    FLAG1 = int(read1[1])
    bits1 = power_find(FLAG1)
    FLAG2 = int(read2[1])
    bits2 = power_find(FLAG2)
    # overlap format: ID1 ID2 POS1 POS2 ORD ORI1 ORI2 PERC1 PERC2 LEN1 LEN2 TYPE1 TYPE2
    id1 = read1[0]
    id2 = read2[0]
    pos1 = str(pos)
    pos2 = "0"
    order = "-"
    ori1 = "-" if 16 in bits1 else "+"
    ori2 = "-" if 16 in bits2 else "+"
    seq1 = read1[9]
    seq2 = read2[9]
    perc = int(round(ovlen / min(len(seq1), len(seq2)) * 100))
    assert perc >= 0
    assert perc <= 100
    perc1 = str(perc)
    perc2 = "0"
    len1 = str(ovlen)
    len2 = "0"
    type1 = "s"
    type2 = "s"
#    line = '\t'.join([id1, id2, pos1, pos2, order, ori1, ori2, perc1, perc2, len1, len2, type1, type2])
#    return line
    overlap = [id1, id2, pos1, pos2, order, ori1, ori2, perc1, perc2, len1, len2, type1, type2]
    return overlap

def merge_overlaps(overlap1, overlap2, type1, type2):
    # if type1 == "s":
    #     if overlap1[5] != overlap2[5]:
    #         print "orientations1 don't match"
    #     if overlap1[6] != overlap2[6]:
    #         print "orientations2 don't match"
    overlap = overlap1
    overlap[11] = type1
    overlap[12] = type2
    if type1 == "p" and type2 == "p":
        # take care of ord
        if overlap1[0] != overlap2[0]:
            assert overlap1[0] == overlap2[1]
            overlap[4] = "2"
        else:
            overlap[4] = "1"
    overlap[3] = overlap2[2]
    overlap[8] = overlap2[7]
    overlap[10] = overlap2[9]
    return overlap


def get_overlaps(record, active_reads, pos, min_overlap_len):
    if len(record) == 3:
        record_paired = True
        [ID1, FLAG1, REF1, POS1, MAPQ1, CIGAR1, RNEXT1, PNEXT1, TLEN1, SEQ1, QUAL1] = record[0]
    else:
        record_paired = False
        [ID1, FLAG1, REF1, POS1, MAPQ1, CIGAR1, RNEXT1, PNEXT1, TLEN1, SEQ1, QUAL1] = record
    assert pos == POS1
    overlaps = []
    new_active_reads = []
    count_problems = 0
    for read in active_reads:
        if len(read) == 3:
            read_paired = True
            [ID2, FLAG2, REF2, POS2, MAPQ2, CIGAR2, RNEXT2, PNEXT2, TLEN2, SEQ2, QUAL2] = read[0]
        else:
            read_paired = False
            [ID2, FLAG2, REF2, POS2, MAPQ2, CIGAR2, RNEXT2, PNEXT2, TLEN2, SEQ2, QUAL2] = read
        overlap_pos = POS1 - POS2
        assert overlap_pos >= 0
        overlap_len = min(len(SEQ2)-overlap_pos, len(SEQ1))
        if len(SEQ2) - overlap_pos >= min_overlap_len:
            new_active_reads.append(read)
        [corrected_overlap_pos, corrected_overlap_len] = compute_overlap_pos(POS1, POS2, len(SEQ1), len(SEQ2), CIGAR1, CIGAR2)
        if corrected_overlap_len <= min_overlap_len or corrected_overlap_pos < 0:
            continue
        if record_paired == False and read_paired == False:
            overlap = get_overlap_line(read, record, corrected_overlap_pos, corrected_overlap_len)
            # check orientations
            ori1 = "-" if 16 in power_find(read[1]) else "+"
            ori2 = "-" if 16 in power_find(record[1]) else "+"
            overlap[5] = ori1
            overlap[6] = ori2
            overlaps.append(overlap)
        elif record_paired == True and read_paired == False:
            overlap1 = get_overlap_line(read, record[0], corrected_overlap_pos, corrected_overlap_len)
            POS3 = record[1][3]
            POS4 = read[3]
            SEQ3 = record[1][9]
            SEQ4 = read[9]
            CIGAR3 = record[1][5]
            CIGAR4 = read[5]
            [corrected_overlap_pos2, corrected_overlap_len2] = compute_overlap_pos(POS3, POS4, len(SEQ1), len(SEQ2), CIGAR1, CIGAR2)
            overlap2 = get_overlap_line(read, record[1], corrected_overlap_pos2, corrected_overlap_len2)
            overlap = merge_overlaps(overlap1, overlap2, "s", "p")
            # check orientations
            ori1 = "-" if 16 in power_find(read[1]) else "+"
            ori2 = "-" if record[2] else "+"
            overlap[5] = ori1
            overlap[6] = ori2
            if corrected_overlap_len2 > min_overlap_len and corrected_overlap_pos2 >= 0:
                overlaps.append(overlap)
        elif record_paired == False and read_paired == True:
            overlap1 = get_overlap_line(read[0], record, corrected_overlap_pos, corrected_overlap_len)
            overlap_pos2 = read[1][3] - record[3]
            if overlap_pos2 < 0:
                count_problems += 1
                continue
            POS3 = read[1][3]
            POS4 = record[3]
            SEQ3 = read[1][9]
            SEQ4 = record[9]
            CIGAR3 = read[1][5]
            CIGAR4 = record[5]
            [corrected_overlap_pos2, corrected_overlap_len2] = compute_overlap_pos(POS3, POS4, len(SEQ1), len(SEQ2), CIGAR1, CIGAR2)
            overlap2 = get_overlap_line(record, read[1], corrected_overlap_pos2, corrected_overlap_len2)
            overlap = merge_overlaps(overlap1, overlap2, "s", "p")
            # check orientations
            ori1 = "-" if read[2] else "+"
            ori2 = "-" if 16 in power_find(record[1]) else "+"
            overlap[5] = ori1
            overlap[6] = ori2
            if corrected_overlap_len2 > min_overlap_len and corrected_overlap_pos2 >= 0:
                overlaps.append(overlap)
        else:
            overlap1 = get_overlap_line(read[0], record[0], corrected_overlap_pos, corrected_overlap_len)
            overlap_pos2 = record[1][3] - read[1][3]
            if overlap_pos2 < 0:
                POS3 = read[1][3]
                POS4 = record[1][3]
                SEQ3 = read[1][9]
                SEQ4 = record[1][9]
                CIGAR3 = read[1][5]
                CIGAR4 = record[1][5]
                [corrected_overlap_pos2, corrected_overlap_len2] = compute_overlap_pos(POS3, POS4, len(SEQ3), len(SEQ4), CIGAR3, CIGAR4)
                overlap2 = get_overlap_line(record[1], read[1], corrected_overlap_pos2, corrected_overlap_len2)
                overlap = merge_overlaps(overlap1, overlap2, "p", "p")
                # check orientations
                ori1 = "-" if read[2] else "+"
                ori2 = "-" if record[2] else "+"
                overlap[5] = ori1
                overlap[6] = ori2
            else:
                POS3 = record[1][3]
                POS4 = read[1][3]
                SEQ3 = record[1][9]
                SEQ4 = read[1][9]
                CIGAR3 = record[1][5]
                CIGAR4 = read[1][5]
                [corrected_overlap_pos2, corrected_overlap_len2] = compute_overlap_pos(POS3, POS4, len(SEQ3), len(SEQ4), CIGAR3, CIGAR4)
                overlap2 = get_overlap_line(read[1], record[1], corrected_overlap_pos2, corrected_overlap_len2)
                overlap = merge_overlaps(overlap1, overlap2, "p", "p")
                # check orientations
                ori1 = "-" if read[2] else "+"
                ori2 = "-" if record[2] else "+"
                overlap[5] = ori1
                overlap[6] = ori2
            if corrected_overlap_len2 > min_overlap_len and corrected_overlap_pos2 >= 0:
                overlaps.append(overlap)
    return [overlaps, new_active_reads, count_problems]

def get_key_s(record):
    return record[3]

def get_key_p(record):
    return record[0][3]

def process_sam(ref, sam_records_s, sam_records_p, outfile, min_overlap_len):
    readcount_s = len(sam_records_s)
    readcount_p = len(sam_records_p)

    sorted_records_s = sorted(sam_records_s, key=get_key_s)
    sorted_records_p = sorted(sam_records_p, key=get_key_p)

    k1 = 0
    k2 = 0
    merged_records = []
    while k1 < readcount_s and k2 < readcount_p:
        single = sorted_records_s[k1]
        paired = sorted_records_p[k2]
        pos_s = single[3]
        pos_p = paired[0][3]
        if pos_s <= pos_p:
            merged_records.append([pos_s, single])
            k1 += 1
        else:
            merged_records.append([pos_p, paired])
            k2 += 1
    if k1 < readcount_s:
        assert k2 == readcount_p
        for single in sorted_records_s[k1:]:
            pos_s = single[3]
            merged_records.append([pos_s, single])
    elif k2 < readcount_p:
        assert k1 == readcount_s
        for paired in sorted_records_p[k2:]:
            pos_p = paired[0][3]
            merged_records.append([pos_p, paired])

    readcount = len(merged_records)
    if verbose:
        print "Total number of alignments: ", readcount
        print "... of which singles: ", readcount_s
        print "... of which paired: ", readcount_p
    assert readcount == readcount_s + readcount_p

    active_reads = []
    overlap_types = [0, 0, 0, 0] # [++, +-, -+, --]
    count_problems = 0
    with open(outfile, 'a') as outfile:
        i = 0
        cur_pos = merged_records[0][0]
        overlap_count = 0
        while cur_pos < len(ref) and i < readcount:
            cur_read = merged_records[i][1]
            new_pos = merged_records[i][0]
            assert new_pos >= cur_pos # records have to be sorted
            cur_pos = new_pos
            [overlaps, active_reads, subcount_problems] = get_overlaps(cur_read, active_reads, cur_pos, min_overlap_len)
            count_problems += subcount_problems
            active_reads.append(cur_read)
            for line in overlaps:
                overlap_count += 1
                if line[5] == "+" and line[6] == "+":
                    overlap_types[0] += 1
                elif line[5] == "+" and line[6] == "-":
                    overlap_types[1] += 1
                elif line[5] == "-" and line[6] == "+":
                    overlap_types[2] += 1
                elif line[5] == "-" and line[6] == "-":
                    overlap_types[3] += 1
                outfile.write('\t'.join(line) + '\n')
            i += 1
    if count_problems > 0 and verbose:
        print "# cases where overlap_pos2 < 0: ", count_problems
    if verbose:
        print "Total number of overlaps found: ", overlap_count
        print "... of which ++: ", overlap_types[0]
        print "... of which +-: ", overlap_types[1]
        print "... of which -+: ", overlap_types[2]
        print "... of which --: ", overlap_types[3]


if __name__ == '__main__':
    sys.exit(main())
