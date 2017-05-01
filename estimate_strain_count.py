#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random
import subprocess
from time import clock


__author__ = "Jasmijn Baaijens"

usage = """

Estimate the number of strains in a quasispecies based
on the contig alignments to a reference genome.
Assumes a sam file as input.

"""

# globals
verbose = False
allow_N = False


def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument('--sam', dest='infile_s', type=str)
    parser.add_argument('--ref', dest='reference', type=str)
#    parser.add_argument('--out', dest='outfile', type=str)
#    parser.add_argument('--readcount', dest='readcount', type=int, required=True)
#    parser.add_argument('--min_overlap_len', dest='min_overlap_len', type=int, default=0)
    parser.add_argument('--no_N', dest='allow_N', action='store_false')
    parser.add_argument('--verbose', dest='verbose', action='store_true')
    args = parser.parse_args()

    global verbose, allow_N
    verbose = args.verbose
    allow_N = args.allow_N
    min_overlap_len = 0

    if not (args.infile_s and args.reference):
        print "Specify input and output files."
        parser.print_help()

    output_dir = "strain_count"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    graph_file = output_dir + "/graph.txt"
    try:
        os.remove(graph_file)
    except OSError:
        pass

    cliques_file = "strain_count/cliques.txt"
    try:
        os.remove(cliques_file)
    except OSError:
        pass

    ref_list = []
    ref_dict = {}
    with open(args.reference, 'r') as f:
        lines = f.readlines()
        if len(lines) == 0:
            print "empty reference fasta... exiting."
            sys.exit(1)
        ref_id = ""
        ref_seq = ""
        idx = 0
        for line in lines:
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
        [sam_records, readcount] = read_sam_to_list(args.infile_s)
    else:
        sam_records = []
        readcount = 0

    # split sam records (single-end) per reference genome
    sam_records_s_per_ref = [[] for i in xrange(len(ref_list))]
    for record in sam_records:
        ref_id = record[2]
        ref_idx = ref_dict[ref_id]
        sam_records_s_per_ref[ref_idx].append(record)

    # find overlaps per reference genome and construct conflict graph
    for idx in xrange(len(ref_list)):
        ref_seq = ref_list[idx]
        sam_singles = sam_records_s_per_ref[idx]
        process_sam(ref_seq, sam_singles, graph_file, min_overlap_len, readcount)

    # run quick-cliques on conflict graph
    base_path = os.path.dirname(os.path.abspath(__file__))
    clique_command = base_path + "/quick-cliques/bin/qc --algorithm=degeneracy --input-file=strain_count/graph.txt > strain_count/cliques.txt 2>/dev/null"
    subprocess.check_call(clique_command, shell=True)

    # compute largest clique size
    max_clique_size = 0
    with open(cliques_file, 'r') as f:
        for line in f:
            clique = line.rstrip().split(' ')
            if all(node.isdigit() for node in clique):
                clique_size = len(clique)
                if clique_size > max_clique_size:
                    max_clique_size = clique_size
    print "Lower bound on the number of strains in this sample is %d." % max_clique_size
    return

#-------------------------------

def check_overlap(seq1, seq2, pos):
    mismatch_count = 0
    overlap_len = min(len(seq1)-pos, len(seq2))
    for i in xrange(overlap_len):
        base1 = seq1[pos + i]
        base2 = seq2[i]
        if allow_N and 'N' in [base1, base2]:
            continue
        if base1 != base2:
            mismatch_count += 1
    return mismatch_count

def power_find(n):
    result = []
    binary = bin(n)[:1:-1]
    for x in range(len(binary)):
        if int(binary[x]):
            result.append(2**x)
    return result

def read_sam_to_list(sam):
    records = []
    header = True
    ID_set = set()
    with open(sam, 'r') as f:
        unmapped = 0
        for line in f:
            if header and line[0] == '@':
                continue
            header = False
            aln = line.strip('\n').split('\t')
            [ID, FLAG, REF, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL] = aln[0:11]
            record = [ID, int(FLAG), REF, int(POS), int(MAPQ), CIGAR, RNEXT, int(PNEXT), int(TLEN), SEQ, QUAL]
            ID_set.add(ID)
            flag_decomp = power_find(int(FLAG))
            if 4 not in flag_decomp: # check if read is mapped
                if not 256 in flag_decomp: # ignore secondary aligments
                    records.append(record)
            else:
                unmapped += 1
        if verbose:
            print "Number of singles unmapped: ", unmapped
    return [records, len(ID_set)]


def get_overlap_line(read1, read2, pos, ovlen):
    assert pos >= 0
    assert ovlen >= 0
    # SAM format: ID FLAG REF POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
    FLAG1 = read1[1]
    bits1 = power_find(FLAG1)
    FLAG2 = read2[1]
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


def get_overlaps(record, active_reads, pos, min_overlap_len):
    [ID1, FLAG1, REF1, POS1, MAPQ1, CIGAR1, RNEXT1, PNEXT1, TLEN1, SEQ1, QUAL1] = record
    assert pos == POS1
    overlaps = []
    new_active_reads = []
    count_problems = 0
    for read in active_reads:
        [ID2, FLAG2, REF2, POS2, MAPQ2, CIGAR2, RNEXT2, PNEXT2, TLEN2, SEQ2, QUAL2] = read
        overlap_pos = POS1 - POS2
        assert overlap_pos >= 0
        overlap_len = min(len(SEQ2)-overlap_pos, len(SEQ1))
        if len(SEQ2) - overlap_pos >= min_overlap_len:
            new_active_reads.append(read)
        if overlap_len > min_overlap_len:
            overlap = get_overlap_line(read, record, overlap_pos, overlap_len)
            # check orientations
            ori1 = "-" if 16 in power_find(read[1]) else "+"
            ori2 = "-" if 16 in power_find(record[1]) else "+"
            overlap[5] = ori1
            overlap[6] = ori2
            # check for mismatches
            mismatch_count = check_overlap(SEQ2, SEQ1, int(overlap[2]))
            if mismatch_count > 0: # we only need mismatch-edges to build the conflict graph
                overlaps.append(overlap)
    return [overlaps, new_active_reads, count_problems]

def get_key_s(record):
    return record[3]

def process_sam(ref, sam_records, outfile, min_overlap_len, readcount):
    aln_count = len(sam_records)
    sorted_records = sorted(sam_records, key=get_key_s)
    if verbose:
        print "Total number of alignments: ", aln_count

    active_reads = []
    overlap_types = [0, 0, 0, 0] # [++, +-, -+, --]
    count_problems = 0
    edges = []
    i = 0
    cur_pos = 0
    overlap_count = 0
    while cur_pos < len(ref) and i < aln_count:
        cur_read = sorted_records[i]
        new_pos = cur_read[3]
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
            if line[0] < line[1]:
                edge = line[0:2]
            else:
                edge = [line[1], line[0]]
            if int(line[0]) != int(line[1]) and not edge in edges:
                edges.append(edge)
        i += 1

    with open(outfile, 'w') as f:
        f.write(str(readcount) + '\n')
        f.write(str(2*len(edges)) + '\n')
        for edge in edges:
            node1 = int(edge[0])
            node2 = int(edge[1])
            if node1 < 0 or node1 >= readcount:
                print node1
            if node2 < 0 or node2 >= readcount:
                print node2
            f.write(','.join(edge) + '\n')
            f.write(','.join(edge[::-1]) + '\n')

    if count_problems > 0 and verbose:
        print "# cases where overlap_pos2 < 0: ", count_problems
    if verbose:
        print "Total number of overlaps found: ", overlap_count
        print "... of which ++: ", overlap_types[0]
        print "... of which +-: ", overlap_types[1]
        print "... of which -+: ", overlap_types[2]
        print "... of which --: ", overlap_types[3]
        print "... duplicates: %d" % (overlap_count - len(edges))


if __name__ == '__main__':
    sys.exit(main())
