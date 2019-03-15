#!/usr/bin/python3
import sys, os
import subprocess
import itertools # groupby()
import argparse

# input preparation:
# 1. split fasta over separate files
# 2. perform pairwise alignments of haplotypes: nucmer hap1.fasta hap2.fasta
# 3. extract polymorphic positions: show-snps -C -T -H out.delta > out.show-snps.tsv

def main():
    parser = argparse.ArgumentParser(prog='polymorphic_positions.py',
                description='Evaluate assemblies at polymorphic positions')
    parser.add_argument('-s', '--snps', dest='snps_file', type=str, required=True)
    parser.add_argument('-t', '--truth', dest='truth_file', type=str, required=True)
    parser.add_argument('-c', '--contigs', dest='contig_file', type=str, required=True)
    parser.add_argument('--skip_indels', dest='skip_indels', action='store_true')
    args = parser.parse_args()
    # input data
    skip_indels = args.skip_indels
    # read SNPs and store by index in both haplotypes
    hap2snplist = read_snps(args.snps_file, skip_indels)
    # map contigs to haplotypes (bwa mem)
    sam_file = "{}.truth.sam".format(args.contig_file)
    # print(sam_file)
    if not os.path.exists(sam_file):
        print("\nAligning contigs against ground truth...")
        if not os.path.exists(args.truth_file + ".bwt"):
            subprocess.check_call(
                "bwa index {} 2>/dev/null".format(args.truth_file),
                shell=True)
        subprocess.check_call(
            "bwa mem -L1000 -t 12 {0} {1} > {2} 2>/dev/null".format(
                args.truth_file, args.contig_file, sam_file),
            shell=True)
    # read contigs and truth
    ground_truth = read_fasta(args.truth_file)[0]
    contigs = read_fasta(args.contig_file)[0]
    # trace contig alignments and store errors
    [alignments, unmapped_ids] = read_sam(sam_file)
    print("#contigs: {}".format(len(contigs)))
    print("#alignments: {}".format(len(alignments)))
    print("#unmapped: {}".format(len(unmapped_ids)))
    contigs2aln = {c : [] for c in contigs}
    print("Processing alignments")
    for aln in alignments:
        [contig_id, ref_id] = aln[0:2]
        contig_seq = contigs[contig_id]
        truth_seq = ground_truth[ref_id]
        error_info = check_alignment(aln, contig_seq, truth_seq, skip_indels)
        contigs2aln[contig_id].append(error_info)
    os.remove(sam_file)

    # calculate statistics:
    # - # true polymorhisms in assembly / # true polymorphisms
    # - # true polymorhisms in contigs / # true polymorphisms covered
    # - # assembly errors at polymorphic positions / # assembly errors
    truth2aln = {hap : [] for hap in hap2snplist}
    num_assembly_errors = 0
    num_snp_assembly_errors = 0
    frac_variants_assembled = [0, 0]
    snps_covered = {hap : {snp : False for snp in snplist}
                    for hap, snplist in hap2snplist.items()}
    print("Calculating statistics")
    progress = 0
    ncontigs = len(contigs)
    for contig, aln in contigs2aln.items():
        progress += 1
        if len(aln) == 0:
            continue
        # select alignment with smallest number of error positions in case of
        # multiple alignments
        best_aln = min(aln, key = lambda x: len(x[0]))
        [error_indices, seq_id, ref_id, truth_start, truth_end] = best_aln
        truth2aln[ref_id].append(best_aln)
        num_assembly_errors += len(error_indices)
        for idx in error_indices:
            if idx in hap2snplist[ref_id]:
                num_snp_assembly_errors += 1
        num_variants_contig = 0
        num_variants_truth = 0
        # print(error_indices)
        # print(hap2snplist)
        for snp in hap2snplist[ref_id]:
            if snp < truth_start or snp > truth_end:
                continue
            num_variants_truth += 1
            if snp not in error_indices:
                num_variants_contig += 1
                snps_covered[ref_id][snp] = True
        frac_variants_assembled[0] += num_variants_contig
        frac_variants_assembled[1] += num_variants_truth
        if ncontigs > 10 and progress % int(round(ncontigs/10)) == 0:
            print("#contigs processed:", progress)
    num_snps_covered = sum([sum(snp_dict.values())
                            for hap, snp_dict in snps_covered.items()])
    num_snps = sum([len(snp_dict) for hap, snp_dict in snps_covered.items()])
    print("\nRaw numbers:")
    print("num_snps_covered {}".format(num_snps_covered))
    print("num_snps {}".format(num_snps))
    print("frac_variants_assembled {}".format(frac_variants_assembled))
    print("num_snp_assembly_errors {}".format(num_snp_assembly_errors))
    print("num_assembly_errors {}".format(num_assembly_errors))
    print("\nOutput stats:")
    print("# true polymorphisms in assembly / # true polymorphisms = {}".format(
            num_snps_covered / num_snps
    ))
    print("# true polymorphisms in contigs /  # true polymorphisms covered"
          " = {}".format(frac_variants_assembled[0] / frac_variants_assembled[1]
    ))
    print("# assembly errors at polymorphic positions / # assembly errors"
          " = {}".format(num_snp_assembly_errors / num_assembly_errors)
    )
    print()

    return


def check_alignment(aln, contig, truth, skip_indels):
    # trace alignment and report all positions in truth where errors occur
    [seq_id, ref_id, pos, ori, cigar] = aln
    error_indices = []
    if ori == "-":
        contig = revcomp(contig)
    # split cigar into numbers and characters all separately
    splitcigar = ["".join(x) for _, x in itertools.groupby(cigar,
                    key=str.isdigit)]
    # trace alignment and store SNP positions
    contig_pos = 0
    truth_pos = pos
    start_pos = pos
    truth_start = pos
    i = 0
    while i+1 < len(splitcigar):
        aln_len = int(splitcigar[i])
        aln_type = splitcigar[i+1]
        local_error_indices = []
        if aln_type == "S" or aln_type == "H":
            # print("WARNING: Clipped alignment")
            if i == 0: # front end clipped
                clipped_truth = min(aln_len, pos)
                start_pos -= clipped_truth
                contig_pos += aln_len
            elif i > 0: # back end clipped
                clipped_truth = min(aln_len, len(truth)-truth_pos)
                contig_pos += aln_len
                truth_pos += clipped_truth
            # compare (partially) aligned sequences
            truth_start = truth_pos-clipped_truth
            contig_part = contig[contig_pos-clipped_truth : contig_pos]
            truth_part = truth[truth_pos-clipped_truth : truth_pos]
            local_error_indices = find_SNPs(
                contig_part, truth_part, truth_pos-clipped_truth)
        elif aln_type == "I":
            local_error_indices.append(truth_pos)
            contig_pos += aln_len
        elif aln_type == "D":
            for j in range(aln_len):
                local_error_indices.append(truth_pos+1+j)
            truth_pos += aln_len
        elif aln_type == "M":
            if contig_pos + aln_len > len(contig):
                print("contig {} too short?".format(seq_id))
                print(contig_pos, aln_len, len(contig))
                print(pos, cigar)
            if truth_pos + aln_len > len(truth):
                print("truth too short?")
                print(truth_pos, aln_len, len(truth))
                print(pos, cigar)
            contig_part = contig[contig_pos : contig_pos + aln_len]
            truth_part = truth[truth_pos : truth_pos + aln_len]
            local_error_indices = find_SNPs(contig_part, truth_part, truth_pos)
            truth_pos += aln_len
            contig_pos += aln_len
        else:
            print("ERROR: cigar string not recognized.")
            sys.exit(1)
        error_indices += local_error_indices
        i += 2
    truth_end = truth_pos
    return [error_indices, seq_id, ref_id, truth_start, truth_end]


def find_SNPs(contig, truth, truth_pos):
    error_indices = []
    assert(len(contig) == len(truth))
    for i in range(len(contig)):
        b1 = contig[i].upper()
        b2 = truth[i].upper()
        if b1 != b2 and b2 != 'N':
            error_indices.append(i+1+truth_pos)
    return error_indices


def read_snps(snps_file, skip_indels):
    # snps_file format:
    # [P1]  [SUB]   [P2]    [BUFF]  [DIST]  [FRM]   [TAGS]
    # 13216	G	A	13216	261	13216	1	1	chr22.0.15649320.35205115	chr22.0.17366122.33487568
    # 13477	.	G	13478	31	13477	1	1	chr22.0.15649320.35205115	chr22.0.17366122.33487568
    hap2snplist = {}
    with open(snps_file, 'r') as f:
        for line in f:
            line = line.rstrip().split()
            [pos1, base1, base2, pos2] = line[:4]
            [hap1, hap2] = line[-2:]
            if skip_indels and '.' in [base1, base2]:
                continue
            try:
                hap2snplist[hap1].append(int(pos1))
                hap2snplist[hap2].append(int(pos2))
            except KeyError:
                hap2snplist[hap1] = [int(pos1)]
                hap2snplist[hap2] = [int(pos2)]
    return hap2snplist


def read_sam(filename):
    # returns a list of all alignments, where each alignment is presented in the
    # following format: [seq_id, ref_id, pos, ori, cigar]
    print("Reading SAM-file")
    aln_list = []
    unmapped_ids = []
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == "@":
                continue
            line = line.rstrip('\n').split()
            [seq_id, flag, ref_id, pos] = line[0:4]
            cigar = line[5]
            bits_flag = power_find(int(flag))
            if 4 in bits_flag: # unmapped contig
                unmapped_ids.append(seq_id)
                continue
            elif 16 in bits_flag: # reversed contig
                ori = "-"
            else:
                ori = "+"
            aln_list.append([seq_id, ref_id, int(pos)-1, ori, cigar]) # SAM-format is 1-based; switching to 0-based here
    return [aln_list, unmapped_ids]


def read_fasta(filename, read_ab=False):
    # returns ID to sequence dict
    id2seq = {}
    ab_est = {}
    with open(filename, 'r') as f:
        seq_id = ""
        seq = ""
        for line in f:
            if line[0] == '>':
                if seq_id != "" and seq != "":
                    id2seq[seq_id] = seq
                    if read_ab:
                        ab_est[seq_id] = ab # estimated abundance
                seq_id = line.lstrip('>').rstrip('\n').split()[0]
                if read_ab:
                    try:
                        ab = float(line.lstrip('>').rstrip('\n').split()[-1].lstrip('frequency='))
                    except ValueError:
                        print("WARNING: could not read abundance estimates from fasta")
                        read_ab = False
                seq = ""
            else:
                seq += line.rstrip('\n')
        # add final entry
        if seq_id != "" and seq != "":
            id2seq[seq_id] = seq
            if read_ab:
                ab_est[seq_id] = ab # estimated abundance
    return id2seq, ab_est


def power_find(n):
    try:
        int(n)
    except TypeError:
        print("power_find TypeError")
        print("n = {}".format(n))
        sys.exit(1)
    result = []
    binary = bin(n)[:1:-1]
    for x in range(len(binary)):
        if int(binary[x]):
            result.append(2**x)
    return result


def revcomp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    revcomp = "".join(complement.get(base, base) for base in reversed(seq))
    assert len(seq) == len(revcomp)
    return revcomp


if __name__ == '__main__':
    sys.exit(main())
