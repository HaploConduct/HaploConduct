#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random
import subprocess
from time import clock
from sets import Set


__author__ = "Jasmijn Baaijens"

usage = """

Estimate relative frequencies for contigs produced by SAVAGE.

The frequency estimation procedure has two modes:
1. Compute frequency estimates for given contigs by running Kallisto, then parse
Kallisto output (abundance.tsv) and compute frequency estimates from TPM counts.
2. Quick estimates based on the relative number of reads contributing to a
contig during SAVAGE assembly.

Option 1 obtains most accurate results, but requires Kallisto to be installed.

Note:
For option 1, please make sure that the input fasta contains all assembled
sequences; contigs of sufficient length will be selected after estimating TPM
counts (Kallisto).
For option 2, please make sure to run SAVAGE with the option --use_subreads for
better frequency estimates of stage b/c contigs.

"""

def print_string(filename, string):
    if filename == "":
        print string
    else:
        with open(filename, 'a') as f:
            f.write(string + '\n')

def main():
    parser = ArgumentParser(description=usage)
    general = parser.add_argument_group('general arguments')
    general.add_argument('-c', '--contigs', dest='contigs', type=str, required=True, help='fastq or fasta file containing contigs')
    general.add_argument('-o', '--out', dest='out_file', type=str, default="", help='write output to file; if not specified, output is written to stdout.')
    general.add_argument('-m', '--min_len', dest='min_len', type=int, default=0, help='only consider contigs at least this length')
    general.add_argument('--select_ids', dest='select_ids', type=str, help='select target contigs for frequency estimation by providing a comma-separated list of contig identifiers')
    kallisto_mode = parser.add_argument_group('kallisto-mode arguments')
    kallisto_mode.add_argument('--kallisto', dest='kallisto', action='store_true', help='use Kallisto for improved frequency estimation')
    kallisto_mode.add_argument('--kallisto_path', dest='kallisto_path', type=str, default='kallisto', help='path to kallisto executable; needed only if Kallisto path is not included in PATH')
    kallisto_mode.add_argument('-l', '--fragmentsize', dest='fragmentsize', type=float, required=True, help='Estimated average fragment size')
    kallisto_mode.add_argument('-d', '--stddev', dest='stddev', type=float, required=True, help='Estimated standard deviation of fragment size')
    kallisto_mode.add_argument('-f', '--forward', dest='forward', type=str, help='original forward input reads (before assembly) in fastq format')
    kallisto_mode.add_argument('-r', '--reverse', dest='reverse', type=str, help='original reverse input reads (before assembly) in fastq format')
    quick_mode = parser.add_argument_group('quick-mode arguments')
    quick_mode.add_argument('-s', '--subreads', dest='subreads_file', type=str, help='subreads file produced by savage, corresponding to the fastq file provided')
    quick_mode.add_argument('-k', '--correction', dest='len_correction', type=float, default=0, help='(optional) correction term for computing effective contig length: use c = fragment_size - 2*(read_len-min_overlap_len)')
    args = parser.parse_args()

    if args.out_file != "" and os.path.isfile(args.out_file):
        os.remove(args.out_file)

    if args.kallisto and not args.forward:
        sys.stderr.write("Input arguments missing: Kallisto-mode requires original fastq files.\n")
        sys.stderr.flush()
        print "Please specify forward and reverse fastq files using options -f and -r, or single-end input using only -f"
        sys.exit(1)

    if not args.kallisto and not args.subreads_file:
        print "Input arguments missing: quick-mode requires SAVAGE subread files."
        print "Please specify the corresponding subread file using --subreads"
        sys.exit(1)

    if args.select_ids:
        select_contigs = args.select_ids.split(',')
    else:
        select_contigs = []

    if args.kallisto:
        # run kallisto on contig file
        kallisto_file = run_kallisto(args.kallisto_path, args.contigs, str(args.fragmentsize), str(args.stddev), args.forward, args.reverse)

        # filter for minimum length and translate kallisto TPM counts to frequencies
        freq_dict = process_kallisto_output(kallisto_file, args.min_len, select_contigs)
        # sort dict by estimated frequency
        sorted_freqs = sorted(freq_dict.items(), key=lambda x:x[1][0], reverse=True)
        # write the resulting frequencies to the output file
        print_string(args.out_file, "id\tlength\tfrequency")
        for contig_id, info in sorted_freqs:
            freq = info[0]
            length = info[1]
            print_string(args.out_file, contig_id + "\t%s\t%.3f" % (length, freq))
        print "*** Done ***\n"

    else:
        # Run quick-mode frequency estimation
        # get read lengths from contig fastq
        contig_dict = {}
        with open(args.contigs, 'r') as f:
            k = 2 if args.contigs[-1] == "a" else 4
            i = 0
            for line in f:
                if i%k == 0:
                    ID = line.strip('\n')[1:]
                elif i%k == 1:
                    seq = line.strip('\n')
                    if select_contigs == [] or ID in select_contigs:
                        if len(seq) >= args.min_len:
                            contig_dict[ID] = seq
                i += 1
        contig_count = len(contig_dict)
        total_len = 0
        for ID, seq in contig_dict.iteritems():
            total_len += len(seq)
        print_string(args.out_file, "#contigs: %d" % contig_count)
        eff_total_len = total_len + contig_count * (1-args.len_correction)
        if contig_count == 0:
            print_string(args.out_file, "WARNING: NO CONTIGS OF SUFFICIENT LENGTH")
            average_len = 0
            eff_average_len = 0
        else:
            average_len = total_len/contig_count
            eff_average_len = eff_total_len/contig_count
        print_string(args.out_file, "total length: %d" % total_len)
        print_string(args.out_file, "total effective length: %d" % eff_total_len)
        print_string(args.out_file, "average length: %d" % average_len)
        print_string(args.out_file, "average effective length: %d" % eff_average_len)

        # create a dict mapping original reads to contigs and vice versa
        contigs2originals = {}
        originals2contigs = {}

        with open(args.subreads_file, 'r') as f:
            for line in f:
                line = line.strip('\n').split('\t')
                contig = line[0]
                if contig not in contig_dict:
                    continue
                subreads_info = line[1:]
                subreads = []
                for info in subreads_info:
                    [ID, poslist] = info.split(':')
                    originals2contigs[ID] = []

        with open(args.subreads_file, 'r') as f:
            for line in f:
                line = line.strip('\n').split('\t')
                contig = line[0]
                if contig not in contig_dict:
                    continue
                subreads_info = line[1:]
                subreads = []
                for info in subreads_info:
                    [ID, poslist] = info.split(':')
                    originals2contigs[ID] += [contig]
                    subreads += [ID]
                contigs2originals[contig] = subreads

        total_subreads_used = 0
        for original, contigs in originals2contigs.iteritems():
            if len(contigs) > 0:
                total_subreads_used += 1
        print_string(args.out_file, "total subread count: %d" % total_subreads_used)
        tmp_freqs = []
        tmp_reads = []
        tmp_lengths = []
        tmp_eff_lengths = []
        for read in contigs2originals:
            seq = contig_dict[read]
            weighted_count = 0
            for subread in contigs2originals[read]:
                if subread in originals2contigs:
                    weighted_count += 1/len(originals2contigs[subread])
            eff_len = len(seq) - args.len_correction + 1
            freq = (weighted_count/total_subreads_used)*(1/eff_len)
            if len(contig_dict[read]) > args.min_len:
                tmp_freqs += [freq]
                tmp_reads += [read]
                tmp_lengths += [len(seq)]
                tmp_eff_lengths += [eff_len]

        print_string(args.out_file, "id\tlength\tfrequency")
        i = 0
        for freq in tmp_freqs:
            read = tmp_reads[i]
            length = tmp_lengths[i]
            eff_len = tmp_eff_lengths[i]
            normalized_freq = freq/float(sum(tmp_freqs))
            print_string(args.out_file, "%s\t%d\t%.3f" % (read, length, normalized_freq))
            i += 1


def process_kallisto_output(abundance_file, min_len, select_contigs):
    print "*** Processing Kallisto output ***\n"
    data = []
    id_list = []
    len_list = []
    with open(abundance_file, 'r') as f:
        i = 0
        for line in f:
            splitline = line.strip('\n').split('\t')
            if i == 0:
                id_idx = splitline.index('target_id')
                len_idx = splitline.index('length')
                tmp_idx = splitline.index('tpm')
            else:
                length = splitline[len_idx]
                tpm = splitline[tmp_idx]
                contig_id = splitline[id_idx]
                if select_contigs != [] and contig_id not in select_contigs:
                    continue
                elif int(length) > min_len:
                    data.append(tpm)
                    id_list.append(contig_id)
                    len_list.append(length)
            i += 1
    tpm_list = [float(x) for x in data]
    total_tpm = sum(tpm_list)
    freq_dict = {}
    for i in xrange(len(tpm_list)):
        tpm = tpm_list[i]
        contig_id = id_list[i]
        freq = tpm/total_tpm
        length = len_list[i]
        freq_dict[contig_id] = [freq, length]
    return freq_dict


def run_kallisto(kallisto, contigs, fragmentsize, stddev, forward, reverse=""):
    # create output directory
    subprocess.call(['mkdir', '-p', 'frequencies'])
    # index construction
    contigs_name, extension = os.path.splitext(contigs)
    index_file = 'frequencies/' + contigs_name + '.idx'
    print "\n*** Running Kallisto index construction ***"
    subprocess.check_call([kallisto, 'index', '-i', index_file, contigs])
    # estimate abundances
    print "*** Running Kallisto abundance quantification ***"
    if reverse:
        subprocess.check_call([kallisto, 'quant', '-i', index_file, '-o', 'frequencies/kallisto', '-b', '100', '-l', fragmentsize, '-s', stddev, forward, reverse])
    else:
        subprocess.check_call([kallisto, 'quant', '-i', index_file, '-o', 'frequencies/kallisto', '-b', '100', '-l', fragmentsize, '-s', stddev, '--single', forward])
    abundance_file = 'frequencies/kallisto/abundance.tsv'
    return abundance_file


if __name__ == '__main__':
    sys.exit(main())
