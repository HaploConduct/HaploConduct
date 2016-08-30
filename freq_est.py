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

usage = """%prog [options]

Estimate relative frequencies for the input contigs.

"""

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument('--fas', dest='contigs_file', type=str, required=True)
    parser.add_argument('--subreads', dest='subreads_file', type=str, required=True)
    parser.add_argument('--split-subreads', dest='split_subreads_file', type=str, required=False)
#    parser.add_argument('--out', dest='out_file', type=str, required=True)
    parser.add_argument('--min_len', dest='min_len', type=int, default=0)
#    parser.add_argument('--min_output_len', dest='min_output_len', type=int, default=0)
#    parser.add_argument('--genome_len', dest='genome_len', type=int)
    args = parser.parse_args()
    
    # get read lengths from contig fastq
    contig_dict = {}
    with open(args.contigs_file, 'r') as f:
        k = 2 if args.contigs_file[-1] == "a" else 4
        i = 0
        for line in f:
            if i%k == 0:
                ID = line.strip('\n')[1:]
            elif i%k == 1:
                seq = line.strip('\n')
                if len(seq) >= args.min_len:
                    contig_dict[ID] = seq
            i += 1
    contig_count = len(contig_dict)
    total_len = 0
    for ID, seq in contig_dict.iteritems():
        total_len += len(seq)
    print "#contigs: ", contig_count
    if contig_count == 0:
        print "WARNING: NO CONTIGS OF SUFFICIENT LENGTH"
        average_len = 0
    else:
        average_len = total_len/contig_count
    print "total length: ", total_len
    print "average length: ", average_len
    
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
            
    # if the contigs result from split data, also use the split contig files
    originals2initials = {}
    initials2originals = {}
    testcount = 0
    if args.split_subreads_file:
        files = (args.split_subreads_file).split(',')
        contig_id = 0
        file_num = 0
        for splitfile in files:
            with open(splitfile, 'r') as f:
                initials = Set()
                for line in f:
                    line = line.strip('\n').split('\t')
                    subreads_info = line[1:]
                    subreads = []
                    for info in subreads_info:
                        [initial_ID, poslist] = info.split(':')
                        new_ID = str(int(initial_ID) + 200000*file_num)
                        if new_ID in initials2originals:
                            if len(initials2originals[new_ID]) > 0:
                                print initials2originals[new_ID]
                        initials2originals[new_ID] = []
                        initials.add(new_ID)
#                print "initials subcount %d" % len(initials)
                testcount += len(initials)
                    
            with open(splitfile, 'r') as f:
                for line in f:
                    line = line.strip('\n').split('\t')
                    subreads_info = line[1:]
                    subreads = []
                    for info in subreads_info:
                        [initial_ID, poslist] = info.split(':')
                        new_ID = str(int(initial_ID) + 200000*file_num)
                        initials2originals[new_ID] += [str(contig_id)]
                        subreads += [new_ID]
                    originals2initials[str(contig_id)] = subreads
                        
                    contig_id += 1
            file_num += 1
        
        print "testcount: %d" % testcount            
        print "largest contig id: %d" % contig_id
                    
        # update the contig dicts
        initials_used = [0 for i in xrange(100000000)]
        
        new_contigs2originals = {}
        for contig, subreads in contigs2originals.iteritems():
            initial_subreads = Set()
            for subread in subreads:
                initials = originals2initials[subread]
                for read in initials:
                    initial_subreads.add(read)
                    initials_used[int(read)] = 1
            assert len(initial_subreads) > 0
            new_contigs2originals[contig] = [read for read in initial_subreads]
            
        new_originals2contigs = {}
        for initial, originals in initials2originals.iteritems():
            contig_set = Set()
            for original in originals:
                if original in originals2contigs:
                    contigs = originals2contigs[original]
                    for contig in contigs:
                        contig_set.add(contig)
            if len(contig_set) > 0:
                new_originals2contigs[initial] = contig_set
            elif initials_used[int(initial)] == 1:
                print initial
        
        print "initials used: ", initials_used.count(1)    
        contigs2originals = new_contigs2originals
        originals2contigs = new_originals2contigs        
            
    # estimate contig frequencies
#    if args.genome_len:    
#        genome_len = args.genome_len
#    else:
#        genome_len = total_len
    total_subreads_used = 0
    for original, contigs in originals2contigs.iteritems():
        if len(contigs) > 0:
            total_subreads_used += 1
    print "total subread count: ", total_subreads_used
    tmp_freqs = []
    tmp_reads = []
    tmp_lengths = []
#    print "absolute frequencies:"
#    print "(ID: freq, length)"
    for read in contigs2originals:
        seq = contig_dict[read]
        weighted_count = 0
        for subread in contigs2originals[read]:
            if subread in originals2contigs:
                weighted_count += 1/len(originals2contigs[subread])
#        freq = (weighted_count/total_subreads_used)*(genome_len/len(seq))
        freq = (weighted_count/total_subreads_used)*(1/len(seq))
        if len(contig_dict[read]) > args.min_len:
#            print "%s: %f, %d" % (read, freq, len(seq))
            tmp_freqs += [freq]
            tmp_reads += [read]
            tmp_lengths += [len(seq)]
            
    print "*"
    
#    if len(tmp_freqs) == 2: # for div-vs-ratio experiments  
#        max_freq = max(tmp_freqs)/float(sum(tmp_freqs))*100
#        print "\t\t\t\tmax normalized frequency: %.2f" % max_freq
#        print "*"
#    else:
    print "normalized frequencies:"
    print "(ID: freq, length)"
    i = 0
    for freq in tmp_freqs:
        read = tmp_reads[i]
        length = tmp_lengths[i]
        normalized_freq = freq/float(sum(tmp_freqs))*100
        print "%s: %.2f, %d" % (read, normalized_freq, length)
        i += 1
    print "*"
        
        
if __name__ == '__main__':
    sys.exit(main())
