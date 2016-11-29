from __future__ import division
import os
import sys

# SNAKEFILE FOR SPLITTING DATA

# path to config file
configfile: "savage_config.yaml"

if not config["DE_NOVO"]:
    assert os.path.isfile(config["REF"])

# if config["SINGLE_END_READS"]:
#     FASTQ_TYPES = ["singles"]
#     INPUT_FILES = [config["SINGLES_FILE"]]
# else:
#     FASTQ_TYPES = ["paired1", "paired2"]
#     INPUT_FILES = [config["PAIRED1_FILE"], config["PAIRED2_FILE"]]

INPUT_FILES = [config["SINGLES_FILE"], config["PAIRED1_FILE"], config["PAIRED2_FILE"]]
ALL_TYPES = ["singles", "paired1", "paired2"]
SPLIT_RANGE = range(config["SPLIT_NUM"])
PATH = os.getcwd()
#ALIGNMENT = "singles.sam" if config["SINGLE_END_READS"] else "paired.sam"

#------------------

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

#------------------

rule all:
    input:
        expand("frequencies_stage_{x}.txt", x=["b", "c"]),
        expand("benchmarks/stage_a/patch{i}_preprocessing.time.txt", i=SPLIT_RANGE),
        expand("benchmarks/stage_a/patch{i}_main.time.txt", i=SPLIT_RANGE)

rule prepare_input_fastq:
    input:
        expand("{filename}", filename=INPUT_FILES)
    params:
        savage_dir=config["SAVAGE_DIR"]
#        singles=config["SINGLE_END_READS"]
    output:
        expand("input_fas/{type}.fastq", type=ALL_TYPES)
    run:
        shell("mkdir -p input_fas")
        for fastq_type in ALL_TYPES:
            shell("rm -f input_fas/%s.fastq" % fastq_type)
            shell("touch input_fas/%s.fastq" % fastq_type)
        singles_count = int(file_len(input[0])/4)
        if input[0] != "":
            shell("python2 {params.savage_dir}/scripts/rename_fas.py "
                   "--in {input[0]} --out {output[0]}")
        if input[1] != "" and input[2] != "":
            shell("python2 {params.savage_dir}/scripts/rename_fas.py "
                   "--in {input[1]} --out {output[1]} "
                   "--id_start %d" % singles_count)
            shell("python2 {params.savage_dir}/scripts/rename_fas.py "
                   "--in {input[2]} --out {output[2]} "
                   "--id_start %d" % singles_count)


rule split_fastqs:
    input:
        fastq=expand("input_fas/{type}.fastq", type=ALL_TYPES)
    params:
        split_num=str(config["SPLIT_NUM"]),
        savage_dir=config["SAVAGE_DIR"]
    output:
        expand("stage_a/{type}.{i}.fastq", type=ALL_TYPES, i=SPLIT_RANGE)
    run:
        shell("mkdir -p stage_a")
        shell("python2 {params.savage_dir}/scripts/random_split_fastq.py "
                "--input input_fas/singles.fastq "
                "--output stage_a/singles "
                "--split_num {params.split_num};")
        shell("python2 {params.savage_dir}/scripts/random_split_fastq.py "
                    "--input input_fas/paired1.fastq "
                    "--input2 input_fas/paired2.fastq "
                    "--output stage_a/paired "
                    "--split_num {params.split_num};")

rule rename_fastq:
    input:
        expand("{{path}}/{file}.{{i}}.fastq", file=ALL_TYPES)
    params:
        savage_dir=config["SAVAGE_DIR"]
    output:
        expand("{{path}}/renamed.{file}.{{i, [0-9]+}}.fastq", file=ALL_TYPES)
    run:
        singles_count = int(file_len(input[0])/4)
        shell("python2 {params.savage_dir}/scripts/rename_fas.py "
                   "--in {input[0]} --out {output[0]}")
        shell("python2 {params.savage_dir}/scripts/rename_fas.py "
               "--in {input[1]} --out {output[1]} "
               "--id_start %d" % singles_count)
        shell("python2 {params.savage_dir}/scripts/rename_fas.py "
               "--in {input[2]} --out {output[2]} "
               "--id_start %d" % singles_count)

rule create_patch_dir:
    input:
        expand("stage_a/renamed.{type}.{{i}}.fastq", type=ALL_TYPES)
    output:
        expand("stage_a/patch{{i, [0-9]+}}/input_fas/{type}.fastq", type=ALL_TYPES),
        "stage_a/patch{i}/input_fas/singles.sam" if not config["DE_NOVO"] else "stage_a/patch{i}/input_fas/singles.fastq",
        "stage_a/patch{i}/input_fas/paired.sam" if not config["DE_NOVO"] else "stage_a/patch{i}/input_fas/paired1.fastq"
    params:
        de_novo=config["DE_NOVO"],
        ref=config["REF"]
    run:
        shell("mkdir -p stage_a/patch{wildcards.i}/input_fas")
        for fastq_type in ALL_TYPES:
            shell("mv stage_a/renamed.%s.{wildcards.i}.fastq "
                "stage_a/patch{wildcards.i}/input_fas/%s.fastq"
                % (fastq_type, fastq_type))
        if not params.de_novo:
            shell("bwa index %s > /dev/null 2>&1" % params.ref)
            path = "stage_a/patch{wildcards.i}"
            shell("bwa mem %s %s/input_fas/singles.fastq 1> %s/input_fas/singles.sam 2> /dev/null" % (params.ref, path, path))
            shell("bwa mem %s %s/input_fas/paired1.fastq %s/input_fas/paired2.fastq 1> %s/input_fas/paired.sam 2> /dev/null" % (params.ref, path, path, path))

rule savage_preprocessing:
    input:
        expand("stage_a/patch{{i}}/input_fas/{type}.fastq", type=ALL_TYPES),
        "stage_a/patch{i}/input_fas/singles.sam" if not config["DE_NOVO"] else "stage_a/patch{i}/input_fas/singles.fastq",
        "stage_a/patch{i}/input_fas/paired.sam" if not config["DE_NOVO"] else "stage_a/patch{i}/input_fas/paired1.fastq"
    output:
        "stage_a/patch{i, [0-9]+}/original_overlaps.txt",
        "benchmarks/stage_a/patch{i, [0-9]+}_preprocessing.time.txt"
    params:
        de_novo=config["DE_NOVO"],
        savage_dir=config["SAVAGE_DIR"],
        min_overlap_len=config["MIN_OVERLAP_LEN"],
        ref=config["REF"]
    benchmark:
        "benchmarks/stage_a/patch{i}_preprocessing.txt"
    run:
        os.chdir("stage_a/patch%s" % wildcards.i)
        if params.de_novo: # de novo mode
            shell("/usr/bin/time -v -o %s/benchmarks/stage_a/patch{wildcards.i}_preprocessing.time.txt python2 %s/savage.py --preprocessing --no_stage_a --no_stage_b --no_stage_c --min_overlap_len %s --num_threads %s" % (PATH, params.savage_dir, params.min_overlap_len, config["NUM_THREADS"]))
        # elif params.singles: # ref-guided single-end reads
        #     shell("/usr/bin/time -v -o %s/benchmarks/stage_a/patch{wildcards.i}_preprocessing.time.txt python2 %s/savage.py --preprocessing --no_stage_a --no_stage_b --no_stage_c --min_overlap_len %s --num_threads %s --ref %s --singles input_fas/singles.sam" % (PATH, params.savage_dir, params.min_overlap_len, config["NUM_THREADS"], params.ref))
        else: # ref-guided paired-end reads
            shell("/usr/bin/time -v -o %s/benchmarks/stage_a/patch{wildcards.i}_preprocessing.time.txt python2 %s/savage.py --preprocessing --no_stage_a --no_stage_b --no_stage_c --min_overlap_len %s --num_threads %s --ref %s --singles input_fas/singles.sam --paired input_fas/paired.sam" % (PATH, params.savage_dir, params.min_overlap_len, config["NUM_THREADS"], params.ref))
        os.chdir('../..')

rule savage_stage_a:
    input:
        expand("stage_a/patch{{i}}/input_fas/{type}.fastq", type=ALL_TYPES),
        "stage_a/patch{i}/original_overlaps.txt"
    output:
        fasta="stage_a/patch{i, [0-9]+}/contigs_stage_a.fasta",
        fastq="stage_a/patch{i, [0-9]+}/stage_a/singles.fastq",
        subreads="stage_a/patch{i, [0-9]+}/stage_a/subreads.txt",
        time="benchmarks/stage_a/patch{i, [0-9]+}_main.time.txt"
    params:
        savage_dir=config["SAVAGE_DIR"],
        min_overlap_len=config["MIN_OVERLAP_LEN"],
        threads=config["NUM_THREADS"]
    benchmark:
        "benchmarks/stage_a/patch{i}_main.txt"
    run:
        os.chdir("stage_a/patch%s" % wildcards.i)
        shell("/usr/bin/time -v -o %s/benchmarks/stage_a/patch{wildcards.i}_main.time.txt python2 %s/savage.py --no_stage_b --no_stage_c --overlaps %s/stage_a/patch{wildcards.i}/original_overlaps.txt --min_overlap_len %s --num_threads %s" % (PATH, params.savage_dir, PATH, params.min_overlap_len, params.threads))
        os.chdir('../..')

rule combine_contigs:
    input:
#        expand("stage_a/patch{i}/contigs_stage_a.fasta", i=SPLIT_RANGE),
        expand("stage_a/patch{i}/stage_a/singles.fastq", i=SPLIT_RANGE),
        expand("stage_a/patch{i}/stage_a/subreads.txt", i=SPLIT_RANGE)
    output:
        "stage_a/singles.fastq",
        "stage_a/subreads.txt",
        "contigs_stage_a.fasta"
    params:
        savage_dir=config["SAVAGE_DIR"]
    run:
        shell("rm -f stage_a/combined_singles.fastq")
        shell("rm -f stage_a/subreads.txt")
        total_contigs = 0
        with open('stage_a/subreads.txt', 'w') as new_subreads_file:
            for patch_num in SPLIT_RANGE:
                shell("cat stage_a/patch%d/stage_a/singles.fastq >> stage_a/combined_singles.fastq" % patch_num)
                singles_count = round(file_len('stage_a/patch%d/stage_a/singles.fastq' % patch_num)/4)
                renamed2originals = {}
                with open('stage_a/%s.%d.fastq' % (ALL_TYPES[0], patch_num), 'r') as f1:
                    with open('stage_a/%s.%d.fastq' % (ALL_TYPES[1], patch_num), 'r') as f2:
                        i = 0
                        for line in f1:
                            if i%4 == 0:
                                old_id = line.strip('\n')[1:]
                                new_id = str(round(i/4))
                                renamed2originals[new_id] = old_id
                            i += 1
                        assert i%4 == 0
                        for line in f2:
                            if i%4 == 0:
                                old_id = line.strip('\n')[1:]
                                new_id = str(round(i/4))
                                renamed2originals[new_id] = old_id
                            i += 1
                with open('stage_a/patch%d/stage_a/subreads.txt' % patch_num, 'r') as f3:
                    for line in f3:
                        split_line = line.strip('\n').split('\t')
                        contig_id = split_line[0]
                        if int(contig_id) >= singles_count: # paired-end contig
                            continue
                        new_contig_id = str(int(contig_id) + total_contigs)
                        new_line = [new_contig_id]
                        for subread_info in split_line[1:]:
                            [ID, poslist] = subread_info.split(':')
                            new_subread_id = renamed2originals[ID]
                            new_info = new_subread_id + ':' + poslist
                            new_line.append(new_info)
                        new_subreads_file.write('\t'.join(new_line) + '\n')
                total_contigs += singles_count
        # now rename the merged fastq and convert to fasta
        shell("python2 {params.savage_dir}/scripts/rename_fas.py --in stage_a/combined_singles.fastq --out stage_a/singles.fastq")
        shell("python2 {params.savage_dir}/scripts/fastq2fasta.py stage_a/singles.fastq contigs_stage_a.fasta")

rule savage_stage_b:
    input:
        "stage_a/singles.fastq",
        "stage_a/subreads.txt",
        "contigs_stage_a.fasta"
    output:
        "stage_b/singles.fastq",
        "stage_b/subreads.txt",
        "contigs_stage_b.fasta",
        "stage_b/time.txt"
    params:
        savage_dir=config["SAVAGE_DIR"],
        min_overlap_len=config["MIN_OVERLAP_LEN"]
    benchmark:
        "benchmarks/stage_b.txt"
    run:
        shell("/usr/bin/time -v -o stage_b/time.txt python2 %s/savage.py --no_stage_a --no_stage_c --min_overlap_len %s --num_threads %s --use_subreads=%s --remove_branches=%s" % (params.savage_dir, params.min_overlap_len, config["NUM_THREADS"], config["USE_SUBREADS"], config["REMOVE_BRANCHES"]))

rule savage_stage_c:
    input:
        "stage_b/singles.fastq",
        "stage_b/subreads.txt",
        "contigs_stage_b.fasta"
    output:
        "stage_c/singles.fastq",
        "stage_c/subreads.txt",
        "contigs_stage_c.fasta",
        "stage_c/time.txt"
    params:
        savage_dir=config["SAVAGE_DIR"],
        min_overlap_len=config["MIN_OVERLAP_LEN"],
        merge_contigs=config["MERGE_CONTIGS"]
    benchmark:
        "benchmarks/stage_c.txt"
    run:
        shell("/usr/bin/time -v -o stage_c/time.txt python2 %s/savage.py --no_stage_a --no_stage_b --min_overlap_len %s --num_threads %s --merge_contigs %s --use_subreads=%s --remove_branches=%s" % (params.savage_dir, params.min_overlap_len, config["NUM_THREADS"], params.merge_contigs, config["USE_SUBREADS"], config["REMOVE_BRANCHES"]))

rule frequency_estimation:
    input:
        fastq="stage_{x}/singles.fastq",
        subreads="stage_{x}/subreads.txt"
    params:
        savage_dir=config["SAVAGE_DIR"],
        min_len=config["MIN_LEN_FREQ_EST"]
    output:
        "frequencies_stage_{x, [a-c]}.txt"
    shell:
        "python2 {params.savage_dir}/freq_est.py --fas {input.fastq} "
        "--subreads {input.subreads} --min_len {params.min_len} "
        "--out {output}"
