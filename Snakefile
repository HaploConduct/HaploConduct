from __future__ import division
import os
import sys

### SNAKEFILE FOR RUNNING SAVAGE ###

# path to config file
configfile: "savage_config.yaml"

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

def analyze_input(filename):
    total_len = 0
    longest_seq = 0
    i = 0
    with open(filename) as f:
        for line in f:
            if i%4 == 1: # sequence line in fastq
                l = len(line.strip('\n'))
                total_len += l
                if l > longest_seq:
                    longest_seq = l
            i += 1
        assert i % 4 == 0 # fastq
    seq_count = i/4
    return [seq_count, total_len, longest_seq]

#------------------

ALL_TYPES = ["singles", "paired1", "paired2"]
SPLIT_RANGE = range(config["SPLIT_NUM"])
PATH = os.getcwd()

if config["SINGLES_FILE"] != "":
    [seq_count, total_len, longest_seq] = analyze_input(config["SINGLES_FILE"])
else:
    [seq_count_1, total_len_1, longest_seq_1] = analyze_input(config["PAIRED1_FILE"])
    total_len = 2*total_len_1 # approximately
    longest_seq = 2*longest_seq_1 # approximately
AV_READ_LEN = total_len/seq_count
SAVAGE = config["SAVAGE_EXE"] + " --average_read_len %s" % AV_READ_LEN

if config["SINGLES_FILE"] != "" and config["PAIRED1_FILE"] != "" and config["PAIRED2_FILE"] != "":
    INPUT_FILES = [config["SINGLES_FILE"], config["PAIRED1_FILE"], config["PAIRED2_FILE"]]
    savage_input = " -s " + config["SINGLES_FILE"] + " -p1 " + config["PAIRED1_FILE"] + " -p2 " + config["PAIRED2_FILE"]
elif config["SINGLES_FILE"] != "":
    INPUT_FILES = [config["SINGLES_FILE"]]
    savage_input = " -s " + config["SINGLES_FILE"]
else:
    INPUT_FILES = [config["PAIRED1_FILE"], config["PAIRED2_FILE"]]
    savage_input = " -p1 " + config["PAIRED1_FILE"] + " -p2 " + config["PAIRED2_FILE"]

OVERLAPS_INPUT = expand("stage_a/patch{i}/input_fas/{type}.fastq", i=SPLIT_RANGE, type=ALL_TYPES)

if config["DE_NOVO"]:
    PREP_INPUT = INPUT_FILES
else:
    PREP_INPUT = INPUT_FILES + [config["REF"]]
    OVERLAPS_INPUT += expand("stage_a/patch{i}/input_fas/singles.sam", i=SPLIT_RANGE)
    OVERLAPS_INPUT += expand("stage_a/patch{i}/input_fas/paired.sam", i=SPLIT_RANGE)

#------------------

rule all:
    input:
        expand("contigs_stage_{x}.fasta", x=["a", "b", "c"]),
        expand("frequencies_stage_{x}.txt", x=["a", "b", "c"])

# prepares all input files for SAVAGE: splitting the data and aligning if necessary
rule savage_preprocessing:
    input:
        expand("{input}", input=PREP_INPUT)
    output:
        expand("{input}", input=OVERLAPS_INPUT)
    params:
        de_novo=config["DE_NOVO"],
        min_overlap_len=config["MIN_OVERLAP_LEN"],
        ref=config["REF"],
        split_num=config["SPLIT_NUM"],
        savage_exe=SAVAGE + savage_input
    run:
        if params.de_novo: # de novo mode
            shell("%s --split %s --no_overlaps --no_stage_a --no_stage_b --no_stage_c --min_overlap_len %s --num_threads %s" % (params.savage_exe, params.split_num, params.min_overlap_len, config["NUM_THREADS"]))
        else: # ref-guided paired-end reads
            shell("%s --split %s --no_overlaps --no_stage_a --no_stage_b --no_stage_c --min_overlap_len %s --num_threads %s --ref %s" % (params.savage_exe, params.split_num, params.min_overlap_len, config["NUM_THREADS"], params.ref))

# compute all overlaps for overlap graph construction
rule savage_overlaps:
    input:
        expand("{input}", input=OVERLAPS_INPUT)
    output:
        expand("stage_a/patch{i}/original_overlaps.txt", i=SPLIT_RANGE),
        "benchmarks/stage_a_overlaps.time.txt"
    params:
        de_novo=config["DE_NOVO"],
        min_overlap_len=config["MIN_OVERLAP_LEN"],
        ref=config["REF"],
        split_num=config["SPLIT_NUM"],
        savage_exe=SAVAGE
    benchmark:
        "benchmarks/stage_a_overlaps.txt"
    run:
        if params.de_novo: # de novo mode
            shell("/usr/bin/time -v -o benchmarks/stage_a_overlaps.time.txt %s --split %s --no_preprocessing --no_stage_a --no_stage_b --no_stage_c --min_overlap_len %s --num_threads %s" % (params.savage_exe, params.split_num, params.min_overlap_len, config["NUM_THREADS"]))
        else: # ref-guided paired-end reads
            shell("/usr/bin/time -v -o benchmarks/stage_a_overlaps.time.txt %s --split %s --no_preprocessing --no_stage_a --no_stage_b --no_stage_c --min_overlap_len %s --num_threads %s --ref %s" % (params.savage_exe, params.split_num, params.min_overlap_len, config["NUM_THREADS"], params.ref))

# run SAVAGE stage a: error correction and initial contig formation
rule savage_stage_a:
    input:
        expand("{input}", input=PREP_INPUT),
        expand("stage_a/patch{i}/input_fas/{type}.fastq", i=SPLIT_RANGE, type=ALL_TYPES),
        expand("stage_a/patch{i}/original_overlaps.txt", i=SPLIT_RANGE)
    output:
        fasta="contigs_stage_a.fasta",
        fastq="stage_a/singles.fastq",
        subreads="stage_a/subreads.txt",
        time="benchmarks/stage_a_main.time.txt"
    params:
        de_novo=config["DE_NOVO"],
        min_overlap_len=config["MIN_OVERLAP_LEN"],
        ref=config["REF"],
        split_num=config["SPLIT_NUM"],
        savage_exe=SAVAGE
    benchmark:
        "benchmarks/stage_a_main.txt"
    run:
        if config["REMOVE_BRANCHES"] == 0:
            shell("/usr/bin/time -v -o benchmarks/stage_a_main.time.txt %s --split %s --no_preprocessing --no_overlaps --no_stage_b --no_stage_c --min_overlap_len %s --num_threads %s --keep_branches" % (params.savage_exe, params.split_num, params.min_overlap_len, config["NUM_THREADS"]))
        elif config["REMOVE_BRANCHES"] == 1:
            shell("/usr/bin/time -v -o benchmarks/stage_a_main.time.txt %s --split %s --no_preprocessing --no_overlaps --no_stage_b --no_stage_c --min_overlap_len %s --num_threads %s" % (params.savage_exe, params.split_num, params.min_overlap_len, config["NUM_THREADS"]))
        else:
            print("REMOVE_BRANCHES must be either 0 or 1")

# run SAVAGE stage b: iterative contig extension
rule savage_stage_b:
    input:
        "stage_a/singles.fastq",
        "stage_a/subreads.txt",
        "contigs_stage_a.fasta"
    output:
        "stage_b/singles.fastq",
        "stage_b/subreads.txt",
        "contigs_stage_b.fasta",
        "benchmarks/stage_b.time.txt"
    params:
        de_novo=config["DE_NOVO"],
        min_overlap_len=config["MIN_OVERLAP_LEN"],
        ref=config["REF"],
        split_num=config["SPLIT_NUM"],
        savage_exe=SAVAGE
    benchmark:
        "benchmarks/stage_b.txt"
    run:
        if config["REMOVE_BRANCHES"] == 0:
            shell("/usr/bin/time -v -o benchmarks/stage_b.time.txt %s --split %s --no_overlaps --no_stage_a --no_stage_c --min_overlap_len %s --num_threads %s --keep_branches" % (params.savage_exe, params.split_num, params.min_overlap_len, config["NUM_THREADS"]))
        elif config["REMOVE_BRANCHES"] == 1:
            shell("/usr/bin/time -v -o benchmarks/stage_b.time.txt %s --split %s --no_overlaps --no_stage_a --no_stage_c --min_overlap_len %s --num_threads %s" % (params.savage_exe, params.split_num, params.min_overlap_len, config["NUM_THREADS"]))
        else:
            print("REMOVE_BRANCHES must be either 0 or 1")

# run SAVAGE stage c: master strain assembly
rule savage_stage_c:
    input:
        "stage_b/singles.fastq",
        "stage_b/subreads.txt",
        "contigs_stage_b.fasta"
    output:
        "stage_c/singles.fastq",
        "stage_c/subreads.txt",
        "contigs_stage_c.fasta",
        "benchmarks/stage_c.time.txt"
    params:
        de_novo=config["DE_NOVO"],
        min_overlap_len=config["MIN_OVERLAP_LEN"],
        ref=config["REF"],
        split_num=config["SPLIT_NUM"],
        savage_exe=SAVAGE,
        merge_contigs=config["MERGE_CONTIGS"]
    benchmark:
        "benchmarks/stage_c.txt"
    run:
        if config["REMOVE_BRANCHES"] == 0:
            shell("/usr/bin/time -v -o benchmarks/stage_c.time.txt %s --split %s --no_overlaps --no_stage_a --no_stage_b --min_overlap_len %s --num_threads %s --merge_contigs %s --keep_branches" % (params.savage_exe, params.split_num, params.min_overlap_len, config["NUM_THREADS"], params.merge_contigs))
        elif config["REMOVE_BRANCHES"] == 1:
            shell("/usr/bin/time -v -o benchmarks/stage_c.time.txt %s --split %s --no_overlaps --no_stage_a --no_stage_b --min_overlap_len %s --num_threads %s --merge_contigs %s" % (params.savage_exe, params.split_num, params.min_overlap_len, config["NUM_THREADS"], params.merge_contigs))
        else:
            print("REMOVE_BRANCHES must be either 0 or 1")

# apply SAVAGE frequency estimation (quick-mode)
rule frequency_estimation:
    input:
        fastq="stage_{x}/singles.fastq",
        fasta="contigs_stage_{x}.fasta",
        subreads="stage_{x}/subreads.txt"
    params:
        kallisto=config["KALLISTO"],
        kallisto_path=config["KALLISTO_PATH"],
        savage_exe=config["SAVAGE_EXE"],
        min_overlap_len=config["MIN_OVERLAP_LEN"],
        min_len=config["MIN_LEN_FREQ_EST"]
    output:
        "frequencies_stage_{x, [a-c]}.txt"
    run:
        correction = max(0, 500 - 2*params.min_overlap_len) # assuming a fragment size of 500
        freq_est_exe = '/'.join((params.savage_exe).split('/')[:-1] + ["freq_est"])
        if params.kallisto == 1:
            shell("%s --kallisto -c {input.fasta} -m {params.min_len} -f %s -r %s --kallisto_path {params.kallisto_path} -o {output}" % (freq_est_exe, config["FORWARD"], config["REVERSE"]))
        else:
            shell("%s -c {input.fastq} -s {input.subreads} -m {params.min_len} -k %s -o {output}" % (freq_est_exe, correction))


# rule prepare_input_fastq:
#     input:
#         expand("{filename}", filename=INPUT_FILES)
#     params:
#         savage_dir=config["SAVAGE_DIR"]
# #        singles=config["SINGLE_END_READS"]
#     output:
#         expand("input_fas/{type}.fastq", type=ALL_TYPES)
#     run:
#         shell("mkdir -p input_fas")
#         for fastq_type in ALL_TYPES:
#             shell("rm -f input_fas/%s.fastq" % fastq_type)
#             shell("touch input_fas/%s.fastq" % fastq_type)
#         if len(INPUT_FILES) == 1 or len(INPUT_FILES) == 3:
#             singles_count = int(file_len(input[0])/4)
#             shell("python2 {params.savage_dir}/scripts/rename_fas.py "
#                    "--in {input[0]} --out {output[0]}")
#             if len(INPUT_FILES) == 3:
#                 shell("python2 {params.savage_dir}/scripts/rename_fas.py "
#                        "--in {input[1]} --out {output[1]} "
#                        "--id_start %d" % singles_count)
#                 shell("python2 {params.savage_dir}/scripts/rename_fas.py "
#                        "--in {input[2]} --out {output[2]} "
#                        "--id_start %d" % singles_count)
#         else:
#             singles_count = 0
#             shell("python2 {params.savage_dir}/scripts/rename_fas.py "
#                    "--in {input[0]} --out {output[1]} "
#                    "--id_start %d" % singles_count)
#             shell("python2 {params.savage_dir}/scripts/rename_fas.py "
#                    "--in {input[1]} --out {output[2]} "
#                    "--id_start %d" % singles_count)
#
#
# rule split_fastqs:
#     input:
#         fastq=expand("input_fas/{type}.fastq", type=ALL_TYPES)
#     params:
#         split_num=str(config["SPLIT_NUM"]),
#         savage_dir=config["SAVAGE_DIR"]
#     output:
#         expand("stage_a/{type}.{i}.fastq", type=ALL_TYPES, i=SPLIT_RANGE)
#     run:
#         shell("mkdir -p stage_a")
#         shell("python2 {params.savage_dir}/scripts/random_split_fastq.py "
#                 "--input input_fas/singles.fastq "
#                 "--output stage_a/singles "
#                 "--split_num {params.split_num};")
#         shell("python2 {params.savage_dir}/scripts/random_split_fastq.py "
#                     "--input input_fas/paired1.fastq "
#                     "--input2 input_fas/paired2.fastq "
#                     "--output stage_a/paired "
#                     "--split_num {params.split_num};")
#
# rule rename_fastq:
#     input:
#         expand("{{path}}/{file}.{{i}}.fastq", file=ALL_TYPES)
#     params:
#         savage_dir=config["SAVAGE_DIR"]
#     output:
#         expand("{{path}}/renamed.{file}.{{i, [0-9]+}}.fastq", file=ALL_TYPES)
#     run:
#         singles_count = int(file_len(input[0])/4)
#         shell("python2 {params.savage_dir}/scripts/rename_fas.py "
#                    "--in {input[0]} --out {output[0]}")
#         shell("python2 {params.savage_dir}/scripts/rename_fas.py "
#                "--in {input[1]} --out {output[1]} "
#                "--id_start %d" % singles_count)
#         shell("python2 {params.savage_dir}/scripts/rename_fas.py "
#                "--in {input[2]} --out {output[2]} "
#                "--id_start %d" % singles_count)
#
# rule create_patch_dir:
#     input:
#         expand("stage_a/renamed.{type}.{{i}}.fastq", type=ALL_TYPES)
#     output:
#         expand("stage_a/patch{{i, [0-9]+}}/input_fas/{type}.fastq", type=ALL_TYPES),
#         "stage_a/patch{i}/input_fas/singles.sam" if not config["DE_NOVO"] else "stage_a/patch{i}/input_fas/singles.fastq",
#         "stage_a/patch{i}/input_fas/paired.sam" if not config["DE_NOVO"] else "stage_a/patch{i}/input_fas/paired1.fastq"
#     params:
#         de_novo=config["DE_NOVO"],
#         ref=config["REF"]
#     run:
#         shell("mkdir -p stage_a/patch{wildcards.i}/input_fas")
#         for fastq_type in ALL_TYPES:
#             shell("mv stage_a/renamed.%s.{wildcards.i}.fastq "
#                 "stage_a/patch{wildcards.i}/input_fas/%s.fastq"
#                 % (fastq_type, fastq_type))
#         if not params.de_novo:
#             shell("bwa index %s > /dev/null 2>&1" % params.ref)
#             path = "stage_a/patch{wildcards.i}"
#             shell("bwa mem %s %s/input_fas/singles.fastq 1> %s/input_fas/singles.sam 2> /dev/null" % (params.ref, path, path))
#             shell("bwa mem %s %s/input_fas/paired1.fastq %s/input_fas/paired2.fastq 1> %s/input_fas/paired.sam 2> /dev/null" % (params.ref, path, path, path))
#
#
# rule combine_contigs:
#     input:
# #        expand("stage_a/patch{i}/contigs_stage_a.fasta", i=SPLIT_RANGE),
#         expand("stage_a/patch{i}/stage_a/singles.fastq", i=SPLIT_RANGE),
#         expand("stage_a/patch{i}/stage_a/subreads.txt", i=SPLIT_RANGE)
#     output:
#         "stage_a/singles.fastq",
#         "stage_a/subreads.txt",
#         "contigs_stage_a.fasta"
#     params:
#         savage_dir=config["SAVAGE_DIR"]
#     run:
#         shell("rm -f stage_a/combined_singles.fastq")
#         shell("rm -f stage_a/subreads.txt")
#         total_contigs = 0
#         with open('stage_a/subreads.txt', 'w') as new_subreads_file:
#             for patch_num in SPLIT_RANGE:
#                 shell("cat stage_a/patch%d/stage_a/singles.fastq >> stage_a/combined_singles.fastq" % patch_num)
#                 singles_count = round(file_len('stage_a/patch%d/stage_a/singles.fastq' % patch_num)/4)
#                 renamed2originals = {}
#                 with open('stage_a/%s.%d.fastq' % (ALL_TYPES[0], patch_num), 'r') as f1:
#                     with open('stage_a/%s.%d.fastq' % (ALL_TYPES[1], patch_num), 'r') as f2:
#                         i = 0
#                         for line in f1:
#                             if i%4 == 0:
#                                 old_id = line.strip('\n')[1:]
#                                 new_id = str(round(i/4))
#                                 renamed2originals[new_id] = old_id
#                             i += 1
#                         assert i%4 == 0
#                         for line in f2:
#                             if i%4 == 0:
#                                 old_id = line.strip('\n')[1:]
#                                 new_id = str(round(i/4))
#                                 renamed2originals[new_id] = old_id
#                             i += 1
#                 with open('stage_a/patch%d/stage_a/subreads.txt' % patch_num, 'r') as f3:
#                     for line in f3:
#                         split_line = line.strip('\n').split('\t')
#                         contig_id = split_line[0]
#                         if int(contig_id) >= singles_count: # paired-end contig
#                             continue
#                         new_contig_id = str(int(contig_id) + total_contigs)
#                         new_line = [new_contig_id]
#                         for subread_info in split_line[1:]:
#                             [ID, poslist] = subread_info.split(':')
#                             new_subread_id = renamed2originals[ID]
#                             new_info = new_subread_id + ':' + poslist
#                             new_line.append(new_info)
#                         new_subreads_file.write('\t'.join(new_line) + '\n')
#                 total_contigs += singles_count
#         # now rename the merged fastq and convert to fasta
#         shell("python2 {params.savage_dir}/scripts/rename_fas.py --in stage_a/combined_singles.fastq --out stage_a/singles.fastq")
#         shell("python2 {params.savage_dir}/scripts/fastq2fasta.py stage_a/singles.fastq contigs_stage_a.fasta")
