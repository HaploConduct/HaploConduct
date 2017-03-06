# SAVAGE: Strain Aware VirAl GEnome assembly

[![license](https://img.shields.io/badge/license-GPL%20v3.0-blue.svg)](http://www.gnu.org/licenses/)

## Synopsis

SAVAGE is a computational tool for reconstructing individual
haplotypes of intra-host virus strains (a viral quasispecies) without
the need for a high quality reference genome. SAVAGE makes use of
either FM-index based data structures or ad-hoc consensus reference
sequence for constructing overlap graphs from patient sample data.
In this overlap graph, nodes represent reads and/or contigs, while
edges reflect that two reads/contigs, based on sound statistical
considerations, represent identical haplotypic sequence.
Following an iterative scheme, a new overlap assembly algorithm that
is based on the enumeration of statistically well-calibrated groups
of reads/contigs then efficiently reconstructs the individual
haplotypes from this overlap graph.

For a more detailed description, please see the [SAVAGE preprint](http://biorxiv.org/content/early/2017/01/21/080341).

## Motivation

A viral quasispecies, the ensemble of viral strains populating an
infected person, can be highly diverse. For optimal assessment of
virulence, pathogenesis and therapy selection, determining the
haplotypes of the individual strains can play a key role. As many
viruses are subject to high mutation and recombination rates, high-
quality reference genomes are often not available at the time of a
new disease outbreak. SAVAGE is able to reconstruct individual
haplotypes of the viral quasispecies without the need for a high-
quality reference genome.

## Dependencies

The SAVAGE program is a combination of C++ code (which has to be
compiled) and Python scripts (version 2.6 or later, but not Python
3.x). The C++ part requires several boost libraries (boost::timer,
boost::system, and boost::program_options) and it needs a compiler
that supports [OpenMP](http://openmp.org/wp/) (such as `g++`).

For maximal clique enumeration, SAVAGE depends on the [quick-cliques package](https://github.com/darrenstrash/quick-cliques) which is already included.

For suffix-prefix overlaps (de novo mode), SAVAGE depends on the
[SFO package](https://www.cs.helsinki.fi/group/suds/sfo/) which is
also included.

For reference-guided assembly, SAVAGE depends on the [bwa mem](http://bio-bwa.sourceforge.net/) aligner. Please make sure this tool is installed.

Stages b and c of the algorithm also require `blastn` and `makeblastdb`
from the ncbi C++ Toolkit.

Optionally, the preprocessing of the input data and the different
algorithm stages can all be run in a single command using [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) to execute the
Snakefile from the SAVAGE root directory.

We provide a simple algorithm for relative frequency estimation of the assembled
contigs. For improved frequency estimation, we recommend using [Kallisto](https://pachterlab.github.io/kallisto/).
This is optional, for more information please see the frequency estimation
section below.

## Installation

Please download and install:

* [ncbi-blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [bwa mem](http://bio-bwa.sourceforge.net/)
* [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) (optional)
* [Kallisto](https://pachterlab.github.io/kallisto/) (optional)

*Each of these tools can also be installed using [Bioconda](https://bioconda.github.io/),
a distribution of bioinformatics software realized as a channel for the
versatile Conda package manager.*

Then download the [latest](https://bitbucket.org/jbaaijens/savage/downloads?tab=tags) version of **SAVAGE**, enter the repository and type `make`.

## Manual

Currently, SAVAGE is designed for Illumina MiSeq sequencing reads
that have been pre-processed by [PEAR](http://sco.h-its.org/exelixis/web/software/pear/) (for merging self-overlapping
pairs). Thus, the input consists of (a combination of) single-end and/or paired-end reads.
**Note:** SAVAGE expects paired-end reads stored in *two fastq files with both
sequences in forward orientation* (instead of the more standard forward-reverse)
because this how PEAR outputs non-merged reads. Please find a more detailed
description of the SAVAGE input files below.

SAVAGE offers two modes: SAVAGE-ref and SAVAGE-de-novo. SAVAGE-ref is the
more efficient option, while SAVAGE-de-novo gives (slightly) better results.
Instructions for using these modes can be found below.

The most important parameter to set is the minimal overlap length *M* during
overlap graph construction. For best results, please read the section on
parameters carefully.

### Input

SAVAGE takes as input single-end and/or paired-end reads, which should be stored
in separate fastq files (e.g. singles.fastq, paired1.fastq, paired2.fastq). The paired-end reads are assumed to be in *two fastq files with both
sequences in forward orientation* (instead of the more standard forward-reverse)
because this how PEAR outputs non-merged reads. These files can be specified
using the following program options:

* `-s` or `--input_s`
Input fastq file containing single-end reads.
* `-p1` or `--input_p1`
Input fastq file containing /1 paired-end reads.
* `-p2` or `--input_p2`
Input fastq file containing /2 paired-end reads.
* `--ref` (optional)
Fasta file containing the reference genome to be used for reference-guided mode.
*Please make sure to enter the full path to this fasta file*

The algorithm has a built-in error correction step, but it does assume the input
data to be cleaned (i.e., adapter sequences should be removed and low quality
bases should be trimmed). This can be done using [cutadapt](https://pypi.python.org/pypi/cutadapt).


### Assembly parameters

* `--min_overlap_len`
By default this parameter is set to 60% of the average length of the sequencing
reads used as input for SAVAGE. The user can manually change this threshold
using the `-m` or `--min_overlap_len` parameter.
Increasing the minimal overlap length speeds up the algorithm and leads to a
lower mismatch rate in the final contigs. It is recommended to set the minimal
overlap length to be larger than the (expected) largest repetitive element in
the target genomes.
However, it also results in a lower fraction of the target genomes being
reconstructed. We advise you to consider this trade-off when setting this
parameter. We have been working with a minimal overlap length of 100-300bp,
depending on the read length.
* `--split`
In case of (ultra-)deep sequencing data, exceeding a coverage of
1000x, we advise to split the data into patches of coverage between 500x and
1000x and run SAVAGE Stage a on each patch individually. After specifying the
number of patches, SAVAGE takes care of the splitting and recombining. Choose
the number of patches using such that 500 < read_coverage/patch_num < 1000.
* `--merge_contigs`
By default this is set to 0, meaning that in stage c, the final assembly step,
we allow overlaps with a mismatch rate of 0% (i.e. exact overlaps). By
increasing this threshold, e.g. to 0.01, virus strains which differ by less
than 1% will be merged into one or more master strains, possibly leading to
longer contigs and a less fragmented assembly (higher N50).
* `--ignore_subreads`
By default this is set to False. When using the ignore_subreads flag, you choose
not to use subread information from previous stages in the current stage(s).
This will speed up the algorithm, but at the cost of less accurate abundance
estimates.
* `--sfo_mm`
This parameter is only relevant when running SAVAGE in de novo mode. It specifies
the -e parameter for running SFO. By default it is equal to 50, meaning that
the overlaps output by SFO allow up to 2% mismatches. This accounts for 1%
sequencing errors. Increasing this parameter will slow down the algorithm, while
decreasing will lead to a possibly incomplete overlap graph.
* `--overlap_len_stage_c`
For Stage c of the algorithm, the final assembly step, it is possible to specify
a different minimum overlap length using this option.
* `--contig_len_stage_c`
By default, only contigs of at least 100 bp in length are considered for stage c
assembly. The user can adjust this threshold by setting the
`contig_len_stage_c` parameter. From the final stage c output, it is usually a
good idea to consider only contigs of sufficient length, e.g. 500 bp.


### SAVAGE-ref

This mode constructs the overlap graph of the reads using the
pairwise overlaps induced from a read-to-reference alignment.
It takes as input a reference genome along with the read set,
and uses bwa mem to compute alignments after splitting the data.

SAVAGE has proven to work not only on high-quality, well-curated
reference genomes, but also on ad-hoc consensus genomes constructed
from the read set itself using a de novo assembler (such as VICUNA).
Therefore, even in ref-mode, SAVAGE can be used to compute very
accurate de novo viral quasispecies assemblies.

To run SAVAGE-ref, complete the following steps:

1. Download or assemble a reference genome `reference.fasta`.
2. Create a directory `my_directory` where you want to store the
    results and enter this directory: `cd my_directory`.
3. Now run SAVAGE with the option `--ref /path/to/reference.fasta`:
```python
savage --ref /path/to/reference.fasta --split patch_num --min_overlap_len M --s singles.fastq --p1 paired1.fastq --p2 paired2.fastq
```

### SAVAGE-de-novo

This mode does not rely on any reference sequence; instead, it
computes all approximate suffix-prefix overlaps among the reads
using the included package [SFO](https://www.cs.helsinki.fi/group/suds/sfo/)
from Valimaki et al.

To run SAVAGE-de-novo, complete the following steps:

1. Create a directory `my_directory` where you want to store the
    results and enter this directory: `cd my_directory`.
2. Now run SAVAGE:
```python
savage --split patch_num --min_overlap_len M --s singles.fastq --p1 paired1.fastq --p2 paired2.fastq
```


### Algorithm stages

The algorithm proceeds in three stages:
* Stage a has the original reads as input and contigs
* Stage b has these contigs as input and maximally extended contigs
  as output
* Stage c (OPTIONAL) merges maximized contigs into master strain
  contigs.
By default, SAVAGE runs all three stages, but the output of Stage b
(`contigs_stage_b.fasta`) is considered as the viral quasispecies
assembly.

SAVAGE performs error correction on the reads in the first iteration
of Stage a, which has been optimized for a coverage of 500x-1000x.


### Frequency estimation

We provide a frequency estimation procedure which computes for
every contig an estimate of its relative frequency. This procedure
needs a fasta/fastq file of contigs and the corresponding subreads
file (`subreads.txt`) as output by SAVAGE. This method works best
when applied to Stage b or Stage c contigs, because those sequences
have been maximally extended.

The frequency estimation can be applied as follows:

```python
freq_est --contigs contigs.fastq --subreads subreads.txt --min_len 1000
```

The parameter `--min_len` (or `-m`) can be set to any non-negative integer;
only contigs of at least this length will be considered during
frequency estimation.

For improved frequency estimation, we recommend to use [Kallisto](https://pachterlab.github.io/kallisto/).
This requires Kallisto to be installed on your machine. Kallisto uses
the original paired-end sequencing reads *(i.e., before running PEAR)*
from two separate files (`forward.fastq`, `reverse.fastq`). Run the
following command to estimate relative frequencies in Kallisto mode:
```python
freq_est --kallisto --kallisto_path /path/to/kallisto/executable -f forward.fastq -r reverse.fastq --contigs contigs.fastq --min_len 1000
```


### Example

The directory `example/` contains a small working example, with both
single- and paired-end reads contained in `input_fas/` and an HIV-1 HXB2
reference strain (`hiv-ref.fasta`).

To run the example, cd into `example/` and

* for SAVAGE-ref execute:

```python
savage --ref /path/to/example/hiv-ref.fasta -s input_fas/singles.fastq -p1 input_fas/paired1.fastq -p2 input_fas/paired2.fastq -m 200 --split 1
```

* for SAVAGE-de-novo execute:

```python
savage -s input_fas/singles.fastq -p1 input_fas/paired1.fastq -p2 input_fas/paired2.fastq -m 200 --split 1
```

Then apply the frequency estimation procedure on Stage b contigs:

```python
../freq_est --contigs stage_b/singles.fastq --subreads stage_b/subreads.txt --min_len 1000
```

The desired output of this final procedure (also on Stage a/c contigs) is given in the files
`frequencies_stage_*.ref.txt` and `frequencies_stage_*.denovo.txt`.


### High coverage data sets

In case of (ultra-)deep sequencing data, exceeding a coverage of
1000x, we advise to split the data into patches of coverage between 500x and
1000x and run SAVAGE Stage a on each patch individually. The resulting contigs
can then be combined for a joint Stage b and subsequently also Stage c.

### Snakemake workflow
We provide a **Snakefile** to run experiments with SAVAGE more easily with [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home).
All SAVAGE parameters need to be specified in the config file `savage_config.yaml`.
Please edit the config file and copy it to the directory where you wish to store
all SAVAGE output, then run
```
snakemake -s /path/to/SAVAGE_ROOT/Snakefile
```
where /path/to/SAVAGE_ROOT/ specifies the path to the SAVAGE root directory.
When re-executing a workflow, Snakemake will only run the stages of the algorithm
for which the input has actually changed, thus avoiding unnecessary computations.


## General remarks

SAVAGE expects as input single- and/or paired-end reads. For the
paired-end reads, it assumes they are stored both on the same strand
(hence resulting in F-F alignments) as output by [PEAR](http://sco.h-its.org/exelixis/web/software/pear/), unless specified otherwise. If your reads are
stored in forward-reverse orientations, make sure to use the option `--revcomp` when calling SAVAGE.


## Citation

If you are using SAVAGE, please cite our paper: *De novo viral quasispecies assembly using overlap graphs*,  
**J. Baaijens, A. Zine El Aabidine, E. Rivals, and A. Schoenhuth**,  
[doi:10.1101/080341](https://doi.org/10.1101/080341)


## Contact   

In case of any questions or issues, please contact Jasmijn Baaijens:
<lastname> AT cwi DOT nl
