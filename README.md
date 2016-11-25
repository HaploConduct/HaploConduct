# SAVAGE: Strain Aware VirAl GEnome assembly

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

Stages b and c of the algorithm also require `blastn` and `makeblastdb`
from the ncbi C++ Toolkit.

Optionally, the different preprocessing of the input data and the different
algorithm stages can all be run in a single command using [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) to execute the
Snakefile from the SAVAGE root directory.

## Installation

Please download and install:

* [ncbi-blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) (OPTIONAL)

*These tools can also be installed using [Bioconda](https://bioconda.github.io/),
a distribution of bioinformatics software realized as a channel for the
versatile Conda package manager*

Then download and enter the SAVAGE repository and type `make`.

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

### SAVAGE-ref

This mode constructs the overlap graph of the reads
using the pairwise overlaps induced from a read-to-reference
alignment. Hence, it takes as input a reference genome and alignment
files in [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf).

SAVAGE has proven to work not only on high-quality, well-curated
reference genomes, but also on ad-hoc consensus genomes constructed
from the read set itself using a de novo assembler (such as VICUNA).
Therefore, even in ref-mode, SAVAGE can be used to compute very
accurate de novo viral quasispecies assemblies.

To run SAVAGE-ref, complete the following steps:

1. Download or assemble a reference genome `reference.fasta`.
2. Create a directory `my_directory` where you want to store the
    results, and add a folder `input_fas/` containing
    ```
    singles.fastq,
    paired1.fastq,
    paired2.fastq.
    ```
    Make sure that the read identifiers are numerical, with IDs 0 to
    num\_singles-1 for the single-end reads and IDs num\_singles
    to num\_singles+num\_pairs-1 for the paired-end reads.
3. Align the single- and paired-end reads seperately to this reference,
    for example using [bwa-mem](http://bio-bwa.sourceforge.net/), obtaining
    `singles.sam` and `paired.sam`.
4. Now you're ready to run SAVAGE: enter `my_directory` and run
```python
python savage.py --min_overlap_len M --ref reference.fasta --singles singles.sam --paired paired.sam
```

### SAVAGE-de-novo

This mode does not rely on any reference sequence; instead, it
computes all approximate suffix-prefix overlaps among the reads
using the included package [SFO](https://www.cs.helsinki.fi/group/suds/sfo/)
from Valimaki et al.

SAVAGE-de-novo is very easy to run, because it doesn't require any prior
information:

1. Create a directory `my_directory` where you want to store the
    results, and add a folder `input_fas/` containing
    ```
    singles.fastq,
    paired1.fastq,
    paired2.fastq.
    ```
    Make sure that the read identifiers are numerical, with IDs 0 to
    num\_singles-1 for the single-end reads and IDs num\_singles
    to num\_singles+num\_pairs-1 for the paired-end reads.
2. Now you're ready to run SAVAGE: enter `my_directory` and run
```python
python savage.py --min_overlap_len M
```


### Algorithm stages

The algorithm proceeds in three stages:
* Stage a has the original reads as input and contigs
* Stage b has these contigs as input and maximally extended contigs
  as output
* Stage c (OPTIONAL) merges maximized contigs into master strain
  sequences.
By default, SAVAGE runs all three stages, but the output of Stage b
(`contigs_stage_b.fasta`) is considered as the viral quasispecies
assembly.

SAVAGE performs error correction on the reads in the first iteration
of Stage a, which has been optimized for a coverage of 500x-1000x.


### Assembly parameters

* `--min_overlap_len`
By default this is set to 150bp, which has proven to work well for Illumina
miseq (2x250bp) reads. Increasing this threshold speeds up the algorithm and
leads to a lower mismatch rate in the final contigs.
However, it also results in a lower fraction of the target genomes being
reconstructed. We advise you to consider this trade-off when setting this
parameter. We have been working with a minimal overlap length of 100-200bp.
* `--merge_contigs`
By default this is set to 0.01, meaning that in stage c (master strain assembly)
we allow overlaps with a mismatch rate up to 1%. By increasing this threshold,
you are likely to collapse several strains into one or more master strains, but
possibly leading to longer contigs and hence a higher N50.
* `--use_subreads`
By default this is set to True. When setting use_subreads to False, you choose
not to use subread information from previous stages in the current stage(s).
This will speed up the algorithm, but at the cost of less accurate abundance
estimates.
* `--sfo_mm`
This parameter is only relevant when running SAVAGE in de novo mode. It specifies
the -e parameter for running SFO. By default it is equal to 50, meaning that
the overlaps output by SFO allow up to 2% mismatches. This accounts for 1%
sequencing errors. Increasing this parameter will slow down the algorithm, while
decreasing will lead to a possibly incomplete overlap graph.


### Frequency estimation

We provide a frequency estimation procedure which computes for
every contig an estimate of its relative frequency. This procedure
needs a fasta/fastq file of contigs and the corresponding subreads
file (`subreads.txt`) as output by SAVAGE. This method works best
when applied to Stage b or Stage c contigs, because those sequences
have been maximally extended.

The frequency estimation can be applied as follows:

```python
python freq_est.py --fas contigs.fastq --subreads subreads.txt --min_len 1000
```

The parameter `--min_len` can be set to any non-negative integer;
only contigs of at least this length will be considered during
frequency estimation.


### Example

The directory `example/` contains a small working example, with both
single- and paired-end reads contained in `input_fas/` and alignments
`singles.sam`, `paired.sam` to an HIV-1 HXB2 reference strain
(`hiv-ref.fasta`).

To run the example, cd into `example/` and

* for SAVAGE-ref execute:

```python
python savage.py --ref hiv-ref.fasta --singles singles.sam --paired paired.sam
```

* for SAVAGE-de-novo execute:

```python
python savage.py
```

Then apply the frequency estimation procedure on Stage b contigs:

```python
python ../freq_est.py --fas stage_b/singles.fastq --subreads stage_b/subreads.txt --min_len 1000
```

The desired output of this final procedure (also on Stage a/c contigs) is given in the files
`frequencies_stage_*.ref.txt` and `frequencies_stage_*.denovo.txt`.


### High coverage data sets

In case of (ultra-)deep sequencing data, exceeding a coverage of
1000x, we advise to split the data into patches of coverage between 500x and
1000x and run SAVAGE Stage a on each patch individually. The resulting contigs
can then be combined for a joint Stage b and subsequently also Stage c. We
provide a **Snakefile** to run the entire procedure. As input it takes your
fastq-files: single and/or paired-end reads, where the paired-end reads are
provided in separate files in forward-forward orientation. The file paths need
to be specified in the config file `savage_config.yaml`, along with the other
SAVAGE parameters. Please edit the config file and copy it to the directory
where you wish to store all SAVAGE output, then run
```
snakemake -s /path/to/SAVAGE_DIR/Snakefile
```
where /path/to/SAVAGE_DIR specifies the path to the SAVAGE root directory.
This pipeline uses [bwa-mem](http://bio-bwa.sourceforge.net/) when choosing reference-guided mode, which can also be installed using [Bioconda](https://bioconda.github.io/).


## General remarks

SAVAGE expects as input single- and/or paired-end reads. For the
paired-end reads, it assumes they are stored both on the same strand
(hence resulting in F-F alignments) as output by [PEAR](http://sco.h-its.org/exelixis/web/software/pear/). If your reads are
stored in forward-reverse orientations, this means you first need to take
reverse complements of your /2 sequences.


## Contact   

In case of any questions or issues, please contact Jasmijn Baaijens:
<lastname> AT cwi DOT nl
