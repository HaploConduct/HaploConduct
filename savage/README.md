# SAVAGE: Strain Aware VirAl GEnome assembly

[![license](https://img.shields.io/badge/license-GPL%20v3.0-blue.svg)](http://www.gnu.org/licenses/)

Current version: 0.4.0

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

For a more detailed description, please see our [Genome Research  paper](http://genome.cshlp.org/content/early/2017/04/10/gr.215038.116) or the [preprint](http://biorxiv.org/content/early/2017/01/21/080341) on BioRxiv.

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

## Installation and dependencies

The easiest and recommended way to install SAVAGE is through Conda, using the [Bioconda](https://bioconda.github.io/) channel. Alternatively, SAVAGE can also be installed from source code. The installation procedure for either option is described below.

*Please note that SAVAGE is built for linux-based systems only.*

### Installation via Conda
This is the easiest and recommended way to install SAVAGE. The Conda package manager allows you to install packages without needing any root privileges. First make sure to setup Bioconda by following the steps below:

- Install the **Miniconda Python2** distribution, see [here](https://conda.io/miniconda.html).
- Make sure that the miniconda2 bin directory is added to your PATH.
- Setup the bioconda channel as well as the other channels bioconda depends on; it is important to add them in this order:
```
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

More information about installing Bioconda can be found [here](https://bioconda.github.io/).

Once Bioconda is properly setup, SAVAGE and all its dependencies can be installed through **one simple command**:
```
conda install savage
```

### Installation from source

If installation via Conda is not an option for you, try installation from source.
The SAVAGE program is a combination of C++ code (which has to be
compiled) and Python scripts (version 2.6 or later, but not Python
3.x). The C++ part requires several boost libraries (boost::timer,
boost::system, and boost::program_options) and it needs a compiler
that supports [OpenMP](http://openmp.org/wp/) (such as `g++`).

For maximal clique enumeration, SAVAGE depends on the [quick-cliques package](https://github.com/darrenstrash/quick-cliques) which is already included.

As of version 0.4.0, SAVAGE uses the [rust-overlaps package](https://github.com/jbaaijens/rust-overlaps) for
computing suffix-prefix overlaps in de novo mode.

For reference-guided assembly, SAVAGE depends on the [bwa mem](http://bio-bwa.sourceforge.net/) aligner.

Stages b and c of the algorithm may require `blastn` and `makeblastdb`
from the ncbi C++ Toolkit.

We provide a simple algorithm for relative frequency estimation of the assembled
contigs. For improved frequency estimation and frequency-based filtering during assembly, we recommend using [Kallisto](https://pachterlab.github.io/kallisto/).
This is optional, for more information please see the frequency estimation and filtering sections below.

To summarize, please download and install:

* [rust-overlaps](https://github.com/jbaaijens/rust-overlaps)
* [ncbi-blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [bwa mem](http://bio-bwa.sourceforge.net/)
* [Kallisto](https://pachterlab.github.io/kallisto/)

*Each of these tools can also be installed using [Bioconda](https://bioconda.github.io/),
a distribution of bioinformatics software realized as a channel for the
versatile Conda package manager.*

Once all dependencies are installed, download the [latest release](https://github.com/HaploConduct/HaploConduct/releases) of the HaploConduct package, enter the repository and type `make`.


## Manual

Currently, SAVAGE is designed for Illumina MiSeq sequencing reads
that have been pre-processed by [PEAR](http://sco.h-its.org/exelixis/web/software/pear/) (for merging self-overlapping
pairs). Thus, the input consists of (a combination of) single-end and/or paired-end reads.
**Important:** SAVAGE expects paired-end reads stored in ***two fastq files with
both sequences in forward orientation*** (instead of the more standard
forward-reverse) because this how PEAR outputs non-merged reads. Please find a
more detailed description of the SAVAGE input files below.

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

The two most important parameters when running SAVAGE are `--split` and `--min_overlap_len`, which should be chosen with great care. The first is always required when running SAVAGE, while the second will be estimated from the input data if not specified otherwise.

Before running SAVAGE, please read through the basic parameter settings described below.

* `--split` **(required)**
In case of (ultra-)deep sequencing data, exceeding a coverage of
1000x, we advise to split the data into patches of coverage between 500x and
1000x and run SAVAGE Stage a on each patch individually. After specifying the
number of patches, SAVAGE takes care of the splitting and recombining. Choose
the number of patches using such that 500 < read_coverage/patch_num < 1000.
* `-m`, `--min_overlap_len`
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
* `-t`, `--num_threads`
By default the number of threads is set to 1. Computationally-heavy parts of SAVAGE are multi-threaded and will profit from an increased value of `-t`. Especially in de novo mode it is highly recommended to use a larger number of threads.
* `--merge_contigs`
By default this is set to 0, meaning that in stage c, the final assembly step,
we allow overlaps with a mismatch rate of 0% (i.e. exact overlaps). By
increasing this threshold, e.g. to 0.01, virus strains which differ by less
than 1% will be merged into one or more master strains, possibly leading to
longer contigs and a less fragmented assembly (higher N50).
* `--sfo_mm`
This parameter is only relevant when running SAVAGE in de novo mode. It specifies
the error rate allowed when computing approximate suffix-prefix overlaps. By default it is equal to 50, meaning that up to 2% mismatches is allowed in the overlaps. This accounts for 1%
sequencing errors. Increasing this parameter will slow down the algorithm, while
decreasing will lead to a possibly incomplete overlap graph.
* `--overlap_len_stage_c`
For Stage c of the algorithm, the final assembly step, it is possible to specify
a different minimum overlap length using this option. By default this parameter is set to 100, but depending on the data it can pay off to decrease this parameter further.
* `--contig_len_stage_c`
By default, only contigs of at least 100 bp in length are considered for stage c
assembly. The user can adjust this threshold by setting the
`contig_len_stage_c` parameter. From the final stage c output, it is usually a
good idea to consider only contigs of sufficient length, e.g. 500 bp.
* `--ignore_subreads`
By default this is set to False. When using the ignore_subreads flag, you choose
not to use subread information from previous stages in the current stage(s).
This will speed up the algorithm, but at the cost of less accurate abundance
estimates.
* `--max_tip_len`
Maximum extension length for a sequence to be called a tip in the overlap graph. By default this parameter is set to the average length of the input sequences. If you want to disable tip removal, set `--max_tip_len=0`. In general this will lead to a more fragmented assembly and is therefore not recommended.
* `--no_filtering`
By default, SAVAGE contigs are filtered after stages b and c based on Kallisto frequency estimates: all zero-abundance contigs are removed from the assembly. To disable this filtering procedure, add the `--no_filtering` flag.

### SAVAGE-de-novo

This mode does not rely on any reference sequence; instead, it
computes all approximate suffix-prefix overlaps among the reads
using the [rust-overlaps](https://github.com/sirkibsirkib/rust-overlaps) package. This is computationally more challenging than reference-guided
overlaps (see SAVAGE-ref) but avoids any bias towards a reference genome.

To run SAVAGE-de-novo, complete the following steps:

1. Create a directory `my_directory` where you want to store the
    results and enter this directory: `cd my_directory`.
2. Now run SAVAGE:
```python
savage --split patch_num --min_overlap_len M --s singles.fastq --p1 paired1.fastq --p2 paired2.fastq
```

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


### Algorithm stages

The algorithm proceeds in three stages:
* Stage a has the original reads as input and contigs. In this stage, SAVAGE performs error correction on the input reads. This has been optimized for a coverage of 500x-1000x, hence **it is very important to set the** `--split` **parameter correctly**.
* Stage b has these contigs as input and maximally extended contigs
  as output.
* Stage c attempts to merge maximized contigs into haplotypes. The parameter
`--merge_contigs` allows you to specify the mismatch rate allowed within overlaps in this final stage.
By increasing this threshold, e.g. to 0.01, virus strains which differ by less
than 1% will be merged into one or more *master strains*, possibly leading to
longer contigs and a less fragmented assembly (higher N50).

By default, SAVAGE runs all three stages and the output of Stage c
(`contigs_stage_c.fasta`) is considered as the viral quasispecies
assembly.


### Counting virus strains

When a sample contains strains with conserved regions larger than the read length, SAVAGE cannot connect variants across this region. In general, SAVAGE is very conservative: if it is not sure how to extend contigs, then it doesn't.
As a consequence, SAVAGE may not be able to produce full-length haplotypes. For this reason, we provide an estimated (lower bound) on the number of strains in the sample.
To apply this procedure, add the `--count_strains` flag to your SAVAGE command.

**Note**: to compute the estimated strain count a **reference genome** is required.
The reference genome should be given through the parameter `--ref`.
If you want to run SAVAGE in *de novo* mode it is required to do this separately. First run SAVAGE without specifying a reference genome and without the `--count_strains` flag. When the assembly is completed, run SAVAGE only for strain counting: on the command line, add flags `--no_assembly` and `--count_strains` and specify your reference genome with `--ref /path/to/your/reference.fasta`.


### Frequency estimation

We provide a frequency estimation procedure which computes for
every contig an estimate of its relative frequency. This procedure
needs a fasta/fastq file of contigs and the corresponding subreads
file (`subreads.txt`) as output by SAVAGE. This method works best
when applied to Stage b or Stage c contigs, because those sequences
have been maximally extended.

The frequency estimation can be applied as follows:

```python
python freq_est.py --contigs contigs.fastq --subreads subreads.txt --min_len 1000
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
freq_est --kallisto -f forward.fastq -r reverse.fastq --contigs contigs.fastq --min_len 1000
```

When using the Bioconda installation, please download `freq_est.py` manually from this repository.


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


### High coverage data sets

In case of (ultra-)deep sequencing data, exceeding a coverage of
1000x, we advise to split the data into patches of coverage between 500x and
1000x. This is easily achieved by setting the `--split` parameter in your SAVAGE command. This will run SAVAGE Stage a on each patch individually, followed by a joint Stage b on the full set of contigs obtained from all patches, and subsequently also Stage c.

Make sure to choose the number of patches in which to split your data such that `500 < read_coverage/patch_num < 1000`. Then set `--split=patch_num` and SAVAGE will take care of the rest.

### Snakemake workflow
The earliest versions of SAVAGE provided a **Snakefile** to run experiments with [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home). *However, this
approach is no longer maintained*; the complete SAVAGE workflow is now taken care of by `savage.py`. It is therefore recommended to run SAVAGE directly instead of using the Snakefile.

If you do wish to use the Snakefile, make sure that all SAVAGE parameters are specified in the config file `savage_config.yaml`.
Please edit the config file and copy it to the directory where you wish to store
all SAVAGE output, then run
```
snakemake -s /path/to/SAVAGE_ROOT/Snakefile
```
where /path/to/SAVAGE_ROOT/ specifies the path to the SAVAGE root directory.
Note that the most recent options of SAVAGE are not supported by the Snakefile.


## General remarks

SAVAGE expects as input single- and/or paired-end reads. For the
paired-end reads, it assumes they are stored both on the same strand
(hence resulting in F-F alignments) as output by [PEAR](http://sco.h-its.org/exelixis/web/software/pear/), unless specified otherwise. If your reads are
stored in forward-reverse orientations, make sure to use the option `--revcomp` when calling SAVAGE.


## Citation

If you are using SAVAGE, please cite our paper: *De novo viral quasispecies assembly using overlap graphs*,  
**J. Baaijens, A. Zine El Aabidine, E. Rivals, and A. Schoenhuth**,  
Genome Res. 2017. 27: 835-848,
[doi:10.1101/gr.215038.116](https://doi.org/10.1101/gr.215038.116).


## Contact   

In case of any questions or issues, please contact Jasmijn Baaijens:
baaijens AT cwi DOT nl
