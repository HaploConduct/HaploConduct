# POLYTE: POLYploid genome fitTEr

## Synopsis

POLYTE is a method for reconstructing haplotigs for diploid and polyploid
genome, specifically designed for a low-coverage setting where the ploidy of
the organism is known. POLYTE follows an iterative scheme where in each
iteration reads or contigs are joined, based on their interplay in terms of
an underlying haplotype-aware overlap graph. Along the iterations,
contigs grow while preserving their haplotype identity. POLYTE has been shown
to produce very accurate and complete assemblies, reaching high target genome
reconstruction values at extremely low error rates.

## Installation and dependencies

POLYTE uses a core algorithm implemented in C++11, managed by several workflow scripts implemented in Python2.
The C++ part requires several boost libraries (boost::timer,
boost::system, and boost::program_options) and it needs a compiler
that supports [OpenMP](http://openmp.org/wp/) (such as `g++`).

For maximal clique enumeration, HaploConduct depends on the [quick-cliques package](https://github.com/darrenstrash/quick-cliques) which is already included.

For suffix-prefix overlap computations, HaploConduct uses the [rust-overlaps package](https://github.com/jbaaijens/rust-overlaps) for
computing suffix-prefix overlaps in de novo mode.

For reference-guided assembly, SAVAGE depends on the [bwa mem](http://bio-bwa.sourceforge.net/) aligner.

To summarize, please **download and install the following dependencies**:

* [rust-overlaps](https://github.com/jbaaijens/rust-overlaps)
* [bwa mem](http://bio-bwa.sourceforge.net/)

*Each of these tools can also be installed using [Bioconda](https://bioconda.github.io/),
a distribution of bioinformatics software realized as a channel for the
versatile Conda package manager.*

Once all dependencies are installed, download the [latest release](https://github.com/HaploConduct/HaploConduct/releases) of the HaploConduct package, enter the repository and type `make`.


## Manual

### Input

POLYTE is designed for Illumina paired-end sequencing reads, which should be
stored in separate fastq files.
These files can be specified using the following program options:

* `-p1` or `--input_p1`
Input fastq file containing /1 paired-end reads.
* `-p2` or `--input_p2`
Input fastq file containing /2 paired-end reads.
* `-s` or `--input_s`
Input fastq file containing any unpaired reads.

The algorithm has a built-in error correction step, but it does assume the input
data to be cleaned (i.e., adapter sequences should be removed and low quality
bases should be trimmed). For example, this can be done using [cutadapt](https://pypi.python.org/pypi/cutadapt).

### Parameters

* `--hap_cov` **(required)**
* `--insert_size` **(required)**
* `--stddev` **(required)**
* `-t`, `--num_threads`
* `-m`, `--min_overlap_len`
* `-m_EC`, `--min_overlap_len_EC`

### Diploid mode

* `--diploid`
* `--diploid_contig_len`
* `--diploid_overlap_len`

### Output

The output contigs are stored in `contigs.fasta` and, when using diploid mode (optional), in `contigs_diploid.fasta`.

### Processing genomic regions larger than 100kb

### Example


## Contact   

In case of any questions or issues, please contact Jasmijn Baaijens:
<lastname> AT cwi DOT nl
