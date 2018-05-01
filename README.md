# HaploConduct

[![license](https://img.shields.io/badge/license-GPL%20v3.0-blue.svg)](http://www.gnu.org/licenses/)


## Synopsis

HaploConduct is a package designed for reconstruction of individual haplotypes
from next generation sequencing data, in particular Illumina. Currently,
HaploConduct consists of two methods: [SAVAGE](https://github.com/HaploConduct/HaploConduct/tree/master/savage) and
[POLYTE](https://github.com/HaploConduct/HaploConduct/tree/master/polyte).

### SAVAGE: Strain Aware VirAl GEnome assembly

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


### POLYTE: POLYploid genome fitTEr

POLYTE is a method for reconstructing haplotigs for diploid and polyploid
genome, specifically designed for a low-coverage setting where the ploidy of
the organism is known. POLYTE follows an iterative scheme where in each
iteration reads or contigs are joined, based on their interplay in terms of
an underlying haplotype-aware overlap graph. Along the iterations,
contigs grow while preserving their haplotype identity. POLYTE has been shown
to produce very accurate and complete assemblies, reaching high target genome
reconstruction values at extremely low error rates.


## Installation and dependencies

*Note: the easiest and recommended way to install SAVAGE is through bioconda,
see also the [SAVAGE manual](https://github.com/HaploConduct/HaploConduct/tree/master/savage).
Installation through bioconda takes care of all dependencies and allows you to
skip the steps below. We hope to make POLYTE available on bioconda soon.*

SAVAGE and POLYTE make use of the same C++ codebase but have completely
different workflows (implemented in Python2).
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


## Input

Both methods are designed for Illumina paired-end sequencing reads. SAVAGE
assumes typical viral sequencing data consisting of at least 10.000x coverage;
POLYTE requires only 10x coverage per haplotype and achieves optimal performance
at 15-20x coverage per haplotype.


## Usage

*Note: the HaploConduct workflows are currently implemented in Python2, we plan
to move to Python3 in the very near future.*

For SAVAGE, run `python savage.py` and for POLYTE run `python polyte.py`.
Detailed user instructions can be found in the respective [SAVAGE](https://github.com/HaploConduct/HaploConduct/tree/master/savage) and
[POLYTE](https://github.com/HaploConduct/HaploConduct/tree/master/polyte) manuals.


## Contact   

In case of any questions or issues, please contact Jasmijn Baaijens:
baaijens AT cwi DOT nl
