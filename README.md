# HaploConduct

[![license](https://img.shields.io/badge/license-GPL%20v3.0-blue.svg)](http://www.gnu.org/licenses/)


## Synopsis

HaploConduct is a package designed for reconstruction of individual haplotypes
from next generation sequencing data, in particular Illumina. Currently,
HaploConduct consists of two methods: SAVAGE and POLYTE.

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
genomes. The principles of SAVAGE have been adjusted to a low-coverage setting,
where the ploidy of the organism is known. POLYTE follows an iterative scheme
where in each iteration reads or contigs are joined, based on their interplay in
terms of an underlying haplotype-aware overlap graph. Along the iterations,
contigs grow while preserving their haplotype identity. POLYTE has been shown
to produce very accurate and complete assemblies, reaching high target genome
reconstruction values at extremely low error rates.


## Input

Both methods are designed for Illumina paired-end sequencing reads. SAVAGE
assumes typical viral sequencing data consisting of at least 10.000x coverage;
POLYTE requires only 10x coverage per haplotype and achieves optimal performance
at 15-20x coverage per haplotype.


## Usage

For SAVAGE, run `python savage.py` and for POLYTE run `python polyte.py`.
Detailed user instructions can be found in the respective [SAVAGE](https://github.com/HaploConduct/HaploConduct/tree/master/savage) and
[POLYTE](https://github.com/HaploConduct/HaploConduct/tree/master/polyte) manuals.


## Contact   

In case of any questions or issues, please contact Jasmijn Baaijens:
<lastname> AT cwi DOT nl
