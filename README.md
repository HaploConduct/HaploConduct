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

## Installation

Please download and install:

* [boost](http://www.boost.org/) C++ libraries
* [ncbi-blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

Then download and enter the SAVAGE repository and type `make`. 

## Manual

Currently, SAVAGE is designed for Illumina MiSeq sequencing reads 
that have been pre-processed by [PEAR](http://sco.h-its.org/exelixis/web/software/pear/) (for merging self-overlapping 
pairs). Thus, the input consists of both single-end and paired-end 
reads.

SAVAGE offers two modes: SAVAGE-ref and SAVAGE-de-novo. SAVAGE-ref is the
more efficient option, while SAVAGE-de-novo gives slightly better results.
Instructions for using these modes can be found below.

### SAVAGE-ref

This mode constructs the overlap graph of the reads 
using the pairwise overlaps induced from a read-to-reference 
alignment. Hence, it takes as input a reference genome and alignment
files in SAM format (https://samtools.github.io/hts-specs/SAMv1.pdf). 

SAVAGE has proven to work not only on high-quality, well-curated 
reference genomes, but also on ad-hoc consensus genomes constructed 
from the read set itself using a de novo assembler (such as VICUNA).
Therefore, even in ref-mode, SAVAGE can be used to compute very 
accurate de novo viral quasispecies assemblies.

To run SAVAGE-ref, complete the following steps:

1. Download or assemble a reference genome `reference.fasta`.
2. Create a directory `my_directory` where you want to store the 
    results, and add a folder `pear_reads/` containing
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
python savage.py --ref reference.fasta --singles singles.sam --paired paired.sam
```

### SAVAGE-de-novo

This mode does not rely on any reference sequence; instead, it 
computes all approximate suffix-prefix overlaps among the reads 
using the included package SFO from Valimaki et al.
(https://www.cs.helsinki.fi/group/suds/sfo/). 

SAVAGE-de-novo is very easy to run, since it doesn't use any prior
information:

1. Create a directory `my_directory` where you want to store the 
    results, and add a folder `pear_reads/` containing
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
python savage.py
```

## Frequency estimation

We provide a frequency estimation procedure which computes for 
every contig an estimate of its relative frequency. This procedure 
needs a fasta/fastq file of contigs and the corresponding subreads 
file (`subreads.txt`) as output by SAVAGE. This method works best 
when applied to Stage b or Stage c contigs, because those sequences 
have been maximally extended. 

The frequency estimation can be applied as follows:

```python
python freq_est.py --fas contigs.fastq --subreads subreads.txt --min_len 5000
```

The parameter `--min_len` can be set to any non-negative integer; 
only contigs of at least this length will be considered during 
frequency estimation.


## Example

The directory `example/` contains a small working example, with both 
single- and paired-end reads contained in `pear_reads/` and alignments
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
python ../freq_est.py --fas stage_b/singles.fastq --subreads stage_b/subreads.txt --min_len 0
```
    
The desired output of this final procedure is given in the files 
`frequencies.ref.txt` and `frequencies.denovo.txt`.


## General remarks

SAVAGE expects as input single- and/or paired-end reads. For the 
paired-end reads, it assumes they are stored both on the same strand 
(hence resulting in F-F alignments) as output by [PEAR](http://sco.h-its.org/exelixis/web/software/pear/).

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

### High coverage data sets

In case of (ultra-)deep sequencing data, exceeding a coverage of 
1000x, it is currently required to split the data into patches of 
500x. A pipeline for processing such data sets will appear soon; 
until then, proceed as follows:

Split the data set into patches of 500x-1000x coverage and, on each of 
these subsets, run only SAVAGE Stage a:

```python
python savage.py --no_stage_b --no_stage_c
```
                        
Then concatenate the resulting contig files (`stage_a/singles.fastq`) 
into a file `combined_contigs.fastq`. Now enter `my_directory` and 
run SAVAGE Stages b and c:

```python
python savage.py --no_stage_a --contigs combined_contigs.fastq
```
                        
To apply frequency estimation on the resulting contigs (Stage b or c)
it is necessary to provide the subreads files (`subreads.txt`) for each
of the Stage a results. Use the option `--split-subreads` as follows:

```python
python ../freq_est.py --fas stage_b/singles.fastq --subreads stage_b/subreads.txt --min_len 0 --split-subreads /path/to/subreads1.txt,...,/path/to/subreadsk.txt
```
                          
The subreads files must be given as a list of paths to subreads files, 
separated by a ',' (no spaces!). It is very important that the subreads 
files are provided in the same order in which the corresponding fastq 
files were concatenated after Stage a.


## Contact   

In case of any questions or issues, please contact Jasmijn Baaijens: 
<lastname> AT cwi DOT nl
