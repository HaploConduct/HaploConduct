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

POLYTE uses a core algorithm implemented in C++11, managed by several workflow scripts implemented in Python2 (using scipy libraries).
The C++ part requires several boost libraries (boost::timer,
boost::system, and boost::program_options) and it needs a compiler
that supports [OpenMP](http://openmp.org/wp/) (such as `g++`).

For maximal clique enumeration, HaploConduct depends on the [quick-cliques package](https://github.com/darrenstrash/quick-cliques) which is already included.

For suffix-prefix overlap computations, HaploConduct uses the [rust-overlaps package](https://github.com/jbaaijens/rust-overlaps) for
computing suffix-prefix overlaps in de novo mode.

For reference-guided assembly, SAVAGE depends on [samtools](http://samtools.sourceforge.net) and the [bwa mem](http://bio-bwa.sourceforge.net/) aligner.

To summarize, please **download and install the following dependencies**:

* [rust-overlaps](https://github.com/jbaaijens/rust-overlaps)
* [bwa mem](http://bio-bwa.sourceforge.net/)
* [samtools](http://samtools.sourceforge.net)
* python2.7 + scipy

*Each of these tools can also be installed using [Bioconda](https://bioconda.github.io/),
a distribution of bioinformatics software realized as a channel for the
versatile Conda package manager. This comes down to one simple command, creating a conda environment that has all required dependencies:
`conda install --name haploconduct-deps python=2.7 scipy bwa samtools rust-overlaps`
Then activate the environment with `source activate haploconduct-deps` and you're ready to go!*

Once all dependencies are installed, download the [latest release](https://github.com/HaploConduct/HaploConduct/releases) of the HaploConduct package, enter the repository and type `make`.


## Manual

Start POLYTE by calling `haploconduct` in `polyte` mode:
```
haploconduct polyte
```

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

In addition, POLYTE requires the user to provide the following information about
the input data:

* `--hap_cov` **(required)**
The average per-haplotype coverage (i.e. sequencing depth). POLYTE needs at least 10x per-haplotype coverage to perform assembly, best results are obtained with 15-20x per-haplotype coverage.
* `--insert_size` **(required)**
The average insert size for the paired-end Illumina reads (i.e. 2 x read length + internal segment size).
* `--stddev` **(required)**
The standard deviation corresponding to the insert size distribution.

The algorithm has a built-in error correction step, but it does assume the input
data to be cleaned (i.e., adapter sequences should be removed and low quality
bases should be trimmed). For example, this can be done using [cutadapt](https://pypi.python.org/pypi/cutadapt).


### Parameters

POLYTE allows the user to tune the overlap length parameters manually. The initial error correction step requires a higher overlap length `-m_EC` than subsequent steps. By default `-m_EC` is set to 60% of the read length: for 250bp paired-end reads, this amounts to a minimal overlap length of 150bp for error correction. In subsequent steps, the minimal overlap length is set using `-m` (default=50bp).

* `-m_EC`, `--min_overlap_len_EC`
Minimal overlap length for error correction, default=0.6*readlength.
* `-m`, `--min_overlap_len`
Minimal overlap length after error correction, default=50.

Computing pairwise suffix-prefix overlaps is a time-consuming process, but it is highly parallelizable. Therefore, the user can speedup the algorithm by specifying the number of threads allowed (default=1).

* `-t`, `--num_threads`
Maximum number of threads allowed for overlap computations, default=1.


### Diploid mode (OPTIONAL)

For assembling diploid genomes, POLYTE has an optional diploid mode which can be activated by adding the `--diploid` flag when calling POLYTE. In diploid mode, the standard assembly algorithm is followed by an aditional merging step, exploiting the fact that there can only be two possible haplotypes. There are two additional parameters that the user can (optionally) specify when using diploid mode: the minimal overlap length and the minimal contig length.

* `--diploid`
Activate diploid mode.
* `--diploid_contig_len`
Minimal contig length required for contigs in diploid step.
* `--diploid_overlap_len`
Minimal overlap length used in diploid assembly step; by default equal to assembly `-m`.

*Note that in diploid mode, branches can be resolved based on very little information; this is more risky than default branch reduction and can possibly lead to (minor) misassemblies.*


### Output

The output contigs are stored in `contigs.fasta` and, when using diploid mode (optional), in `contigs_diploid.fasta`.


### Processing genomic regions larger than 100kb

Large regions, such as whole human chromosomes, require an approach that divides read set into bins, such that each bin can be processed individually. Currently, this binning is done based on alignments to a reference genome.

To process genomic regions larger than 100kb, call haploconduct in `polyte-split` mode instead of `polyte`. This pipeline takes the same input and parameters, along with three additional parameters:
* `--ref` **(required)**
Reference genome in fasta format, used for binning the reads')
* `--split_size`
Size of regions into which the reads are binned, default=10000.
* `--split_overlap`
Size of the overlap between regions over which the reads are divided, default=1000.


### Example

The directory `example/` contains a small working example. To give POLYTE a
test-run, cd into `example/` and execute:

```
../../haploconduct polyte -p1 input/forward.fastq -p2 input/reverse.fastq --diploid --hap_cov 14 --insert_size 486.6 --stddev 146.7
```


### Citation

If you are using POLYTE, please cite our paper: *Overlap graph-based generation of haplotigs for diploids and polyploids*,  
**J.A. Baaijens and A. Schoenhuth**,  
Bioinformatics 2019; 35(21): 4281--4289,
[doi:10.1101/378356](https://doi.org/10.1101/378356).


## Contact   

Please report any questions or issues [here](https://github.com/HaploConduct/HaploConduct/issues).
