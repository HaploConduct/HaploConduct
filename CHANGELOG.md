# Change log

## [0.2.1] - 2017-01-25
### Added
- Option --revcomp to enable processing of forward-reverse paired-end reads instead
of only forward-forward as output by PEAR

### Changed
- Print error messages to stderr instead of stdout.

## [0.2.0] - 2017-01-17
### Added
- Frequency estimation using Kallisto.
- Command line argument --overlap_len_stage_c to adjust the minimum overlap length required in Stage c
- Command line argument --contig_len_stage_c to set a minimum length on the input sequences (contigs) in Stage c

### Changed
- Fixed reference fasta file reading to allow for fasta wrapping.
- Updated wrapper script to make processing of large data sets possible without the Snakefile. Executing SAVAGE using Snakemake is still very practical but no longer compulsory.
- Updated Snakefile to work with updated wrapper script

## [0.1.0] - 2016-12-21
Initial pre-release.
### Comments
- Large data sets need to be processed using Snakemake together with the Snakefile provided.
