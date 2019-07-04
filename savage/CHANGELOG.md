# Change log

## [0.4.1] - 2019‑03‑15
### Added
- Command line argument -o/--outdir
### Changed
Minor bugfixes, improved error handling.

## [0.4.0] - 2017-07-17
### Added
- Frequency-based filtering using Kallisto
- Estimate strain count with --count_strains flag
- User options for specifying max_tip_len and min_clique_size
### Changed
- Replaced SFO tool by rust-overlaps for faster de novo overlaps, also allowing exact overlap computations in stages b and c instead of running blast
- Take cigar strings into account when computing overlaps from read-to-reference alignments
- Keep paired-end contigs as single-end after stage a (instead of removal)
- Fixed bug when -m is odd

## [0.3.0] - 2017-03-06
### Added
- Remove contigs which are fully included in another contig without any mismatches
- Remove all tip sequences from the read set and store them, together with all inclusions, in a separate fastq file `removed_tip_sequences.fastq` per stage
- More parameter options in the config file when using Snakemake
### Changed
- Increase edge threshold (overlap score) in stages b and c
- Higher threshold for placing 'N's in contigs: phred < 20
- Improved default settings and small bug fixes

## [0.2.2] - 2017-02-14
### Changed
- Fixed bug when using only paired-end reads OR single-end reads
- Allow for non-integer read identifiers
- Disregard contigs with more than 5% uncalled bases (N's)
- Write contig sequences corresponding to tips in the overlap graph to a separate fastq file

## [0.2.1] - 2017-01-25
### Added
- Option --revcomp to enable processing of forward-reverse paired-end reads instead of only forward-forward as output by PEAR

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
