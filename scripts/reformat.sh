#!/bin/bash
#reformat in=<infile> out=<outfile>

function usage(){
echo "
Written by Brian Bushnell
Last modified October 23, 2015

Description:  Reformats reads to change ASCII quality encoding, interleaving, file format, or compression format.
Optionally performs additional functions such as quality trimming, subsetting, and subsampling.
Supports sam, fastq, fasta, fasta+qual, scarf, gzip, zip.

Usage:  reformat.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2>

in2 and out2 are for paired reads and are optional.
If input is paired and there is only one output file, it will be written interleaved.


Other parameters and their defaults:

ow=f                    (overwrite) Overwrites files that already exist.
app=f                   (append) Append to files that already exist.
zl=4                    (ziplevel) Set compression level, 1 (low) to 9 (max).
int=f                   (interleaved) Determines whether INPUT file is considered interleaved.
fastawrap=70            Length of lines in fasta output.
fastareadlen=0          Set to a non-zero number to break fasta files into reads of at most this length.
fastaminlen=1           Ignore fasta reads shorter than this.
qin=auto                ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto               ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
qfake=30                Quality value used for fasta to fastq reformatting.
qfin=<.qual file>       Read qualities from this qual file, for the reads coming from 'in=<fasta file>'
qfin2=<.qual file>      Read qualities from this qual file, for the reads coming from 'in2=<fasta file>'
qfout=<.qual file>      Write qualities from this qual file, for the reads going to 'out=<fasta file>'
qfout2=<.qual file>     Write qualities from this qual file, for the reads coming from 'out2=<fasta file>'
outsingle=<file>        (outs) If a read is longer than minlength and its mate is shorter, the longer one goes here.

Processing Parameters:

verifypaired=f          (vpair) When true, checks reads to see if the names look paired.  Prints an error message if not.
verifyinterleaved=f     (vint) sets 'vpair' to true and 'interleaved' to true.
allowidenticalnames=f   (ain) When verifying pair names, allows identical names, instead of requiring /1 and /2 or 1: and 2:
tossbrokenreads=f       (tbr) Discard reads that have different numbers of bases and qualities.  By default this will be detected and cause a crash.
ignorebadquality=f      (ibq) Fix out-of-range quality values instead of crashing with a warning.
addslash=f              Append ' /1' and ' /2' to read names, if not already present.  Please include the flag 'int=t' if the reads are interleaved.
spaceslash=t            Put a space before the slash in addslash mode.
underscore=f            Change whitespace in read names to underscores.
rcomp=f                 (rc) Reverse-compliment reads.
rcompmate=f             (rcm) Reverse-compliment read 2 only.
changequality=t         (cq) N bases always get a quality of 0 and ACGT bases get a min quality of 2.
quantize=f              Quantize qualities to a subset of values like NextSeq.  Can also be used with comma-delimited list, like quantize=0,8,13,22,27,32,37
tuc=f                   (touppercase) Change lowercase letters in reads to uppercase.
uniquenames=f           Make duplicate names unique by appending _<number>.
remap=                  A set of pairs: remap=CTGN will transform C>T and G>N.
                        Use remap1 and remap2 to specify read 1 or 2.
iupacToN=f              (itn) Convert non-ACGTN symbols to N.
monitor=f               Kill this process if it crashes.  monitor=600,0.01 would kill after 600 seconds under 1% usage.
tossjunk=f              Discard reads with invalid characters as bases.
fixjunk=f               Convert invalid bases to N.
fixheaders=f            Convert nonstandard header characters to standard ASCII.
recalibrate=f           (recal) Recalibrate quality scores.  Must first generate matrices with CalcTrueQuality.
maxcalledquality=41     Quality scores capped at this upper bound.
mincalledquality=2      Quality scores of ACGT bases will be capped at lower bound.
trimreaddescription=f   (trd) Trim the names of reads after the first whitespace.

Histogram output parameters:

bhist=<file>            Base composition histogram by position.
qhist=<file>            Quality histogram by position.
qchist=<file>           Count of bases with each quality value.
aqhist=<file>           Histogram of average read quality.
bqhist=<file>           Quality histogram designed for box plots.
lhist=<file>            Read length histogram.
gchist=<file>           Read GC content histogram.
gcbins=100              Number gchist bins.  Set to 'auto' to use read length.
gcplot=f                Add a graphical representation to the gchist.

Histograms for sam files only (requires sam format 1.4 or higher):

ehist=<file>            Errors-per-read histogram.
qahist=<file>           Quality accuracy histogram of error rates versus quality score.
indelhist=<file>        Indel length histogram.
mhist=<file>            Histogram of match, sub, del, and ins rates by read location.
idhist=<file>           Histogram of read count versus percent identity.
idbins=100              Number idhist bins.  Set to 'auto' to use read length.

Sampling parameters:

reads=-1                Set to a positive number to only process this many INPUT reads (or pairs), then quit.
skipreads=-1            Skip (discard) this many INPUT reads before processing the rest.
samplerate=1            Randomly output only this fraction of reads; 1 means sampling is disabled.
sampleseed=-1           Set to a positive number to use that prng seed for sampling (allowing deterministic sampling).
samplereadstarget=0     (srt) Exact number of OUTPUT reads (or pairs) desired.
samplebasestarget=0     (sbt) Exact number of OUTPUT bases desired.
                        Important: srt/sbt flags should not be used with stdin, samplerate, qtrim, minlength, or minavgquality.

Trimming and filtering parameters:

qtrim=f                 Trim read ends to remove bases with quality below trimq.
                        Values: t (trim both ends), f (neither end), r (right end only), l (left end only), w (sliding window).
trimq=6                 Regions with average quality BELOW this will be trimmed.
minlength=0             (ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter.
mlf=0                   (mlf) Reads shorter than this fraction of original length after trimming will be discarded.
maxlength=0             If nonzero, reads longer than this after trimming will be discarded.
breaklength=0           If nonzero, reads longer than this will be broken into multiple reads of this length.  Does not work for paired reads.
requirebothbad=t        (rbb) Only discard pairs if both reads are shorter than minlen.
minavgquality=0         (maq) Reads with average quality (after trimming) below this will be discarded.
maqb=0                  If positive, calculate maq from this many initial bases.
chastityfilter=f        (cf) Reads with names  containing ' 1:Y:' or ' 2:Y:' will be discarded.
barcodefilter=f         Remove reads with unexpected barcodes if barcodes is set, or barcodes containing 'N' otherwise.  
                        A barcode must be the last part of the read header.
barcodes=               Comma-delimited list of barcodes or files of barcodes.
maxns=-1                If 0 or greater, reads with more Ns than this (after trimming) will be discarded.
minconsecutivebases=0   (mcb) Discard reads without at least this many consecutive called bases.
forcetrimleft=0         (ftl) If nonzero, trim left bases of the read to this position (exclusive, 0-based).
forcetrimright=0        (ftr) If nonzero, trim right bases of the read after this position (exclusive, 0-based).
forcetrimright2=0       (ftr2) If positive, trim this many bases on the right end.
forcetrimmod=5          (ftm) If positive, trim length to be equal to zero modulo this number.
mingc=0                 Discard reads with GC content below this.
maxgc=1                 Discard reads with GC content above this.

Sam and bam processing options:

mappedonly=f            Toss unmapped reads.
unmappedonly=f          Toss mapped reads.
primaryonly=f           Toss secondary alignments.  Set this to true for sam to fastq conversion.
requiredbits=0          (rbits) Toss sam lines with any of these flag bits unset.  Similar to samtools -f.
filterbits=0            (fbits) Toss sam lines with any of these flag bits set.  Similar to samtools -F.
stoptag=f               Set to true to write a tag indicating read stop location, prefixed by YS:i:
sam=                    Set to 'sam=1.3' to convert '=' and 'X' cigar symbols (from sam 1.4+ format) to 'M'.

Cardinality estimation:
cardinality=f           (loglog) Count unique kmers using the LogLog algorithm.
loglogk=31              Use this kmer length for counting.
loglogbuckets=1999      Use this many buckets for counting.

Shortcuts: 
The # symbol will be substituted for 1 and 2.  The % symbol in out will be substituted for input name minus extensions.
For example:
reformat.sh in=read#.fq out=%.fa
...is equivalent to:
reformat.sh in1=read1.fq in2=read2.fq out1=read1.fa out2=read2.fa

Java Parameters:
-Xmx                    This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                        -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

Supported input formats are fastq, fasta, fast+qual, scarf, and bread (BBMap's native format)
Supported output formats are fastq, fasta, fast+qual, bread, sam, and bam (bam only if samtools is installed)
Supported compression formats are gz, zip, and bz2
To read from stdin, set 'in=stdin'.  The format should be specified with an extension, like 'in=stdin.fq.gz'
To write to stdout, set 'out=stdout'.  The format should be specified with an extension, like 'out=stdout.fasta'

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

#Old sam options - these probably only work with samv1.4+ input.  TODO: test.
#build=<integer>         Assign a genome's build id.  You can index like this: bbmap.sh ref=<file> build=1
#sam=1.4                 Set to 1.4 to write Sam version 1.4 cigar strings, with = and X, or 1.3 to use M.
#md=f                    Set to true to write MD tags.
#xs=f                    Set to 'ss', 'fs', or 'us' to write XS tags for RNAseq using secondstrand, firststrand,
#                        or unstranded libraries.  Needed by Cufflinks.  JGI mainly uses 'firststrand'.
#idtag=t                 Set to true to write a tag indicating percent identity, prefixed by YI:f:

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx200m"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
}
calcXmx "$@"

function reformat() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module unload samtools
		module load oracle-jdk/1.7_64bit
		module load pigz
		module load samtools
	fi
	local CMD="java $EA $z -cp $CP jgi.ReformatReads $@"
	echo $CMD >&2
	eval $CMD
}

reformat "$@"
