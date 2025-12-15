# expos

Spatial information extraction for variant sites from alignment data,
plus comparison to background via monte-carlo simulation.
Useful for inspecting and flagging false postive variants.

Core functionality all present but niceties and guard rails are not.

More extensive documentation TODO!

## Installation

```bash
 # from the cloned repo
 mkdir build
 cd build
 cmake ..
 cmake --build .
 ./expos --help
```

## Usage

```
EXtract POSitional data and statistics from alignment at VCF variant sites. Annotated VCF to stdout.

Usage:
  expos [OPTION...] <VCF> <ALN>

  -h, --help          Print usage
  -i, --include arg   Only operate on VCF records with this value present
                      in FILTER. e.g. -i PASS. May be passed multiple
                      times.
  -e, --exclude arg   Only operate on VCF records without this value
                      present in FILTER. May be passed multiple times.
  -t, --tsv arg       Write a tsv of extended statistics to file specified.
  -n, --normal arg    Alignment for use as additional background data for
                      simulation
  -r, --ref arg       Alignment Reference Fasta for optionally adding
                      template kolmogorov complexity to statistics.
  -u, --uncompressed  output uncompressed VCF
```

## Annotations made

These are the header lines from an output VCF describing the INFO fields added. The [] notation is used to indicate which element of the array holds the data in question where the INFO field added is an array.

```
##INFO=<ID=KC,Number=1,Type=Integer,Description="Kolmogorov Complexity of region spanned by supporting templates, scaled by x100">
##INFO=<ID=TM1NN,Number=3,Type=Float,Description="[0]Median nearest neighbour distance of template endpoints from read pairs supporting variant; [1]log2 ratio effect size and [2]P-value against background, from monte-carlo simulation">
##INFO=<ID=QM1NN,Number=3,Type=Float,Description="[0]Median nearest neighbour distance of variant query position; [1]log2 ratio effect size and [2]P-value against background, from monte-carlo simulation">
##INFO=<ID=MLAS,Number=3,Type=Float,Description="[0]Median read-Length-normalised Alignment Score (AS) of reads supporting variant;[1]log2 ratio effect size and [2]P-value against background, from monte-carlo simulation">
```

## Example

```bash
 # example pipeline - Add some soft flags in the FILTER column
 # (or alternately, subset entirely with bcftools view instead of filter)

# line by line:
# 1: pipe VCF producing program to expos stdin
# 2: read from stdin, calculate statistics, output uncompressed VCF to stdout
# 3: statisically-backed flagging on distribution/clustering stats, comparing variant-supporting data to background
# 4: statistically-backed flagging on alignment score, comparing variant-supporting data to background
# 5: conservative heuristic flagging poor alignment score in regions of low complexity, and write to disk
<some vcf producing command> |
./path/to/expos -u --ref ref.fa - my.bam |
bcftools filter --mode + -s SPATIAL -e'(INFO/QM1NN[1] <= -1.0 & INFO/QM1NN[2] <= 0.05) | (INFO/TM1NN[1] <= -1.0 & INFO/TM1NN[2] <= 0.05)' |
bcftools filter --mode + -s LOW_AS -e'(INFO/MLAS[0] < 0.93 & INFO/MLAS[1] < 0 & INFO/MLAS[2] <= 0.05)' |
bcftools filter --mode + -s LOW_CMPLX -e'(INFO/MLAS[0] < 0.93 & INFO/KC < 150)' > my.flagged.vcf
```
