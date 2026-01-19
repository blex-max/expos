# expos

Statistically-backed VCF flagging calculating effect sizes
and p-values for spatial properties of somatic mutations via
monte-carlo simulation. Compares mutant reads (supporting observations)
with the set of all reads (background), optionally including reads
from a matched normal.

Useful for inspecting and flagging false postive variants caused
by a variety of processes most commonly associated with a high
number of PCR cycles. Builds on Ellis et al. 2021, amongst others.

For both SNVs and indels. MNV handling logic is present,
but largely untested.

Core functionality present but niceties and guard rails are not.
Please report any bugs!

#### TODO
 - report correct position in error messages - chrom and position are both
   off by one! (this affects error messages only, and has no impact on the output)
 - More extensive documentation
 - Possibly remove p-value for MLAS/ASRD, as it can be misleading
 - Add clip fraction/CLPM
 - record command used in VCF header

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
EXtract POSitional data and statistics from alignment at
VCF variant sites, and encode them as INFO fields to VCF.
Requires the presence of .(b/cr)ai indexes of the same name
as the relevant alignment. Annotated VCF to stdout. See
README or output VCF header for descriptions of fields
added.

Usage:
  expos [OPTION...] <VCF/BCF (- for stdin)> <ALN.(b/cr)am>

  -h, --help              Print usage
  -i, --include arg       Only operate on VCF records with this value
                          present in FILTER. e.g. -i PASS. May be passed
                          multiple times.
  -e, --exclude arg       Only operate on VCF records without this value
                          present in FILTER. May be passed multiple times.
  -f, --flag-include arg  Only consider reads with these bits set in the
                          SAM flag. Applies to both target and background
                          alignment data. Default: 3
  -F, --flag-exclude arg  Do not consider reads with these bits set in the
                          SAM flag. Applies to both target and background
                          alignment data. Default: 3852
  -t, --tsv arg           Write a tsv of extended statistics to file
                          specified.
  -n, --normal arg        Alignment for use as additional background data
                          for simulation
  -r, --ref arg           Alignment Reference Fasta for optionally adding
                          template complexity to statistics.
  -u, --uncompressed      output uncompressed VCF
```
basic usage then looks like:
```bash
expos my.vcf my.bam
```


## Assessments made

Assessments of spatial clustering of mutant bases have been used to filter false-positive mutations which may otherwise be difficult to consistently detect (Ellis, Peter Cambpell, GATK ReadPosRankSum, others). The underlying hypothesis is that mutant reads should be drawn from the same spatial distribution as reference reads; if mutant reads differ significantly from reference reads, then the spatial process producing those reads deviates from the spatial process producing the reference reads. This may indicate that a non-biological process, or sequencing artefact, is responsible for the mutant reads since it would not be expected that mutant reads exhibit a unique preference for a particular region.

expos implements a nearest-neighbour algorithm (Cover, 1967) on two spatial properties of the set of mutant reads; the query position of the mutation on each read, and the endpoints of the inferred template from which each read was amplified. For each property, expos finds the distance to the single closest neighbour for each read, and reports the median of the set of these distances. Since the unit of these metrics is in sequence bases, they are readily interpretable as descriptive statistics – what is the average distance in bases to the closest neighboring observation?

The underlying hypothesis is that mutant reads should be drawn from the same spatial distribution as reference reads; if mutant reads differ significantly from reference reads, then the spatial process producing those reads deviates from the spatial process producing the reference reads. This may indicate that a non-biological process, or sequencing artefact, is responsible for the mutant reads since it would not be expected that mutant reads exhibit a unique preference for a particular region. These statistics are an advancement on the heuristic approach to this kind of artefact detection put forward by Ellis et al., and GATK's ReadPosRankSumTest

These are the header lines from an output VCF describing the INFO fields added. The [] notation is used to indicate which element of the array holds the data in question where the INFO field added is an array.

```
##INFO=<ID=QM1NN,Number=3,Type=Float,Description="[0]Median nearest neighbour distance of variant query position; [1]log2 ratio effect size and [2]two-sided P-value against background, from monte-carlo simulation">
##INFO=<ID=TM1NN,Number=3,Type=Float,Description="[0]Median nearest neighbour distance of template endpoints from read pairs supporting variant; [1]log2 ratio effect size and [2]two-sided P-value against background, from monte-carlo simulation">
##INFO=<ID=MLAS,Number=3,Type=Float,Description="[0]Median read-Length-normalised Alignment Score (AS) of reads supporting variant;[1]delta (supporting - background) effect size and [2]two-sided P-value against background, from monte-carlo simulation">
##INFO=<ID=RCMPLX,Number=1,Type=Integer,Description="Complexity (Lempel-Ziv estimated entropy rate) of region spanned by supporting templates, scaled by x100">
```

log2 effect sizes scale such no effect is 0, -1 means the statistic is 1/2 on supporting data compared background, -2 1/4, whereas an effect size of 1 means the statistic 2x on supporting data compared to background, 2 4x. Practically this means that effect sizes below 0 indicate tighter clustering of observations as compared to background.

The effect size for median length-normalised alignment score is simply reported as difference between the statisic as calculated on the supporting data and the mean of all simulated calculations. Therefore a lower alignment score on supporting reads as compared to background data results in a negative effect size. Note that the presence of a variant will by definition lower alignment score so significant effects do not necessarily indicate a spurious variant, particularly in the case of larger indels.

MLAS is equivalent to ASRD as may be familiar to some users.

## Example

```bash
 # example pipeline - Add some soft flags in the FILTER column
 # (or alternately, subset entirely with bcftools view instead of filter)

# line by line:
# 1: pipe VCF producing program to expos stdin.
# 2: calculate statistics with expos, reading VCF from stdin (-), output uncompressed VCF to stdout.
# 3, 4: statisically-backed flagging on distribution/clustering stats;
# flagging variants where observations are at least 2x as tightly clustered or spread as the background data,
# and the difference is statistically significant (P <= 0.05).
# 5: statistically-backed flagging on siginificant drops in alignment score on supporting reads in regions
# of low reference complexity;
# 6: heuristic/rule-of-thumb flagging on poor alignment score, and write to disk.
<some vcf producing command> |
./path/to/expos -u --ref ref.fa - my.bam |
bcftools filter --mode + -s CLUSTER -e'(INFO/QM1NN[1] <= -1.0 & INFO/QM1NN[2] < 0.05) | (INFO/TM1NN[1] <= -1.0 & INFO/TM1NN[2] < 0.05)' |
bcftools filter --mode + -s SPREAD -e'(INFO/QM1NN[1] >= 1 & INFO/QM1NN[2] < 0.05) | (INFO/TM1NN[1] >= 1 & INFO/TM1NN[2] < 0.05)' |
bcftools filter --mode + -s LOW_CMPLX -e'(INFO/MLAS[1] < 0 & INFO/MLAS[2] < 0.05 & INFO/RCMPLX < 150)' |
bcftools filter --mode + -s LOW_AS -e'(INFO/MLAS[0] < 0.93)' > my.flagged.vcf
```
