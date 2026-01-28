# expos

Statistically-backed VCF flagging calculating effect sizes
and p-values for spatial properties of somatic mutations via
monte-carlo simulation. Compares mutant reads (supporting observations)
with the set of all reads (background), optionally including reads
from a matched normal.

Useful for inspecting and flagging false postive variants caused
by a variety of processes most commonly associated with a high
number of PCR cycles. Builds on Ellis et al. 2021, GATK ReadPosRankSum, amongst others.

For both SNVs and small indels. MNV handling logic is present,
but largely untested.

Alpha Software!
Core functionality present but niceties and guard rails are not.
Please report any bugs and ask any questions!

#### definite TODOs
 - report correct position in error messages - chrom and position are both
   off by one! (this affects error messages only, and has no impact on the output)
 - More extensive documentation
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
      --normal-only       Use only reads from the provided normal as
                          background data, excluding non-supporting reads
                          from the sample
  -r, --ref arg           Alignment Reference Fasta for optionally adding
                          reference complexity to statistics.
      --seed arg          Set random seed. Default: 24601
  -u, --uncompressed      output uncompressed VCF
```
basic usage then looks like:
```bash
expos my.vcf my.bam
```


## Assessments made

Assessments of spatial clustering of mutant bases have been used to filter false-positive mutations which may otherwise be difficult to consistently detect (Ellis 2021, Cambpell Lab, GATK ReadPosRankSum, others). The underlying hypothesis is that mutant reads should be drawn from the same spatial distribution as reference reads; if mutant reads differ significantly from reference reads, then the spatial process producing those reads deviates from the spatial process producing the reference reads. This may indicate that a non-biological process, or sequencing artefact, is responsible for the mutant reads since it would not be expected that mutant reads exhibit a unique preference for a particular region.

expos implements a nearest-neighbour algorithm (Cover, 1967) on two spatial properties of the set of mutant reads; the query position of the mutation on each read, and the endpoints of the inferred template from which each read was amplified. For each property, expos finds the distance to the single closest neighbour for each read, and reports the upper quartile of the set of these distances. Since the unit of these metrics is in sequence bases, they are readily interpretable as descriptive statistics – what is the average distance in bases to the closest neighboring observation? The inclusion of empirical two-tailed P values, and log2-fold change effect sizes, extend the metrics beyond descriptive statistics and allow for defensible flagging of variants on a sound statistical basis.

These are the header lines from an output VCF describing the INFO fields added. The [] notation is used to indicate which element of the array holds the data in question where the INFO field added is an array.

```bash
##INFO=<
  ID=Q1NN_UQ,
  Number=5,
  Type=Float,
  Description="""
  Array detailing the upper quartile of nearest neighbour distances
  between mutant query positions and monte-carlo simulation results:
  [0]calculated statistic;
  [1]log2 ratio effect size from comparisons to simulation against all reads;
  [2]two-sided P-value from comparisons to simulation against all reads;
  [3]log2 ratio effect size from comparisons to simulation against uniform distribution;
  [4]two-sided P-value from comparisons to simulation against uniform distribution
  """>

##INFO=<
  ID=T1NN_UQ,
  Number=3,
  Type=Float,
  Description="""
  Array detailing the upper quartile of nearest neighbour distances
  between endpoints of supporting templates and monte-carlo simulation results:
  [0]calculated statistic;
  [1]log2 ratio effect size from comparisons to simulation against all reads;
  [2]two-sided P-value from comparisons to simluation against all reads
  """>

##INFO=<
  ID=RCMPLX,
  Number=1,
  Type=Integer,
  Description="""
  Complexity (Lempel-Ziv estimated entropy rate)
  of region spanned by supporting templates, scaled by x100
  """>

##INFO=<
  ID=MLAS,
  Number=2,
  Type=Float,
  Description="""
  Array of median read-length normalised alignment scores:
  [0]of reads supporting variant,
  [1]of all queried reads covering the variant location in the sample alignment
  """>
```

log2-fold change scales such that no effect is 0, -1 means the statistic is 1/2 on supporting data compared background, -2 1/4, whereas an effect size of 1 means the statistic is 2x on supporting data compared to background, 2 4x. Practically this means that effect sizes below 0 indicate tighter clustering of observations as compared to background.

MLAS[0] is equivalent to ASRD as may be familiar to some users.

## Example

This is fairly non-specific example showing the breadth of what
one might do with the information encoded by expos - it's not
strictly a recommendation, though it is statistically defensible.

```bash
 # example pipeline - Add some soft flags in the FILTER column
 # (or alternately, subset entirely with bcftools view instead of filter)

# command by command:
# 1: pipe VCF producing program to expos stdin.
# 2: calculate statistics with expos, reading VCF from stdin (-), output uncompressed VCF to stdout.
# note that for brevity no normal is provided, but providing a normal can add a lot of statistical power
# if an appropriate normal is available.
# 3, 4, 5: statisically-backed flagging on distribution/clustering stats;
# flagging variants where observations are at least 2x as tightly clustered or spread as the background data,
# and the difference is statistically significant (P <= 0.05).
# 6: heuristic/rule-of-thumb on poor alignment score on supporting reads in regions
# of low reference complexity;
# 7: heuristic/rule-of-thumb flagging on poor alignment score
# and > write to disk.
./path/to/expos -u --ref ref.fa my.vcf my.bam |
bcftools filter -Ov \
  --mode + \
  -s QPOS_CLUSTER \
  -e'(INFO/Q1NN_UQ[1] <= -1.0 & INFO/Q1NN_UQ[2] < 0.05)' |
bcftools filter -Ov \
  --mode + \
  -s TEMPLATE_CLUSTER \
  -e'(INFO/T1NN_UQ[1] <= -1.0 & INFO/T1NN_UQ[2] < 0.05)' |
bcftools filter -Ov \
  --mode + \
  -s QPOS_NON_UNIFORM \
  -e'(INFO/Q1NN_UQ[3] <= -1.0 & INFO/Q1NN_UQ[4] < 0.05)' |
bcftools filter -Ov \
  --mode + \
  -s POOR_ALN_REG \
  -e'(INFO/MLAS[1] < 0.93 & INFO/RCMPLX < 150)' |
bcftools filter -Oz \
  --mode + \
  -s LOW_SUPPORTING_AS \
  -e'(INFO/MLAS[0] < 0.93)' > my.flagged.vcf.gz
```

A more targeted approach can inform you as to particular
scenarios that may be strongly associated with false postive variants:

```bash
./path/to/expos -u --ref ref.fa my.vcf my.bam |
bcftools filter -Oz \
  --mode + \
  -s LOW_CMPLX_CLUSTER \
  -e'INFO/Q1NN_UQ[1] <= -1.0 & INFO/Q1NN_UQ[2] < 0.05 & INFO/RCMPLX < 150' > my.flagged.vcf.gz
```
at the cost of missing more generic variants with spurious looking spatial properties

P-values and effect sizes can be modified
```bash
# relaxed p-val, very large effect size (8x as clustered)
# an example of the concept, again not a recommendation per se
./path/to/expos -u --ref ref.fa my.vcf my.bam |
bcftools filter -Oz \
  --mode + \
  -s QPOS_CLUSTER_2 \
  -e'INFO/Q1NN_UQ[1] <= -3.0 & INFO/Q1NN_UQ[2] < 0.1' > my.flagged.vcf.gz
```

Since the p-values are returned are two-tailed, you can also look
at deviation in the other direction - though it is not intuitively obvious
that this would be associated with a false positive variant
```bash
# at least twice as spread as expected, and statistically significant
./path/to/expos -u --ref ref.fa my.vcf my.bam |
bcftools filter -Oz \
  --mode + \
  -s QPOS_SPREAD \
  -e'INFO/Q1NN_UQ[1] >= 1.0 & INFO/Q1NN_UQ[2] < 0.05' > my.flagged.vcf.gz
```
