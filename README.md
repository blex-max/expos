# expos

[![CI](https://github.com/blex-max/expos/actions/workflows/ci.yml/badge.svg)](https://github.com/blex-max/expos/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-online-blue)](https://blex-max.github.io/expos/)
[![C++23](https://img.shields.io/badge/C%2B%2B-23-blue)](https://en.cppreference.com/w/cpp/23)
[![status: alpha](https://img.shields.io/badge/status-alpha-orange)](https://github.com/blex-max/expos)

**FULL DOCUMENTATION AVAILABLE [HERE](https://blex-max.github.io/expos/)**

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
 - More extensive documentation
 - Add clip fraction/CLPM

## Installation

```bash
 # from the cloned repo
 mkdir build
 cd build
 cmake .. -DCMAKE_BUILD_TYPE=Release
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

  -h, --help                   Print usage
  -l, --expected-read-len arg  Sequencing read length. Default: 150
  -r, --ref arg                Alignment Reference Fasta for optionally
                               adding reference complexity to statistics.
  -n, --normal arg             Alignment for use as additional background
                               data for simulation. May be passed multiple
                               times.
      --normal-only            Use only reads from the provided normal as
                               background data, excluding reads from the
                               sample
  -i, --include-records arg    Only operate on VCF records with this value
                               present in FILTER. e.g. -i PASS. May be
                               passed multiple times.
  -e, --exclude-records arg    Only operate on VCF records without this
                               value present in FILTER. May be passed
                               multiple times.
  -I, --include-reads arg      SAM flag: only include reads with all of
                               these bits set. Set 0 to disable. Default: 3
  -E, --exclude-reads arg      SAM flag: exclude reads with any of these
                               bits set. Default: 3852
  -t, --tsv arg                Write a tsv of extended statistics to file
                               specified.
  -u, --uncompressed           output uncompressed VCF
      --seed arg               Set random seed. Default: 24601
      --uniform                additionally simulate against uniform null
                               model for query position, and add result to
                               --tsv output. For assessment of correlation
                               with simulation against all-reads null.
      --assess-microhomology   additionally assess STR and homopolymer
                               content of reference regions, and add result
                               to --tsv output. For assessment of
                               correlation with drop in LZ.
      --debug                  Enable debug logging
```
basic usage then looks like:
```bash
expos my.vcf my.bam
```


## Assessments made

Assessments of spatial clustering of mutant bases have been used to filter false-positive mutations which may otherwise be difficult to consistently detect (Ellis 2021, Cambpell Lab, GATK ReadPosRankSum, others). The underlying hypothesis is that mutant reads should be drawn from the same spatial distribution as reference reads; if mutant reads differ significantly from reference reads, then the spatial process producing those reads deviates from the spatial process producing the reference reads. This may indicate that a non-biological process, or sequencing artefact, is responsible for the mutant reads since it would not be expected that mutant reads exhibit a unique preference for a particular region.

These are the header lines from an output VCF describing the INFO fields added. The [] notation is used to indicate which element of the array holds the data in question where the INFO field added is an array.

```bash
##INFO=<
  ID=QRK,
  Number=2,
  Type=Float,
  Description="""
  Array detailing monte-carlo simulation results for Ripley's K on mutant query position:
  [1]log2 ratio effect size from comparisons to simulation against all reads;
  [2]two-sided P-value from comparisons to simulation against all reads;
  """>

##INFO=<
  ID=TRK,
  Number=2,
  Type=Float,
  Description="""
  Array detailing Monte-Carlo simulation results for Ripley's K on endpoints of mutant templates:
  [1]log2 ratio effect size from comparisons to simulation against all reads;
  [2]two-sided P-value from comparisons to simluation against all reads
  """>

##INFO=<
  ID=RCMPLX,
  Number=1,
  Type=Integer,
  Description="""
  Mean 100-base window complexity (Lempel-Ziv estimated entropy rate) of
  the reference region spanned by supporting templates, scaled by x100
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

log2-fold change scales such that no effect is 0.0, 1.0 means the statistic is 2x on supporting data compared background, 2.0 == 4x. Practically this means that effect sizes greater than 0 indicate tighter clustering of observations as compared to background.

MLAS[0] is equivalent to ASRD as may be familiar to some users - thresholding on this value may be inadvisable for indels since a decrease in alignment score is confounded with the presence of the indel itself.

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
# 3, 4: statisically-backed flagging on distribution/clustering stats;
# flagging variants where observations are at least 2x as tightly clustered as the background
# and the difference is statistically significant (P <= 0.05).
# 6: heuristic/rule-of-thumb on poor alignment score on supporting reads in regions
# of low reference complexity;
# 7: heuristic/rule-of-thumb flagging on poor alignment score
# and > write to disk.
./path/to/expos -u --ref ref.fa my.vcf my.bam |
bcftools filter -Ov \
  --mode + \
  -s QPOS_CLUSTER \
  -e'(INFO/QRK[0] >= 1.0 & INFO/QRK[1] < 0.05)' |
bcftools filter -Ov \
  --mode + \
  -s TEMPLATE_CLUSTER \
  -e'(INFO/TRK[0] >= 1.0 & INFO/TRK[1] < 0.05)' |
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
  -e'INFO/QRK[0] >= 1.0 & INFO/QRK[1] < 0.05 & INFO/RCMPLX < 150' > my.flagged.vcf.gz
```
at the cost of missing more generic variants with spurious looking spatial properties.

P-values and effect sizes can be modified:
```bash
# relaxed p-val, very large effect size (8x as clustered)
# an example of the concept, again not a recommendation per se
./path/to/expos -u --ref ref.fa my.vcf my.bam |
bcftools filter -Oz \
  --mode + \
  -s QPOS_CLUSTER_2 \
  -e'INFO/QRK[0] >= 3.0 & INFO/QRK[1] < 0.1' > my.flagged.vcf.gz
```

Since the p-values are returned are two-tailed, you can also look
at deviation in the other direction - though it is not intuitively obvious
that this would be associated with a false positive variant.
```bash
# at least twice as spread as expected, and statistically significant
./path/to/expos -u --ref ref.fa my.vcf my.bam |
bcftools filter -Oz \
  --mode + \
  -s QPOS_SPREAD \
  -e'INFO/QRK[0] <= -1.0 & INFO/QRK[1] < 0.05' > my.flagged.vcf.gz
```


## Extras

This repo also contains a daughter tool for estimating the entropy rate of strings
in the same way as is done for assessing reference complexity. Useful if you're interested
in assessing the complexity of sequence data from any source (or any string at all!).

You can turn on compliation of `estimate-entropy` via
```
  cmake .. -DMAKE_DAUGHTER=ON
```
and then proceeding with compilation as normal.
See `./estimate-entropy --help` for usage.
