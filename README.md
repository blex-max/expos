# expos

Spatial information extraction for variant sites from alignment data,
plus comparison to background via monte-carlo simulation.
Useful for inspecting and flagging false postive variants.

Core functionality present but niceties and guard rails are not.

NOT CURRENTLY EXTENDED TO INDELS. The concept applies,
and the logic is largely in place,
but at this stage the program has only been applied to SNVs
and almost certainly won't work for indels. I'll add
indel support in Jan 2026.

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
EXtract POSitional data and statistics from alignment at
VCF variant sites, and encode them as INFO fields to VCF.
Requires the presence of .(b/cr)ai indexes of the same name
as the relevant alignment. Annotated VCF to stdout. See
README or output VCF header for descriptions of fields
added.

Usage:
  expos [OPTION...] <VCF/BCF (- for stdin)> <ALN.(b/cr)am>

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
basic usage then looks like:
```bash
expos my.vcf my.bam
```


## Annotations made

These are the header lines from an output VCF describing the INFO fields added. The [] notation is used to indicate which element of the array holds the data in question where the INFO field added is an array.

```
##INFO=<ID=KC,Number=1,Type=Integer,Description="Kolmogorov Complexity of region spanned by supporting templates, scaled by x100">
##INFO=<ID=TM1NN,Number=3,Type=Float,Description="[0]Median nearest neighbour distance of template endpoints from read pairs supporting variant; [1]log2 ratio effect size and [2]P-value against background, from monte-carlo simulation">
##INFO=<ID=QM1NN,Number=3,Type=Float,Description="[0]Median nearest neighbour distance of variant query position; [1]log2 ratio effect size and [2]P-value against background, from monte-carlo simulation">
##INFO=<ID=MLAS,Number=3,Type=Float,Description="[0]Median read-Length-normalised Alignment Score (AS) of reads supporting variant;[1]delta (supporting - background) effect size and [2]P-value against background, from monte-carlo simulation">
```

log2 effect sizes scale such no effect is 0, -1 means the statistic is 1/2 on supporting data compared background, -2 1/4, whereas an effect size of 1 means the statistic 2x on supporting data compared to background, 2 4x. Practically this means that effect sizes below 0 indicate tighter clustering of observations as compared to background.

The effect size for median length-normalised alignment score is simply reported as difference between the statisic as calculated on the supporting data and the mean of all simulated calculations. Therefore a lower alignment score on supporting reads as compared to background data results in a negative effect size. Note that the presence of a variant will by definition lower alignment score so very small but significant (by p value) effects do not necessarily indicate a spurious variant.

MLAS is equivalent to ASRD as may be familiar to some users.

## Example

```bash
 # example pipeline - Add some soft flags in the FILTER column
 # (or alternately, subset entirely with bcftools view instead of filter)

# line by line:
# 1: pipe VCF producing program to expos stdin.
# 2: calculate statistics with expos, reading VCF from stdin (-), output uncompressed VCF to stdout.
# 3: statisically-backed flagging on distribution/clustering stats;
# flagging variants where observations are at least 2x as tightly clustered as the background data,
# and the difference is statistically significant (P <= 0.05).
# 4: statistically-backed threshold flagging on alignment score;
# flagging variants where MLAS < 0.93, and the difference between supporting reads
# and background is statistically significant (P <= 0.05)
# but ignoring very small effects.
# 5: conservative heuristic/rule-of-thumb flagging on poor alignment score in regions of low complexity, and write to disk.
<some vcf producing command> |
./path/to/expos -u --ref ref.fa - my.bam |
bcftools filter --mode + -s SPATIAL -e'(INFO/QM1NN[1] <= -1.0 & INFO/QM1NN[2] <= 0.05) | (INFO/TM1NN[1] <= -1.0 & INFO/TM1NN[2] <= 0.05)' |
bcftools filter --mode + -s LOW_AS -e'(INFO/MLAS[0] < 0.93 & INFO/MLAS[1] < -0.05 & INFO/MLAS[2] <= 0.05)' |
bcftools filter --mode + -s LOW_CMPLX -e'(INFO/MLAS[0] < 0.93 & INFO/KC < 150)' > my.flagged.vcf
```
