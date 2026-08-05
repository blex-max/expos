# Usage

## CLI

`expos` takes a VCF and annotates variants using the INFO field.
Usage of the tool is best described by the helptext:

```
Usage: expos [--help] [--version] [--seed SEED] [--uncompressed] [--quiet]
             [--skip-filtered] [--additional-background-samples PATH...]
             [--debug] VCF REF ALN

Positional arguments:
  VCF                 input VCF/BCF of variants to annotate
  REF                 indexed reference genome FASTA
  ALN                 indexed alignment (BAM/CRAM) of the sample

Optional arguments:
  -h, --help          shows help message and exits
  -v, --version       prints version information and exits
  --seed SEED         random seed for the Monte-Carlo simulation [default: 24601]
  -u, --uncompressed  write uncompressed VCF (default: bgzip-compressed)
  -q, --quiet         suppress per-record warnings to stderr
  --skip-filtered     only analyse records where FILTER is PASS or . (unset)
  --bg, --additional-background-samples PATH...
                      additional indexed alignment file(s) whose reads are
                      merged into the Monte-Carlo background. Supporting reads
                      are always taken from the primary ALN only.
  --debug             enable debug logging to stderr
```

The three positional arguments are **required** and consumed in the order
`VCF REF ALN`. `VCF` may be `-` to read from stdin. The reference FASTA and
the alignment must both be indexed (`.fai`, and `.bai`/`.crai` of the same
name). The annotated VCF is written to stdout. A basic call would then look
like:

```bash
expos my.vcf ref.fa my.bam > annotated.vcf
```

!!! note "Skipping records"
    `--skip-filtered` passes any record whose `FILTER` is set to a value
    other than `PASS` or `.` (unset) untouched - unanalysed records are
    are not removed from the output.


### Additional Background Samples

One or more further indexed alignments (typically a matched normal) can
be merged into the Monte-Carlo background with `--bg`, e.g. `--bg other1.bam other2.bam`.
Supporting reads are always drawn from the primary `ALN` only; the extra samples contribute to
the background population against which clustering is simulated.

Background samples are only a valid source of statistical power if they were sequenced with the
same protocol as the primary sample. To guard against this, for each record the pileup of each
`--bg` sample is evaluated against the primary sample before inclusion. Therefore:

- its own read lengths must be internally consistent (not itself a mix of very different read
  lengths);
- its median read length must match the primary sample's;
- its median fragment (template) length must match the primary sample's.

A source that fails any of these checks is excluded from the background for that record.
Insufficient reads or templates to evaluate a check also excludes the source,
rather than assuming validity. A warning is emitted to stderr per exclusion (unless `--quiet`),
and the total number of exclusions across the whole run is reported in the end-of-run summary.

The guards can't catch every possible way two populations might differ. In other words they reduce
the risk of an inappropriate background rather than eliminating it.


## Annotations Made

For a full guide to the theory behind these annotations see [Concepts](concepts.md).

These are the header lines from an output VCF describing the INFO fields added. The `[]` notation indicates which element of the array holds the data in question where the INFO field added is an array. Arrays are 0-indexed, matching how `bcftools` addresses them (so `INFO/QRK[0]` is the effect size and `INFO/QRK[1]` the p-value).

```
##INFO=<
  ID=QRK,
  Number=2,
  Type=Float,
  Description="""
  Monte-Carlo results for spatial clustering of mutant query positions:
  [0]standardised effect size (z-score) versus simulation against all reads;
  [1]one-sided p-value.
  Effect sizes greater than ~2.0 with a significant p-value may indicate a spurious variant.
  """>

##INFO=<
  ID=TRK,
  Number=2,
  Type=Float,
  Description="""
  Monte-Carlo results for spatial clustering of mutant template endpoints:
  [0]standardised effect size (z-score) versus simulation against all reads;
  [1]one-sided p-value.
  Effect sizes greater than ~2.0 with a significant p-value may indicate a spurious variant.
  """>

##INFO=<
  ID=MLAS,
  Number=2,
  Type=Float,
  Description="""
  Median read-length-normalised alignment scores:
  [0]of reads supporting the variant;
  [1]of all reads covering the variant site.
  """>

##INFO=<
  ID=RCMPLX,
  Number=1,
  Type=Float,
  Description="""
  Mean 100-base window complexity (Lempel-Ziv 76 entropy rate) of the
  reference region spanned by supporting templates.
  """>

##INFO=<
  ID=EXPOS_SKIP,
  Number=.,
  Type=String,
  Description="""
  Why expos produced no value, as '<scope>:<reason>' tokens.
  Scope is 'record' (whole record skipped) or a statistic ID.
  Reasons: multiallelic, complex, insufficient_support, insufficient_background,
  heterogeneous_read_length, zero_variance,
  no_support, no_background, reference_too_short, reference_has_n.
  """>
```

Effect size is a standardised z-score relative to the simulated null distribution: 0.0 indicates no difference from background, positive values indicate tighter clustering than background (negative values looser), and values beyond ~2.0 are roughly two standard deviations from the mean of the null. The p-value is one-sided, giving the probability that a random subsample of the background would produce a statistic at least as extreme as the observed value. In blunt summary, effect sizes greater than ~2.0 combined with small p-values, especially if found in a low complexity region as indicated by `RCMPLX`, may indicate a spurious variant call.

Every output VCF also carries provenance in its header: a `##source` line (program, version and UTC timestamp) and an `##expos_command` line recording the invocation.

!!! note "MLAS"
    `MLAS[0]` is equivalent to ASRD as may be familiar to some users - thresholding on this value may be inadvisable for indels since a decrease in alignment score is confounded with the presence of the indel itself.

## Missing Values

`expos` attempts to handle invalid records gracefully. When it cannot compute a statistic it writes a VCF missing value (`.`) for that field and records the reason in `EXPOS_SKIP` as one or more `<scope>:<reason>` tokens:

- **Whole-record skips** use the `record` scope. Multiallelic records (`record:multiallelic`) and complex/untypeable records (`record:complex`) are passed through unannotated. Records skipped by `--skip-filtered` are not marked.
- **Per-statistic skips** use the statistic's ID as the scope, e.g. `QRK:insufficient_support`, `TRK:zero_variance`, `RCMPLX:reference_has_n`.

!!! note "QRK"
  `QRK` may be suppressed with `QRK:heterogeneous_read_length` when the primary sample's reads at a locus have markedly uneven lengths, since query position is an offset within a read and a mixed read-length population confounds it. `TRK` has no equivalent per-statistic reason; read-length heterogeneity doesn't apply.

!!! note "insufficent power"
  If insufficient reads are available to perform testing, expos skip will record the `insufficient_background` reason.


## Filtering Pipeline Examples

The following are some examples of how the annotations made by `expos` might be used for variant filtering downstream.

### **Broad Flagging Pipeline**

This is a fairly non-specific example showing the breadth of what
one might do with the information encoded by expos — it's not
strictly a recommendation, though it is statistically defensible.

```bash
# example pipeline - Add some soft flags in the FILTER column
# (or alternately, subset entirely with bcftools view instead of filter)

# command by command:
# 1: pipe VCF producing program to expos stdin.
# 2: calculate statistics with expos, reading VCF from stdin (-), output uncompressed VCF to stdout.
# note that for brevity no additional background sample is provided, but providing one
# (e.g. a matched normal) via --bg can add a lot of statistical power if one is available.
# 3, 4: statisically-backed flagging on distribution/clustering stats;
# flagging variants where observations are at least 4x as tightly clustered as the background
# and the difference is statistically significant (P <= 0.05).
# 6: heuristic/rule-of-thumb on poor alignment score on supporting reads in regions
# of low reference complexity;
# 7: heuristic/rule-of-thumb flagging on poor alignment score
# and > write to disk.
./path/to/expos -u my.vcf ref.fa my.bam |
bcftools filter -Ov \
  --mode + \
  -s QPOS_CLUSTER \
  -e'(INFO/QRK[0] >= 2.0 & INFO/QRK[1] < 0.05)' |
bcftools filter -Ov \
  --mode + \
  -s TEMPLATE_CLUSTER \
  -e'(INFO/TRK[0] >= 2.0 & INFO/TRK[1] < 0.05)' |
bcftools filter -Ov \
  --mode + \
  -s POOR_ALN_REG \
  -e'(INFO/MLAS[1] < 0.93 & INFO/RCMPLX < 1.5)' |
bcftools filter -Oz \
  --mode + \
  -s LOW_SUPPORTING_AS \
  -e'(INFO/MLAS[0] < 0.93)' > my.flagged.vcf.gz
```

### **Targeted Approach**

A more targeted approach can inform you as to particular
scenarios that may be strongly associated with false positive variants. This example
looks specifically for clustered variants in low complexity regions

```bash
./path/to/expos -u my.vcf ref.fa my.bam |
bcftools filter -Oz \
  --mode + \
  -s LOW_CMPLX_CLUSTER \
  -e'INFO/QRK[0] >= 2.0 & INFO/QRK[1] < 0.05 & INFO/RCMPLX < 1.5' > my.flagged.vcf.gz
```

At the cost of missing more generic variants with spurious looking spatial properties.

### **Adjusting Thresholds**

Thresholds for p-value and effect size can be tuned:

```bash
# relaxed p-val, very large effect size (8x as clustered)
# an example of the concept, again not a recommendation per se
./path/to/expos -u my.vcf ref.fa my.bam |
bcftools filter -Oz \
  --mode + \
  -s QPOS_CLUSTER_2 \
  -e'INFO/QRK[0] >= 3.0 & INFO/QRK[1] < 0.1' > my.flagged.vcf.gz
```

## Thresholding on Complexity (`RCMPLX`)

`RCMPLX` is only meaningful relative to what "low complexity" looks like across the genome as a whole. To calibrate a sensible cutoff, `estimate-entropy` was run genome-wide against GRCh38, and separately on the low-complexity regions identified by the GIAB v3.1 genome stratifications (the same homopolymer and tandem-repeat region sets GIAB itself uses to stratify benchmarking performance):

| Region                                              | Min   | Median | Mean  | SD    |
|------------------------------------------------------|------:|-------:|------:|------:|
| Whole genome (100bp windows)                          | 0.199 | 1.927  | 1.919 | 0.150 |
| Homopolymers (&gt;6bp, imperfect &gt;10bp, 5bp slop)   | 0.567 | 1.472  | 1.497 | 0.266 |
| Tandem repeats (&gt;100bp, 5bp slop)                   | 0.185 | 1.472  | 1.447 | 0.424 |
| Homopolymers + tandem repeats (union)                  | 0.250 | 1.466  | 1.464 | 0.294 |
| Outside homopolymers                                   | 0.380 | 1.898  | 1.874 | 0.159 |
| Outside homopolymers + tandem repeats                  | 0.862 | 1.905  | 1.898 | 0.121 |

GIAB's difficult regions cluster tightly around a median of ~1.47, over half a bit below the rest of the genome (median ~1.9-1.93) and with roughly double the spread. This is the basis for the `RCMPLX < 1.5` heuristic used in the filtering examples above: it sits almost exactly at the median of GIAB's difficult-region strata, and comfortably below the bulk of the rest of the genome.

!!! note "Not a hard boundary"
    The two distributions overlap substantially (the "outside" regions still have a minimum as low as 0.38-0.86), so `RCMPLX` is best treated as a continuous risk factor to combine with other statistics, not a clean binary classifier of "difficult" vs "normal" sequence.

