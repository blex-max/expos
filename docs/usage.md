# Usage

## CLI

`expos` takes a VCF and annotates variants using the INFO field.
Usage of the tool is best described by the helptext:

```
Usage: expos [--help] [--version] [--seed SEED] [--flank SIZE]
             [--max-frag-len LEN] [--quiet] [--skip-filtered]
             [--background-sample PATH]...
             VCF REF ALN

Positional arguments:
  VCF                 input VCF/BCF of variants to annotate
  REF                 indexed reference genome FASTA
  ALN                 indexed alignment (BAM/CRAM) of the sample

Optional arguments:
  -h, --help          shows help message and exits
  -v, --version       prints version information and exits
  --seed SEED         random seed for the Monte-Carlo simulation [default: 24601]
  --flank SIZE        Size of reference sequence flanks to retrieve from
                      either side of variant, for use in calcuating
                      reference complexity. It is only necessary to modify
                      from the default value if low complexity windows more
                      distant from variant loci are likely to correlate with
                      artefacts given the sequencing protocol. [default: 250]
  --max-frag-len LEN  Upper bound on fragment size for a fragment to be
                      included in analysis. Useful to avoid confounding
                      template TJAC statistic with ambiguously mapped
                      fragments with improbable length. [default: 2000]
  -q, --quiet         suppress per-record warnings to stderr
  --skip-filtered     only analyse records where FILTER is PASS or . (unset)
  -b, --background-sample PATH
                      additional indexed alignment file/s from which to
                      merge data into Monte-Carlo background. Supporting
                      reads are always taken from the primary ALN only.
```

`VCF` may be `-` to read from stdin. The reference FASTA and
the alignment must both be indexed (`.fai`, and `.bai`/`.crai` of the same
name). The annotated VCF is written to stdout. A basic call would then look like:

```bash
expos my.vcf ref.fa my.bam > annotated.vcf
```

!!! note "Skipping records"
    `--skip-filtered` passes any record whose `FILTER` is set to a value
    other than `PASS` or `.` (unset) untouched - unanalysed records are
    are not removed from the output.


### Additional Background Samples

One or more additional indexed alignments can be merged into the Monte-Carlo background
with `-b`/`--background-sample`, repeating the flag once per file, e.g.
`-b other1.bam -b other2.bam`.
Supporting reads are always drawn from the primary `ALN` only; the extra samples contribute to
the background population against which statistics are simulated.

Background samples are only a valid source of statistical power if they were sequenced with the
same protocol as the primary sample. To guard against this, for each record the pileup of each
`-b`/`--background-sample` sample is evaluated against the primary sample before inclusion.
Therefore at the loci in question:

- read length within the sample must be homogenous;
- median read length must match the primary sample;
- median fragment (template) length must match the primary sample.

A source that fails any of these checks is excluded from the background for that record.
A warning is emitted to stderr per exclusion, unless `--quiet` has been passed.

The guards can't catch every possible way two populations might differ. In other words, they reduce
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
  Spatial clustering of mutant query positions against all reads:
  [0]standardised effect size (z-score) against the exact null;
  [1]one-sided Monte-Carlo p-value.
  Effect sizes greater than ~3.0 with a significant p-value may indicate a spurious variant.
  """>

##INFO=<
  ID=TJAC,
  Number=2,
  Type=Float,
  Description="""
  Graded pairwise overlap of mutant templates against all reads:
  [0]standardised effect size (z-score) against the exact null;
  [1]one-sided Monte-Carlo p-value.
  Effect sizes greater than ~3.0 with a significant p-value may indicate a spurious variant.
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
  Minimum LZ76 motif count over 100-base windows of the reference
  within a flank either side of the REF allele (lower means more
  repetitive; default 250 bases, see --flank; consult the
  expos_command header line for the value used this run).
  """>

##INFO=<
  ID=EXPOS_SKIP,
  Number=.,
  Type=String,
  Description="""
  Why expos produced no value, as '<scope>:<reason>' tokens.
  Scope is 'record' (whole record skipped) or a statistic ID.
  Reasons: not_biallelic, complex, insufficient_support, insufficient_background,
  heterogeneous_read_length, insufficient_reads_for_test, zero_variance,
  no_support, no_background, reference_too_short, reference_has_n.
  """>
```

Effect size is a standardised z-score relative to the null distribution. 0.0 indicates no difference from background, positive values indicate tighter clustering than background, in standard deviations away from the null. The null here is the distribution of the statistic over every equally-sized subsample of the background. Effect size is deterministic and does not move with `--seed`. The p-value is one-sided, giving the probability that a random subsample of the background would produce a statistic at least as extreme as the observed value. *p*-values are estimated by Monte-Carlo sampling; resolution is bounded by the draw count to 0.005. A large effect size combined with a small p-value, especially if found in a low complexity region as indicated by `RCMPLX`, may indicate a spurious variant call. A z-score of 3.0 is approximately the 1% false-positive point for each. It is not constant across support sizes; with only five supporting reads the true false-positive rate at 3.0 is nearer 2%, and low support is common in somatic calling.

Every output VCF is also tagged with metadata in the header: a `##source` line (program, version and UTC timestamp) and an `##expos_command` line recording the invocation.

!!! note "MLAS"
    `MLAS[0]` is equivalent to ASRD as may be familiar to some users - thresholding on this value may be inadvisable for indels since a decrease in alignment score is confounded with the presence of the indel itself.

## Missing Values

`expos` attempts to handle invalid records gracefully. When it cannot compute a statistic it writes a VCF missing value (`.`) for that field and records the reason in `EXPOS_SKIP` as one or more `<scope>:<reason>` tokens:

- **Whole-record skips** use the `record` scope. Records that are not biallelic (`record:not_biallelic`) and complex/untypeable records (`record:complex`) are passed through unannotated. Records skipped by `--skip-filtered` are not marked.
- **Per-statistic skips** use the statistic's ID as the scope, e.g. `QRK:insufficient_support`, `TJAC:zero_variance`, `RCMPLX:reference_has_n`.

!!! note "QRK"
  `QRK` may be suppressed with `QRK:heterogeneous_read_length` when the primary sample reads have markedly uneven lengths at a locus. Query position is an offset within a read and a mixed read-length population would confound statistics calculated.

!!! note "insufficent power"
  A locus with fewer than two supporting observations records `insufficient_support`; one whose background pool is smaller than ten, or smaller than twice the support, records `insufficient_background`.

  The absolute floor of ten total observation is to preserve p-value granularity. Both spatial statistics are sums over pairs within the supporting set, so at the minimum support of two the observed value is a single pair and the null is the distribution over every pair the background can supply. Ten background observations give 45 such pairs and a p-value granular to about 0.02; five would give 10 pairs and a granularity of 0.1.

## Filtering Pipeline Examples

The following are some examples of how the annotations made by `expos` might be used for variant filtering downstream.

The annotations are meant to be read together rather than in isolation. A variant tripping several is therefore much more suspect than one tripping any single test. The examples below add flags separately rather than combining them into one expression so that the count remains visible in the `FILTER` column.

### **Broad Flagging Pipeline**

This is a fairly non-specific example showing the breadth of what
one might do with the information encoded by expos — it's not
strictly a recommendation, though it is statistically defensible.

```bash
# example pipeline - Add some soft flags in the FILTER column
# (or alternately, subset entirely with bcftools view instead of filter)

# command by command:
# 1: pipe VCF producing program to expos stdin.
# 2: calculate statistics with expos, reading VCF from stdin (-); output to stdout.
# note that for brevity no additional background sample is provided, but providing one
# via -b/--background-sample can add a lot of statistical power if one is available.
# 3, 4: statisically-backed flagging on distribution/clustering stats;
# flagging variants whose clustering sits about three standard deviations above
# the background and is statistically significant (P <= 0.05).
# 6: heuristic/rule-of-thumb on poor alignment score on supporting reads in regions
# of low reference complexity;
# 7: heuristic/rule-of-thumb flagging on poor alignment score
# and > write to disk.
./path/to/expos my.vcf ref.fa my.bam |
bcftools filter -Ov \
  --mode + \
  -s QPOS_CLUSTER \
  -e'(INFO/QRK[0] >= 3.0 & INFO/QRK[1] < 0.05)' |
bcftools filter -Ov \
  --mode + \
  -s TEMPLATE_OVERLAP \
  -e'(INFO/TJAC[0] >= 3.0 & INFO/TJAC[1] < 0.05)' |
bcftools filter -Ov \
  --mode + \
  -s POOR_ALN_REG \
  -e'(INFO/MLAS[1] < 0.93 & INFO/RCMPLX < 20)' |
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
./path/to/expos my.vcf ref.fa my.bam |
bcftools filter -Oz \
  --mode + \
  -s LOW_CMPLX_CLUSTER \
  -e'INFO/QRK[0] >= 3.0 & INFO/QRK[1] < 0.05 & INFO/RCMPLX < 20' > my.flagged.vcf.gz
```

At the cost of missing more generic variants with spurious looking spatial properties.

### **Adjusting Thresholds**

Thresholds for p-value and effect size can be tuned:

```bash
# relaxed p-val, and an effect size well beyond the calibrated 1% point
# an example of the concept, again not a recommendation per se
./path/to/expos my.vcf ref.fa my.bam |
bcftools filter -Oz \
  --mode + \
  -s QPOS_CLUSTER_2 \
  -e'INFO/QRK[0] >= 4.5 & INFO/QRK[1] < 0.1' > my.flagged.vcf.gz
```

## Thresholding on Complexity (`RCMPLX`)

Unlike `QRK`/`TJAC`, `RCMPLX` has no resampling null to calibrate against: it
is a property of the reference sequence alone, not of the read data, so
there is no false-positive rate to derive a cutoff from. `RCMPLX` reports
the *minimum* LZ76 phrase count over the 100-base windows tiling the flank,
not a mean — a single genuinely repetitive window near the variant is real,
actionable signal on its own, and a mean over many windows across a wide
flank dilutes exactly that signal almost to invisibility. `RCMPLX < 20` is a
practical operating point rather than a calibrated threshold, justified
instead by real-data corroboration, checked two independent ways:

1. **Downstream alignment quality.** On a real annotated somatic call set,
   binning calls by `RCMPLX` and comparing against markers `RCMPLX` never
   sees (the read-alignment-based `CLPM`/`ASRD`/`ASMD` fields from the
   calling pipeline, and expos's own `MLAS`) shows a real, monotonic
   degradation in `ASRD`/`ASMD`/`MLAS` as `RCMPLX` falls, steepest in the
   bottom two or three deciles of real calls, which is where 20 sits.
   `CLPM` shows no discriminating power either way in this call set.
2. **Independent sequence-structure confirmation.** A separate, purely
   reference-sequence-based tandem-repeat detector (no read data, no LZ76)
   run over the same window agrees: the length of the longest periodic
   repeats found correlates with `RCMPLX` (Spearman ≈ -0.5, monotonic
   across deciles), confirming low `RCMPLX` reflects genuine local
   repetitiveness rather than being an artefact correlated with something
   else.

This threshold was also checked at flank sizes from 200 to 400 bases before
the default settled at 250; the degradation pattern barely moved across
that range, so the cutoff isn't sensitive to the exact flank chosen.

!!! note "Caveat"
    Treat as a continuous risk factor to combine with other statistics, not a clean binary classifier of "difficult" vs "normal" sequence.

