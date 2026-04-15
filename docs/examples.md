# Examples

## **Broad Flagging Pipeline**

This is a fairly non-specific example showing the breadth of what
one might do with the information encoded by expos — it's not
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

## **Targeted Approach**

A more targeted approach can inform you as to particular
scenarios that may be strongly associated with false positive variants:

```bash
./path/to/expos -u --ref ref.fa my.vcf my.bam |
bcftools filter -Oz \
  --mode + \
  -s LOW_CMPLX_CLUSTER \
  -e'INFO/QRK[0] >= 1.0 & INFO/QRK[1] < 0.05 & INFO/RCMPLX < 150' > my.flagged.vcf.gz
```

At the cost of missing more generic variants with spurious looking spatial properties.

## **Adjusting Thresholds**

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

## **Deviation in the Other Direction**

Since the p-values returned are two-tailed, you can also look
at deviation in the other direction — though it is not intuitively obvious
that this would be associated with a false positive variant.

```bash
# at least twice as spread as expected, and statistically significant
./path/to/expos -u --ref ref.fa my.vcf my.bam |
bcftools filter -Oz \
  --mode + \
  -s QPOS_SPREAD \
  -e'INFO/QRK[0] <= -1.0 & INFO/QRK[1] < 0.05' > my.flagged.vcf.gz
```
