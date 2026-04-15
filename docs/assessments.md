# Assessments

Assessments of spatial clustering of mutant bases have been used to filter false-positive mutations which may otherwise be difficult to consistently detect (Ellis 2021, Cambpell Lab, GATK ReadPosRankSum, others). The underlying hypothesis is that mutant reads should be drawn from the same spatial distribution as reference reads; if mutant reads differ significantly from reference reads, then the spatial process producing those reads deviates from the spatial process producing the reference reads. This may indicate that a non-biological process, or sequencing artefact, is responsible for the mutant reads since it would not be expected that mutant reads exhibit a unique preference for a particular region.

These are the header lines from an output VCF describing the INFO fields added. The `[]` notation is used to indicate which element of the array holds the data in question where the INFO field added is an array.

```
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

log2-fold change scales such that no effect is 0.0, 1.0 means the statistic is 2x on supporting data compared to background, 2.0 == 4x. Practically this means that effect sizes greater than 0 indicate tighter clustering of observations as compared to background.

`MLAS[0]` is equivalent to ASRD as may be familiar to some users - thresholding on this value may be inadvisable for indels since a decrease in alignment score is confounded with the presence of the indel itself.
