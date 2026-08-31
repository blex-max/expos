# Statistic validation

Standalone scripts that reproduce the calibration claims made
about `QRK`, `TJAC`, and `RCMPLX` in `docs/usage.md` and the paper. Not part
of the build. Requires `samtools` on `PATH` and `numpy`. The BAM must be
indexed, dupmarked, and carry `MC` tags (`samtools fixmate -m`).

## `null_calibration.py`

Evidence for the `z >= 3.0` cutoff: draws same-size subsamples from real
pileup columns and reports the
resulting effect-size distribution and false-positive rate by support size,
for both statistics.

```bash
./null_calibration.py sample.bam \
    21:20000000-20005000 22:30000000-30005000 \
    21:9500000-9505000 21:10500000-10505000
```

Example output: [`results/null_calibration.low_support.out.txt`](results/null_calibration.low_support.out.txt).

## `vcf_tails.py`

Behaviour on an annotated VCF: coverage, skip reasons, the upper
tail against the calibrated null, and effect size cross-tabulated against an
externally-defined difficult region plus `RCMPLX`.

```bash
expos --quiet calls.vcf ref.fa sample.bam > annotated.vcf
./vcf_tails.py annotated.vcf --stat TJAC --region 21:9000000-12000000
```

Example output: [`results/vcf_tails.germline.tjac.out.txt`](results/vcf_tails.germline.tjac.out.txt),
[`results/vcf_tails.germline.qrk.out.txt`](results/vcf_tails.germline.qrk.out.txt).

## `rcmplx_corroboration.py`

Evidence for the `RCMPLX < 20` cutoff: does the low tail concentrate in loci
with degraded alignment quality?

```bash
expos --quiet calls.vcf ref.fa sample.bam > annotated.vcf
./rcmplx_corroboration.py annotated.vcf sample.bam --cutoff 20 --bins 10
```

Example output: [`results/rcmplx_corroboration.somatic.out.txt`](results/rcmplx_corroboration.somatic.out.txt).

## `periodic_rle_corroboration.py`

A second, independent check on the same cutoff for RCMPLX.
`rcmplx_corroboration.py` shows low RCMPLX co-occurs with worse alignment
quality; this checks that it also co-occurs with actual sequence
repetitiveness, using an unrelated repeat detector (`periodic_rle`, bounded
period 1-5) on the same reference window.

```bash
expos --quiet calls.vcf ref.fa sample.bam > annotated.vcf
./periodic_rle_corroboration.py annotated.vcf ref.fa --flank 250 --bins 10
```

`--flank` must match whatever `expos` was run with, since the script
re-fetches the reference window from the FASTA directly. Reports two
correlations; `run_max` is noisier than `run_top10_mean` since one unrelated
repeat anywhere in a wide flank dominates it regardless of proximity to the
variant.

Example output: [`results/periodic_rle_corroboration.somatic.out.txt`](results/periodic_rle_corroboration.somatic.out.txt).
