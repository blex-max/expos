# Statistic Validation

Empirical backing for the `QRK` and `TJAC` effect-size cutoffs, and for the
replacement of the former `TRK` statistic. Measured against `PD23344a.21-22.bam`,
a dupmarked GRCh37 WGS sample (151bp reads) restricted to chromosomes 21 and 22.
Reproduction scripts are in [`validation/`](https://github.com/blex-max/expos-repo/tree/main/validation).

## Summary

- `TJAC` replaced `TRK` (Ripley's K over template endpoints) because `TRK`
  responded to only one thing — the presence of a near-coincident template pair
  — and **72.7% of real pileup columns contain one**, so it fired almost
  everywhere.
- **3.0 is roughly the 1% false-positive point for both `QRK` and `TJAC`**, not
  the textbook 2.0. Both statistics sum over pairs of observations, which leaves
  their nulls with a heavier right tail than a normal distribution.
- On 21,952 real variants both statistics show a strongly graded enrichment for
  a known-pathological genomic region and for poorly-aligning reads.

## Why `TRK` was replaced

`TRK` counted template pairs within a manhattan radius of 5 over their
endpoints. Over 995 pileup columns in `21:20.0–20.3Mb` (median 68 templates per
column), **72.7%** contained at least one such pair, averaging 1.36 per column —
yet only **0.0595%** of all template pairs qualified. A random draw of 8
templates therefore expects 0.017 such pairs, so a single one in a supporting set
is roughly 60× the null mean: enough to drive the statistic to significance on
its own.

Of the 1,354 qualifying pairs, exactly **one** sat at distance 0. Exact
duplicates essentially never survive duplicate marking, as expected; the pairs
that do survive sit at distances 1–4, where a difference in clipping or trimming
shifted the 5′ key enough for markdup to keep both copies.

Simulation on measured parameters (fragment lengths sampled from 99,263 real
templates, duplicate marking modelled, nuisance pairs offset by the measured
distance histogram), comparing each statistic at **its own** matched
false-positive rate:

| supporting-set scenario | `TJAC` | `TRK` |
| --- | --- | --- |
| markdup escapee present *(want low)* | **9.5%** | 59.5% |
| shared 50bp interval *(want high)* | **51.2%** | 2.5% |
| jittered recurrent break sites *(want high)* | 21.8% | **47.8%** |

`TRK` retains an advantage on recurrent break sites, but by the same mechanism as
its escapee sensitivity — so a `TRK` hit cannot be attributed to one cause or the
other.

!!! note "What `TJAC` responds to"
    `TJAC` measures how much of the supporting set is involved, not whether any
    single pair coincides. With all eight supporting templates near-identical it
    reaches a median effect size of 7.56 and fires 100% of the time; with five of
    eight, 3.43 and 94%; with three of eight, 1.49 and 35%. Summing over pairs
    deliberately dilutes a minority — the same property that stops a lone
    duplicate escapee (one pair in twenty-eight) from firing it.

## Threshold calibration

The null the code implements is precisely "the supporting set is a same-size
subsample of the background". That needs no generative model to reproduce: a null
replicate is *n* observations drawn from a real pileup column. `validation/null_calibration.py`
does this over four regions — two ordinary, two acrocentric — standardising each
column's draws by that column's own null moments, then pooling, since a single
published cutoff faces the mixture over loci.

The arithmetic behaves. For 151bp reads two random query positions fall within
the radius with probability 0.0587, predicting 1.64 pairs among 8 supporting
reads; real columns give 1.69. Variance runs slightly *below* the Poisson
expectation (var/mean 0.93–0.96), which is exactly the 1−p of a sum of Bernoulli
pair indicators. Acrocentric columns behave like ordinary ones.

The tail does not behave. The 1-in-100 point of the null effect size sits well
above 2.33 at every support size, for both statistics:

| supporting observations | `QRK` | `TJAC` | FPR at 2.7 | FPR at 3.0 |
| --- | --- | --- | --- | --- |
| 3 | 2.20 | 2.44 | 0.6% / 0.3% | 0.6% / 0.0% |
| 5 | 3.28 | 3.13 | 2.7% / 2.2% | 2.2% / 1.3% |
| 8 | 3.18 | 3.09 | 1.7% / 1.7% | 1.2% / 1.1% |
| 15 | 2.89 | 2.94 | 1.4% / 1.4% | 0.8% / 0.9% |
| 30 | 2.75 | 2.74 | 1.1% / 1.1% | 0.6% / 0.6% |
| 60 | 2.70 | 2.60 | 1.0% / 0.8% | 0.6% / 0.5% |

*False-positive rates are for the recommended conjunction (effect size above the
cutoff **and** p < 0.05), given as `QRK` / `TJAC`.*

**3.0 is the honest ~1% figure** — worst case 2.2% at five supporting reads,
0.5–1.3% elsewhere — and is what both headers now quote. The two statistics track
each other closely because the null is governed by the pair-sum structure rather
than by what is measured per pair, so one cutoff serves both.

!!! warning "Very low support"
    With three supporting reads there are only three pairs, so `QRK` can score
    0, 1, 2 or 3 and its effect size can only take the values −0.44, 2.01, 4.43
    or 6.86. Nothing lies between. A cutoff of 2.0 then admits 7.5% of pure noise
    while 2.1 admits 2.6% — a threefold swing from moving the threshold by a
    tenth. The p-value is what saves this case, capping the conjunction at 0.6%;
    from five reads upward the two conditions become largely redundant above
    ~2.3, so the p-value stops compensating for a low effect-size cutoff.

!!! note "Correcting an earlier figure"
    An earlier revision of this page gave `TJAC`'s 1% point as 2.72, from a
    simulation of fragments rather than resampling of real pileups. At the
    support size it was calibrated for the true rate is 1.7%. Real columns are
    more structured than simulated fragments; the design of that simulation
    should if anything have made *its* tail the heavier one, so the gap is
    understated rather than exaggerated.

## How many draws the null needs

`DEFAULT_NSIM` is 2,500. This section gives the evidence behind that, and — more
usefully — why no draw count is singled out by it.

**What the draw count controls.** For each variant the question asked is: *if these
supporting reads were merely a random handful of the reads covering this site, how
clustered would they look?* That is answered by drawing a random same-size handful
and measuring it, `nsim` times over. The pile of values that results is the null,
and the reported effect size says how far the real observation sits above the mean
of that pile, in units of its standard deviation.

So the draw count is a precision knob on the null and nothing else. More draws
sharpen the picture, fewer blur it; neither changes what the statistic measures.
That limits what "enough" can mean, and two candidate meanings are worth keeping
apart.

**Could a blurred null break the cutoff?** If a small `nsim` shifted the rate at
which pure noise clears 3.0, then the cutoff would no longer mean what the
calibration above says it means. That is a question about correctness and it has a
yes/no answer. Measured across six support sizes and millions of simulated loci:
**no.** The rate barely moves between 100 draws and an unlimited null.

**Could it change individual calls?** **Yes**, and this is the real cost. Annotate
the same VCF twice with different seeds and a few percent of flagged variants swap
sides.

Both are true at once, which is worth spelling out because it looks like a
contradiction. The error a small null puts on an effect size is random, not
directional: it pushes some loci above 3.0 and others below in roughly equal
numbers. The *number* of flagged loci therefore stays about right while *which*
loci they are shifts. Aggregate calibration survives; individual assignments churn.
The first question is about the count, the second about the identities.

Everything below quantifies the second, since the first is settled. The measurements
come from `validation/nsim_convergence.py` on the null regime and
`validation/nsim_churn.py` on real call sets.

### It does not buy the cutoff

The false-positive rate of the recommended conjunction — effect size above 3.0
**and** p < 0.05 — barely moves with the draw count. Over 0.96M to 5.8M null
replicates per cell, depending on statistic and support size, nothing sits more
than 0.3 points off the unlimited-null rate, and from 1,000 draws upwards nothing
is more than 0.12 off:

| supporting observations | 100 | 500 | 2,500 | unlimited |
| --- | --- | --- | --- | --- |
| 3 | 0.60 / 0.15 | 0.60 / 0.04 | 0.61 / 0.03 | 0.61 / 0.03 |
| 5 | 1.97 / 1.60 | 2.08 / 1.43 | 2.19 / 1.39 | 2.24 / 1.38 |
| 8 | 1.33 / 1.36 | 1.21 / 1.22 | 1.20 / 1.19 | 1.23 / 1.17 |
| 15 | 1.00 / 1.11 | 0.88 / 0.99 | 0.86 / 0.97 | 0.85 / 0.96 |
| 30 | 0.77 / 0.80 | 0.67 / 0.70 | 0.64 / 0.68 | 0.64 / 0.68 |
| 60 | 0.66 / 0.58 | 0.58 / 0.49 | 0.56 / 0.48 | 0.56 / 0.48 |

*Percentages, as `QRK` / `TJAC`. The calibration table above was measured at
20,000 draws, so the last column is the regime it describes.*

A small null smooths the decision across the threshold, and since the null's
density is falling there the smoothing nets slightly *more* false positives — in
every row of the table but one, `QRK` at five supporting observations, which goes
the other way. The largest drift, 0.27 points, is a fraction of the 0.6–2.2 point
spread the cutoff already carries across support sizes, so nothing about the
published 3.0 rests on the draw count. This was the check most likely to constrain
`nsim` from below, and it passed at every value tested.

Nor does the p-value constrain it. Its resolution is 1/(nsim+1), so a p < 0.05
screen remains reachable at all down to 20 draws, and the rates above show it
still doing its job at 100. (A p-value small enough to survive multiple-testing
correction across a whole call set is a different matter, and `nsim` is the wrong
lever for it — 20,000 variants would need a null of some 400,000 draws. The
p-value here is a companion screen, not a corrected test.)

### What it buys is precision in the effect size

Both statistics sum over pairs of observations, which leaves their nulls skewed
and heavy-tailed, and the standard error of a z-score taken against an *estimated*
null grows with the z being estimated:

```
se(z) = sqrt( (1 + z*skew + z^2 * (kurt - 1) / 4) / nsim )
```

for the null's skewness and kurtosis. Measured against the true effect size on
real pileup columns, this holds to within about 10% in every cell, so it can be
used to reason about draw counts that were never run. Note what the `z^2` term
means: precision is worst exactly where the cutoff sits. At the null's centre
`se(z)` is `1/sqrt(nsim)`; at z = 3 it is two to five times that, depending on how
skewed the null is at that support size.

The consequence is easier to read as the width of the band around 3.0 in which the
flag is between 5% and 95% likely — loci decided by the seed rather than by the
data:

| supporting observations | 100 | 500 | 1,000 | 2,500 | 5,000 |
| --- | --- | --- | --- | --- | --- |
| 5 | 1.44 / 1.19 | 0.65 / 0.52 | 0.47 / 0.38 | 0.31 / 0.24 | 0.22 / 0.16 |
| 8 | 1.23 / 1.22 | 0.55 / 0.54 | 0.39 / 0.39 | 0.26 / 0.24 | 0.19 / 0.17 |
| 15 | 1.09 / 1.15 | 0.48 / 0.51 | 0.34 / 0.36 | 0.22 / 0.24 | 0.16 / 0.17 |
| 30 | 0.98 / 1.00 | 0.44 / 0.44 | 0.30 / 0.32 | 0.20 / 0.21 | 0.15 / 0.14 |
| 60 | 0.95 / 0.92 | 0.41 / 0.42 | 0.30 / 0.29 | 0.19 / 0.19 | 0.13 / 0.14 |

*Width in units of effect size, as `QRK` / `TJAC`.*

At 2,500 draws a locus must sit within about ±0.1 of the cutoff before the seed
decides it. At 500 that window is ±0.25, and at 100 it is ±0.6 — wide enough to
swallow the whole interval between the textbook 2.0 and the calibrated 3.0.

There is no knee in any of this. The loci at risk are those within about one
`se(z)` of the cutoff, and their number falls as `1/sqrt(nsim)`, so every
quadrupling of the draw count halves the churn and no value is singled out.
Choosing `nsim` is a cost/precision trade with no threshold to discover, which is
why it is settled below on what the churn costs a reader rather than on any
property of the null.

!!! note "Low support is the safe end, for `QRK`"
    `QRK` counts pairs, so at three supporting observations it can only score 0,
    1, 2 or 3 and its effect size can only take four values roughly 2.4 apart.
    Nothing lands near 3.0, and churn at that support size is 1.7% at 100 draws
    and 0.0% at 2,500 — the lowest of any row. The lattice that makes the cutoff
    awkward at low support (see the warning above) also makes the draw count
    almost irrelevant there. `TJAC` sums a continuous quantity and has no such
    protection.

### What the churn costs on a real call set

The null regime above overstates churn, because under the null every flagged locus
is a marginal one sitting just above the cutoff. On real data much of the flagged
set is genuine artefact at an effect size of 5, 8 or 20, where no plausible draw
count changes the answer. `validation/nsim_churn.py` measures the operational
figure by annotating one VCF repeatedly with builds differing only in
`DEFAULT_NSIM`.

Over 8,000 germline `bcftools` calls, scored against a 50,000-draw reference run
and averaged over three seeds:

| draws | churn | mean \|Δz\| near 3.0 | 95th percentile |
| --- | --- | --- | --- |
| 2,500 | **3.3% / 2.8%** | 0.044 / 0.044 | 0.10 / 0.11 |
| 1,000 | 4.1% / 5.7% | 0.072 / 0.069 | 0.17 / 0.18 |
| 500 | 6.4% / 7.9% | 0.107 / 0.098 | 0.27 / 0.24 |
| 250 | 9.7% / 11.5% | 0.140 / 0.138 | 0.33 / 0.34 |
| 100 | 15.5% / 16.9% | 0.223 / 0.225 | 0.59 / 0.58 |

*`QRK` / `TJAC` throughout. Churn is a percentage of the 161 and 223 flags the
reference produces. `Δz` is measured against the reference over records whose
reference effect size lies in (2.5, 3.5). Runs predate the carried index
permutation in `subsample_wo_replace`, which reseeds the null without changing how
it is drawn, so these rates stand while the individual VCFs do not — see
[`validation/README.md`](https://github.com/blex-max/expos-repo/tree/main/validation).*

At the shipped default about 3% of flags differ — 5 records of 161 for `QRK`, 6 of
223 for `TJAC` — and **near the cutoff** the reported effect size is within about
0.1 of its unlimited-null value 95% of the time. Across the range the churn tracks
`1/sqrt(nsim)`: a 25-fold reduction in draws multiplies it by 4.7 and 6.0 against
the 5.0 exact scaling predicts, which is as close as bases of this size resolve.

The same measurement on 4,645 somatic CaVEMan calls — 99.4% of records scored
rather than 49% — gives seed-to-seed churn at 2,500 draws of 3.9% and 2.3%. Those
levels are **not comparable to the germline table**: the somatic flag rates are far
higher (19.3% of scored records for `QRK`, 44.5% for `TJAC`), and a larger flagged
set is proportionally more full of loci nowhere near the cutoff, which shrinks the
percentage mechanically. What is valid is the comparison *within* the somatic set,
across support bins on their own bases — and there the two statistics part company:

| supporting reads | `QRK` | `TJAC` |
| --- | --- | --- |
| < 8 | 3.2% *(278 flags)* | **4.9%** *(368 flags)* |
| 8–20 | 4.3% *(508)* | 1.8% *(1,533)* |
| 20–40 | 4.0% *(100)* | 1.4% *(147)* |

`QRK` is flat across support, because the lattice described above keeps its
low-support effect sizes far from 3.0. `TJAC` sums a continuous quantity and has no
such protection, so it is roughly three times more seed-sensitive below 8 supporting
reads — exactly what the null-regime harness predicts, where `TJAC` at three
supporting observations still churns 53% at 2,500 draws against `QRK`'s 0.0%.

This is the one place the draw count bites the regime `expos` targets, and it is
still modest in absolute terms: 4.9% is 18 flags of 368, and it is the *only* cell
in either call set where the shipped default exceeds 5%. At 1,000 draws it would be
nearer 8%. It is a reason to read a low-support `TJAC` flag as the least stable
thing the tool reports, not a reason to raise the default — the germline data could
not have shown it at all, having produced 0 and 1 flags in that bin.

!!! note "Reading these two measurements together"
    The null-regime harness runs millions of replicates and settles the *shape* of
    the curve; the end-to-end runs have a few hundred flags each and settle its
    *level*. Adjacent rows of the table above are not separable — a percentage
    point there is one or two records, and the three seeds share one reference run,
    so their effective precision is below three independent estimates. Where the
    two methods overlap they agree: seed-to-seed churn is 1.44× the churn against
    the reference, against the √2 = 1.41 expected if the density of loci near the
    cutoff is locally flat.

### The cost side

Cost is linear in the draw count. When these measurements were taken it was also
linear in *locus depth*, because `subsample_wo_replace` rebuilt an index vector the
size of the entire pileup on every draw however few observations it went on to
sample. That made the regime `expos` targets the expensive one — roughly 130ms per
somatic CaVEMan record against 20ms per germline `bcftools` record at 2,500 draws,
germline loci here being shallower with 56.8% of them skipping the Monte Carlo
entirely on `insufficient_background`. Differencing the wall clock at 2,500 against
1,000 draws put the Monte Carlo at about a third of germline runtime and the large
majority of somatic runtime.

Both figures are now **upper bounds rather than current cost**: that index
permutation is carried between draws instead of rebuilt, so the per-draw term no
longer scales with depth. They also came from the sweep runs rather than a dedicated
benchmark and were not taken on an otherwise-idle machine, so they were only ever
approximate. Re-measuring is the job of a benchmark, not of this page — the point
that survives is directional and unchanged: **the draw count costs more on somatic
call sets than a germline benchmark suggests**, because nearly every record is
scored and the loci are deeper.

That the allocation was worth removing is the general lesson. It bought runtime at
no cost in precision whatsoever, which is strictly better than any trade `nsim`
offers — so the per-draw path is the first place to look when runtime binds, and the
draw count the last.

### The choice

The most useful result here is a negative one: **there is no draw count to
discover.** Churn falls as `1/sqrt(nsim)` without a knee, a plateau or a stopping
point — four times the compute always halves it, at 100 draws and at 100,000 alike.
An "empirically optimal `nsim`" does not exist to be found. What the measurements
can do is price the options; choosing between them is a judgement about how much
churn is tolerable, and that judgement is stated below rather than dressed up as a
finding.

**2,500 is a defensible place to sit,** and the measurements support keeping it.
It puts about 3% of flags at odds with an unlimited null and holds the reported
effect size near the cutoff within about 0.1 of its limiting value 95% of the time,
with one cell — low-support `TJAC` — reaching 4.9%. The loci it churns sit within
~0.1 of a cutoff published as "~3.0" and meant to be read alongside three other
annotations — the regime where a single flag is least actionable anyway.

But the reason is an argument about consequences, not a statistical result, and it
should be read as one. The variants that churn are by construction those within
~0.1 of a cutoff this page publishes as "~3.0", on an annotation meant to be read
alongside three others rather than thresholded on alone. The instability lands
where a single flag is least load-bearing. That is a defensible reading; it is not
the only one.

The counter deserves stating plainly. If ±0.15 of slop on a cutoff specified to one
significant figure is acceptable — and the calibration table's own spread across
support sizes is wider than that — then **1,000 draws is also defensible**: 4–6%
churn for appreciably less runtime, and the draw count is not what limits this
tool's accuracy. Going the other way, 10,000 would halve the churn to ~1.5% for
four times the Monte-Carlo cost. Nothing here rules either out.

What *would* change the answer is a use that the "read the flags together" argument
does not cover: a downstream filter thresholding hard on `QRK` or `TJAC` alone, a
requirement that conclusions be reproducible across seeds and not merely across
reruns of a fixed one, or work concentrated in the low-support `TJAC` corner where
churn reaches 4.9%. Any of those makes 2,500 look thin and 10,000 the safer choice.
Absent them, 2,500 is where the seed stops visibly competing with the data, and that
is the whole of the case for it.

## Real-data behaviour

Full `expos` run over 21,952 germline `bcftools` calls. Both statistics were
computed for ~36% of records; 56.8% were skipped `insufficient_background`,
because the sample-size guard requires twice as many background as supporting
observations and a high-VAF germline variant is supported by nearly every read at
the locus. This is expected rather than a defect, but it makes germline call sets
a poor test bed — the somatic low-VAF case `expos` is built for is unaffected.
Median supporting depth among scored records was 28, so the `n = 30` row above is
the relevant null.

| | `QRK` observed | `TJAC` observed | null at *n* = 30 |
| --- | --- | --- | --- |
| effect size > 2.7 | 6.09% | 5.60% | 1.1% |
| effect size > 3.0 | 5.03% | 4.66% | 0.6% |

Binning `QRK` by effect size and cross-tabulating against quantities it does not
use:

| effect size | n | in `21:9–12Mb` | median `MLAS` | median `RCMPLX` |
| --- | --- | --- | --- | --- |
| < 0 | 3,591 | 10.2% | 0.964 | 1.928 |
| 0 – 1.95 | 3,483 | 10.1% | 0.934 | 1.925 |
| 1.95 – 3.0 | 453 | 14.8% | 0.914 | 1.915 |
| 3.0 – 4 | 178 | 21.9% | 0.868 | 1.882 |
| 4 – 8 | 168 | 25.6% | 0.825 | 1.854 |
| > 8 | 53 | **56.6%** | **0.808** | **1.684** |

Against a baseline of 11.3% of scored loci falling in `21:9–12Mb`, the share
rises monotonically to 56.6% in the extreme tail, and median `MLAS` — the
alignment-score statistic, computed from entirely different inputs — degrades
alongside it. The effect is not confined to that region: excluding it, 4.1% of
loci still exceed 3.0 against a 0.6% null.

`TJAC` behaves the same way on region membership (9.5% rising to 47.4%) and on
`MLAS` (0.954 falling to 0.834). Its top `RCMPLX` bin stays flat at ~1.91 where
`QRK`'s falls to 1.684, but that bin holds 38 loci and the appearance of no
association there does not survive a proper test.

## What the statistics track

Spearman correlations of effect size against three locus properties, over all
scored loci. `MLAS[1]` is used alongside `MLAS[0]` because `MLAS[0]` is computed
from the same supporting reads as the statistic and could correlate through
shared noise; `MLAS[1]` covers all reads at the locus. MAPQ comes straight from
the BAM and is fully independent — `expos` applies no mapping-quality filter, only
the flag mask.

| | `QRK` | `TJAC` |
| --- | --- | --- |
| `MLAS[0]`, supporting reads | −0.219 | −0.176 |
| `MLAS[1]`, all reads at locus | −0.162 | −0.136 |
| fraction of reads with MAPQ 0 | +0.109 | +0.149 |
| `RCMPLX` | −0.121 | **−0.091** |

Every one of these excludes zero at 95%. **`TJAC` does correlate with `RCMPLX`** —
the bin table simply lacked the power to show it. A paired bootstrap over the same
loci puts the `QRK`–`TJAC` gap at 0.033 in ρ for `RCMPLX` and 0.034 for `MLAS[1]`
(slightly stronger for `QRK` in both cases), and finds no distinguishable
difference on either MAPQ measure. The two statistics are far more alike in what
they track than the extreme bins suggest.

`MLAS[0]` correlates more strongly than `MLAS[1]` for both statistics, which is
expected rather than informative: it is computed from the same supporting reads,
so some of that gap is shared noise. `MLAS[1]` is the honest number.

The robust asymmetry is not between the statistics but between the covariates:
**both track alignment quality roughly twice as strongly as reference
complexity**, and the two are separable.

- Inside `21:9–12Mb`, `MLAS[0]` association strengthens sharply (−0.343 / −0.322)
  while `RCMPLX` association disappears entirely (−0.060, CI crossing zero, and
  +0.019).
- Restricted to loci of ordinary complexity (`RCMPLX` > 1.8), the MAPQ-zero
  association survives intact (+0.083 / +0.106), as does `MLAS[1]`
  (−0.114 / −0.098).

`RCMPLX` is a *local* property — Lempel-Ziv entropy of a 100bp window. Mismapping
is driven by homology *elsewhere in the genome*, and a segmental duplication has
perfectly ordinary local entropy: it is unremarkable sequence that happens to
exist in two places. Local entropy is blind to that by construction, while MAPQ
and alignment score are not. So there are two distinct routes to a spurious call,
and `RCMPLX` measures the one the others are worst at seeing.

That separation is a reason to keep `RCMPLX`, not a limitation of it. Its weak
correlation with the clustering statistics is exactly what an independent axis
should look like: of 657 loci flagged by `QRK` or `TJAC` above 3.0, only 15.7%
have `RCMPLX` < 1.5, and of the low-complexity loci only 33% are flagged by
either statistic. The sets barely overlap. The loci `RCMPLX` flags alone are
nonetheless real trouble — among records neither clustering statistic flags,
those with `RCMPLX` < 1.5 have a median `MLAS[1]` of 0.868 against 0.967 for the
rest. `RCMPLX` is also the only annotation `expos` produces that is derived from
the reference rather than the reads, so it cannot agree with the others
circularly, and it is computed for 92.5% of records against 36.7% for the
clustering statistics, which the sample-size guard gates and it does not.

## Why the acrocentric enrichment is meaningful

`21:9–12Mb` in GRCh37 covers the pericentromeric region and proximal short arm of
chromosome 21. The short arms of the five acrocentric chromosomes are built
almost entirely from repeat families shared *between* them — ribosomal DNA
arrays, alpha and beta satellite, large segmental duplications — and in GRCh37
these are largely unassembled or represented by short placeholder models. Reads
originating there have no correct home and are placed into whatever
representative sequence exists, which is why these regions appear in every
standard exclusion set (the ENCODE blacklist, the GIAB low-mappability and
segmental-duplication strata).

That makes enrichment there useful evidence for three reasons. It is
**external**: region membership is a fixed property of the reference assembly,
not derived from the read data or from any `expos` statistic, so unlike agreement
between two statistics computed from the same pileup it cannot be circular. The
**prior is independent**: that calls here are predominantly artefactual is
established knowledge unconnected to this tool. And the relationship is
**graded** rather than merely present — it rises across every bin, which is much
harder to produce by chance than a single tail association.

There is also a mechanism specific to these statistics. Mismapped reads do not
arrive from arbitrary places; they originate from a limited number of true source
loci. Reads from a common source forced into a common destination arrive with
correlated geometry — stacked, with shared boundaries and shared query offsets.
That is precisely what `QRK` and `TJAC` measure. The enrichment is not merely
"fires where calling is difficult" but "fires where the property being quantified
is expected to arise."

## Convergent evidence

The annotations are designed to be read together: each covers a different route
to a spurious call, so a variant that trips several is far more suspect than one
that trips any single test. Counting red flags per locus — `QRK` > 3.0 with
p < 0.05, `TJAC` > 3.0 with p < 0.05, `MLAS[0]` < 0.93, `RCMPLX` < 1.5 — over the
7,757 records carrying all four, and scoring the result against region membership
and MAPQ, neither of which any of the four uses:

| red flags | n | % of set | in `21:9–12Mb` | any MAPQ-0 read |
| --- | --- | --- | --- | --- |
| 0 | 4,970 | 64.1% | 7.5% | 11.4% |
| 1 | 2,161 | 27.9% | 14.3% | 33.5% |
| 2 | 456 | 5.9% | 28.3% | 62.5% |
| 3 | 154 | 2.0% | **39.0%** | **79.2%** |
| 4 | 16 | 0.2% | 37.5% | 62.5% |

Against an 11.3% baseline, three flags puts a locus in the difficult region 39% of
the time — five times the rate for clean loci — and 79.2% of such loci carry at
least one ambiguously-placed read, against 11.4% of unflagged loci.

`bwa` assigns MAPQ 0 when a read has two or more equally-good best alignments. The
read is not discarded: it is placed arbitrarily among the tied candidates and
enters the pileup like any other, so the support it lends a variant may belong to
the paralogous locus instead. That is the mismapping mechanism described above,
observed directly rather than inferred — and reads arriving together from a common
true source are exactly the ones that carry the correlated geometry `QRK` and
`TJAC` measure. In this sample MAPQ 0 accounts for 0.4% of reads in an ordinary
region against 27.5% at `21:9.5–10Mb`.

Two caveats. The four-flag bin holds 16 loci and says nothing on its own; the
apparent plateau there is noise, not a ceiling. And requiring all four is so
restrictive as to be useless operationally — 0.2% of records. The flag *count* is
the usable quantity, not the full conjunction.

## Limitations

- One sample, one region pair, germline call set. `expos` targets somatic low-VAF
  calling, which was not tested here — and that is the regime where support
  counts are lowest and the cutoff least stable.
- `QRK` and `TJAC` correlate with `MLAS` and with each other. That correlation is
  what corroborates them, but they are not independent evidence.
- The end-to-end draw-count figures rest on a few hundred flagged records, so
  neighbouring values of `nsim` are not separable there; only the null-regime
  harness has the replicate count to resolve the curve, and it measures pure noise
  rather than real artefacts. The two agree on shape, which is the basis for
  quoting either.
- Enrichment in a difficult region is consistent with the mismapping mechanism
  above but does not prove it. Without truth data these calls cannot be separated
  into "artefact" and "real variant in a hard region", though for quality-control
  purposes both warrant flagging.
- The calibration samples columns ~100bp apart, so with 151bp reads adjacent
  columns share reads and there are fewer independent observations than the
  column counts suggest. Moments and percentiles are also taken from the same
  draws, whereas the tool estimates its null from 2,500 independent ones. Both
  push the estimate conservative.
- The recurrent-break-site scenario in the `TRK` comparison uses invented
  geometry (three break points, ±3bp). Real fragmentation discreteness was not
  measured, and that is the one scenario where `TRK` outperformed `TJAC`, so it
  is the least trustworthy number here.
