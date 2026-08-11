#include "compute_info_field.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

#include "expos/guards.hpp"
#include "expos/pileup_features.hpp"
#include "expos/skip.hpp"
#include "shared/stats.hpp"

// --- helpers --- //

// Builders keeping the compute_* functions readable.
static StatValue stat_value (double v)
{
  return {v, std::nullopt};
}

static StatValue stat_missing (StatSkipReason reason)
{
  return {std::nullopt, reason};
}

static StatValue stat_or (
    std::optional<double> v, StatSkipReason reasonIfMissing
)
{
  return v ? stat_value (*v) : stat_missing (reasonIfMissing);
}

// --- kernels --- //

// Both Monte-Carlo statistics are a sum, over unordered pairs of the
// sampled observations, of a symmetric function of just those two
// observations. Nothing else about the sample enters. QRK's kernel is
// the indicator |qPos_i - qPos_j| < radius; TJAC's is the Jaccard
// overlap of two templates.

// --- query position spatial clustering (Ripley's K) --- //

// Ripley's K in its general form sums, over every point, that point's
// neighbour count within the radius. Here the unordered pair count is returned
// directly, because the constants (intensity normalisation, etc.)
// all cancel when the count is compared against a same-size Monte-Carlo null.

// Count unordered pairs whose 1-D separation is strictly < radius.
// Sorts pts in place, so callers holding data that must survive the call
// pass a copy. Taking a mutable view rather than a value is what lets the
// Monte-Carlo draws reuse one buffer instead of copying per draw.
static std::size_t count_pairs_within_1d (
    std::span<int32_t> pts, uint64_t radius
)
{
  std::sort (pts.begin(), pts.end());
  std::size_t count = 0;
  std::size_t l = 0;
  for (std::size_t r = 0; r < pts.size(); ++r) {
    // advance the left edge until pts[r] - pts[l] < radius. Bounded by r
    // so radius == 0 yields no pairs rather than underflowing. pts are
    // sorted and non-negative, so the diff is a non-negative uint64_t.
    while (l < r &&
           static_cast<uint64_t> (pts[r] - pts[l]) >= radius) {
      ++l;
    }
    count += r - l;
  }
  return count;
}

// --- template overlap --- //

static size_t overlap (
    const TemplateEndpoints& a, const TemplateEndpoints& b
)
{
  const auto overlapLEdge = std::max (a.first, b.first);
  const auto overlapREdge = std::min (a.second, b.second);
  return (overlapREdge > overlapLEdge)
             ? static_cast<size_t> (overlapREdge - overlapLEdge)
             : 0;
}

static size_t size (const TemplateEndpoints& t)
{
  return static_cast<size_t> (t.second - t.first);
}

// Fraction of the union of two templates that they share. Unlike a
// min(len) denominator this does not saturate when one template is nested
// in a longer one, so the null keeps headroom for real coincidence to show
// against. Reaches 1.0 only for identical intervals.
static double jaccard (
    const TemplateEndpoints& a, const TemplateEndpoints& b
)
{
  const auto olap = static_cast<double> (overlap (a, b));
  const auto unionSz =
      (static_cast<double> (size (a) + size (b))) - olap;
  return olap / unionSz;
}

static double pairwise_jaccard_sum (
    std::span<const TemplateEndpoints> obs
)
{
  // NOTE: it would be possible to chase further
  // perf gains by `tiling`. Jaccard-ing several j
  // for each i, i.e. i + 4 <= n; i += 4. But
  // it would be more for fun than an important
  // major gain at this point.
  double jaccardSum = 0.0;
  // attempting to collapse the nested loop
  // with omp simd collapse(2) does not help
  for (size_t i = 0; i < obs.size(); ++i) {
// but reduction does!
// request vectorisation reduction via
// openmp-simd. Works for clang and gcc
#pragma omp simd reduction(+ : jaccardSum)
    for (size_t j = i + 1; j < obs.size(); ++j) {
      jaccardSum += jaccard (obs[i], obs[j]);
    }
  }
  return jaccardSum;
}

// --- exact null moments --- //

// Because each statistic is a pair-sum, the null's mean and variance
// under a uniform size-n draw without replacement are available in
// closed form from three sums over the background pool. No
// simulation needed for the effect size. This is SRSWOR
// for a pair-sum.
struct KernelTerms {
  double t1;  // kernel summed over every background pair
  double t2;  // kernel value squaring each pair's value first
  double q;  // sum of the squared per-observation row sums
};

struct NullMoments {
  double mean;
  double var;
};

// The three kernel terms are weighted by the probability that a given pair,
// triple or quadruple of background observations all land in one draw
// -- the three ways two pairs can overlap: identical, sharing one
// observation, disjoint.
static NullMoments null_moments (
    size_t nObs, const KernelTerms& k, size_t nBackground
)
{
  const auto n = static_cast<double> (nObs);
  const auto bgN = static_cast<double> (nBackground);
  // No division by zero: size_guard puts nBackground at MIN_BACKGROUND
  // or above. A support below four zeroes p3/p4 through their own
  // numerators, which is right -- a draw of three cannot contain four
  // distinct observations.
  const double p2 = (n * (n - 1)) / (bgN * (bgN - 1));
  const double p3 = p2 * ((n - 2) / (bgN - 2));
  const double p4 = p3 * ((n - 3) / (bgN - 3));

  const double sqT1 = k.t1 * k.t1;
  return {
      p2 * k.t1, (p2 * k.t2) + (p3 * (k.q - (2 * k.t2))) +
                     (p4 * (sqT1 + k.t2 - k.q)) -
                     ((p2 * p2) * sqT1)
  };
}

// Sorts pts in place. The kernel is an indicator, so every pair value
// is 0 or 1 and squaring changes nothing: t2 == t1.
static KernelTerms qpos_null_terms (
    std::span<int32_t> pts, uint64_t radius
)
{
  std::sort (pts.begin(), pts.end());
  const std::size_t n = pts.size();
  std::size_t pairCount = 0;
  std::size_t rowSqSum = 0;
  std::size_t l = 0;
  std::size_t r = 0;
  for (std::size_t c = 0; c < n; ++c) {
    // r is monotone in (increases with) c, so it is carried rather than reset; the
    // clamp keeps the window containing c when the previous point had
    // no right neighbour.
    r = std::max (r, c);
    while (r + 1 < n &&
           static_cast<uint32_t> (pts[r + 1] - pts[c]) <
               radius) {
      ++r;
    }
    while (l < c &&
           static_cast<uint32_t> (pts[c] - pts[l]) >= radius) {
      ++l;
    }
    rowSqSum += (r - l) * (r - l);
    // left-only, so each pair is counted once
    pairCount += c - l;
  }
  const auto t1 = static_cast<double> (pairCount);
  return {t1, t1, static_cast<double> (rowSqSum)};
}

static KernelTerms jaccard_null_terms (
    std::span<const TemplateEndpoints> bg
)
{
  const std::size_t n = bg.size();
  double t1 = 0.0;
  double t2 = 0.0;
  std::vector<double> rowSum (n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    double iAcc = 0.0;
#pragma omp simd reduction(+ : t1, t2, iAcc)
    for (std::size_t j = i + 1; j < n; ++j) {
      const double f = jaccard (bg[i], bg[j]);
      t1 += f;
      t2 += f * f;
      iAcc += f;
      rowSum[j] += f;
    }
    rowSum[i] += iAcc;
  }
  double q = 0.0;
#pragma omp simd reduction(+ : q)
  for (std::size_t i = 0; i < n; ++i) {
    q += rowSum[i] * rowSum[i];
  }
  return {t1, t2, q};
}

// precondition: nm.var > 0
static double z_score (double obsStat, const NullMoments& nm)
{
  return (obsStat - nm.mean) / sqrt (nm.var);
}

// --- statistics --- //

void set_qrk_guard (
    StatContext& sCtx, const PileupFeatures& primaryAll
)
{
  if (!sufficient_reads (primaryAll.readLen.size())) {
    // fail-closed
    sCtx.readLenSuppression =
        StatSkipReason::readLengthUnverified;
    return;
  }
  sCtx.readLenSuppression =
      read_lens_within_tol (primaryAll.readLen)
          ? std::nullopt
          : std::optional{
                StatSkipReason::heterogeneousReadLength
            };
}

static ValuesOrSkip compute_qrk (
    const VariantStatInputs& in, const StatContext& ctx
)
{
  constexpr uint64_t QPOS_RADIUS = 5;
  const auto& obs = in.supporting.qPos;
  const auto& bg = in.all.qPos;

  // Size before quality: a locus with too few reads must report that it
  // was too thin, not a verdict on read-length spread that the guard
  // never had the data to reach.
  if (const auto reason = size_guard (obs.size(), bg.size())) {
    return std::unexpected (*reason);
  }

  // Query position is an offset within the read, so a heterogeneous read
  // population would confound
  if (ctx.readLenSuppression) {
    return std::unexpected (*ctx.readLenSuppression);
  }

  // Exact null moments.
  auto bgCopy = bg;  // qpos_null_terms sorts, hence copy.
  const auto mom = null_moments (
      obs.size(), qpos_null_terms (bgCopy, QPOS_RADIUS),
      bg.size()
  );

  if (!(mom.var > 0)) {
    return std::unexpected (StatSkipReason::zeroVariance);
  }

  // count_pairs_within_1d sorts, so the observed statistic works on a copy
  // and the record's features stay const. One copy per record, not per draw.
  auto obsCopy = obs;
  const double obsStat = static_cast<double> (
      count_pairs_within_1d (obsCopy, QPOS_RADIUS)
  );

  const auto zScore = z_score (obsStat, mom);

  auto draw = [&] {
    return subsample_wo_replace (
        bg, obs.size(), ctx.mc.rng, ctx.mc.qPos
    );
  };
  auto stat = [&] (std::span<int32_t> sample) {
    return static_cast<double> (
        count_pairs_within_1d (sample, QPOS_RADIUS)
    );
  };
  const auto pVal = run_monte_carlo (obsStat, draw, stat);
  return std::vector<StatValue>{
      stat_value (zScore), stat_value (pVal)
  };
}

static ValuesOrSkip compute_tjac (
    const VariantStatInputs& in, const StatContext& ctx
)
{
  const auto& obs = in.supporting.endpoints;
  const auto& bg = in.all.endpoints;

  if (const auto reason = size_guard (obs.size(), bg.size())) {
    return std::unexpected (*reason);
  }

  // Template endpoints track fragment boundaries, not read length,
  // so not gated on read-length homogeneity.

  // Exact null moments.
  const auto mom = null_moments (
      obs.size(), jaccard_null_terms (bg), bg.size()
  );

  if (!(mom.var > 0)) {
    return std::unexpected (StatSkipReason::zeroVariance);
  }

  const auto obsStat = pairwise_jaccard_sum (obs);

  const auto zScore = z_score (obsStat, mom);

  const auto draw = [&] {
    return subsample_wo_replace (
        bg, obs.size(), ctx.mc.rng, ctx.mc.endpoints
    );
  };

  const auto pVal =
      run_monte_carlo (obsStat, draw, pairwise_jaccard_sum);
  return std::vector<StatValue>{
      stat_value (zScore), stat_value (pVal)
  };
}

// The only statistic with independently-missing subfields: two unrelated
// medians, either of which can be absent while the other is reportable.
static ValuesOrSkip compute_mlas (
    const VariantStatInputs& in, const StatContext&
)
{
  return std::vector<StatValue>{
      stat_or (
          percentile (
              in.supporting.normalisedAs, PERCENTILE_MEDIAN
          ),
          StatSkipReason::noSupport
      ),
      stat_or (
          percentile (in.all.normalisedAs, PERCENTILE_MEDIAN),
          StatSkipReason::noBackground
      )
  };
}

static ValuesOrSkip compute_rcmplx (
    const VariantStatInputs& in, const StatContext&
)
{
  constexpr size_t WIN_SZ = 100;
  constexpr size_t WIN_STEP = 10;

  const std::string_view ref = in.refSlice;
  // too-short spans have no full window; masked/ambiguous bases make the
  // complexity meaningless.
  if (ref.size() < WIN_SZ) {
    return std::unexpected (StatSkipReason::referenceTooShort);
  }
  if (ref.find_first_of ("Nn") != std::string_view::npos) {
    return std::unexpected (StatSkipReason::referenceHasN);
  }

  double entropySum = 0.0;
  size_t nWin = 0;
  for (size_t winStart = 0; winStart + WIN_SZ <= ref.size();
       winStart += WIN_STEP) {
    entropySum += entropy_lz76 (ref.substr (winStart, WIN_SZ));
    ++nWin;
  }
  const double meanWindowEntropy =
      entropySum / static_cast<double> (nWin);
  return std::vector<StatValue>{stat_value (meanWindowEntropy)};
}

// --- INFO header definitions --- //

constexpr std::string_view QRK_HEADER =
    "##INFO=<ID=QRK,Number=2,Type=Float,Description=\"Spatial "
    "clustering of mutant query positions against all reads: "
    "[0] standardised effect size (z-score) against the exact "
    "null; [1] one-sided Monte-Carlo p-value. Effect sizes "
    "greater than ~3.0 with a significant p-value may indicate "
    "a spurious variant.\">";
constexpr std::string_view TJAC_HEADER =
    "##INFO=<ID=TJAC,Number=2,Type=Float,Description=\"Graded "
    "pairwise overlap of mutant templates against all reads: "
    "[0] standardised effect size (z-score) against the exact "
    "null; [1] one-sided Monte-Carlo p-value. Effect sizes "
    "greater than ~3.0 with a significant p-value may indicate "
    "a spurious variant.\">";
constexpr std::string_view MLAS_HEADER =
    "##INFO=<ID=MLAS,Number=2,Type=Float,Description=\"Median "
    "read-length-normalised alignment scores: [0] of reads "
    "supporting the "
    "variant; [1] of all reads covering the variant site.\">";
constexpr std::string_view RCMPLX_HEADER =
    "##INFO=<ID=RCMPLX,Number=1,Type=Float,Description=\"Mean "
    "100-base "
    "window complexity (Lempel-Ziv 76 entropy rate) of the "
    "reference "
    "within 400 bases either side of the REF allele.\">";

// --- registry --- //
constexpr std::array<VariantStat, 4> VARIANT_STATS = {{
    {{"QRK", QRK_HEADER, 2}, &compute_qrk},
    {{"TJAC", TJAC_HEADER, 2}, &compute_tjac},
    {{"MLAS", MLAS_HEADER, 2}, &compute_mlas},
    {{"RCMPLX", RCMPLX_HEADER, 1}, &compute_rcmplx},
}};

std::span<const VariantStat> expos_field_registry()
{
  return VARIANT_STATS;
}
