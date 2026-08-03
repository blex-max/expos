#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <optional>
#include <random>
#include <string_view>
#include <vector>

// Pure statistical primitives for the expos analysis engine. Depends only
// on the standard library so it can be unit-tested in isolation. Guards
// (insufficient observations/background) belong to the calling statistic
// layer, not here — these functions assume valid inputs.

// --- reference complexity --- //

// Lempel-Ziv 76 entropy rate (bits per char). Used as a reference
// sequence-complexity measure.
double entropy_lz76 (std::string_view s);

// --- summary statistics --- //

// Arithmetic mean; nullopt for an empty sample.
std::optional<double> mean (const std::vector<double>& v);

// Linear-interpolated percentile at fraction pt (0 < pt < 1); nullopt for
// an empty sample.
template <class T>
  requires std::unsigned_integral<T> || std::floating_point<T>
std::optional<double> percentile (std::vector<T> obs, double pt)
{
  assert (pt > 0 && pt < 1);
  if (obs.empty()) {
    return std::nullopt;
  }
  if (obs.size() == 1) {
    return static_cast<double> (obs[0]);
  }
  std::sort (obs.begin(), obs.end());
  const double rank = static_cast<double> (obs.size() - 1) * pt;
  const double lower = std::floor (rank);
  const double frac = rank - lower;
  const auto lowerI = static_cast<std::size_t> (lower);
  const auto upperI =
      static_cast<std::size_t> (std::ceil (rank));
  if (lowerI == upperI) {
    return static_cast<double> (obs[lowerI]);
  }
  return static_cast<double> (obs[lowerI]) +
         (frac * (static_cast<double> (obs[upperI]) -
                  static_cast<double> (obs[lowerI])));
}

// --- spatial clustering (Ripley's K) --- //

// On "pairs within radius" vs the per-point definition: Ripley's K sums,
// over every point, that point's neighbour count within the radius. That
// sum counts each pair twice (once from each end), so it equals
// 2 * (unordered pairs within radius). The count_* functions below return
// that unordered pair count directly -- there is deliberately no per-point
// loop; ripley_k restores the factor of 2 and the per-point / n.

// Count unordered pairs whose 1-D separation is strictly < radius.
// Sort + two-pointer sweep, O(N log N). This is the hot path for the
// query-position statistic (values are read query positions: int32_t).
std::size_t count_pairs_within_1d (
    std::vector<int32_t> pts, uint64_t radius
);

// Count unordered pairs whose dist(a, b) is strictly < radius, O(N^2).
// For metrics that don't reduce to a 1-D sort, e.g. manhattan distance
// over template-endpoint pairs. dist must return a value comparable to
// radius.
template <class T, class DistFn>
std::size_t count_pairs_within (
    const std::vector<T>& pts, uint64_t radius, DistFn dist
)
{
  std::size_t count = 0;
  const std::size_t n = pts.size();
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      if (dist (pts[i], pts[j]) < radius) {
        ++count;
      }
    }
  }
  return count;
}

// Ripley's K (no edge correction) from an unordered within-radius pair
// count over n points at the given point intensity. Equal to the mean
// per-point neighbour count scaled by 1/intensity. Precondition: n >= 1.
double ripley_k (
    std::size_t pairsWithin, std::size_t n, double intensity
);

// --- monte carlo --- //

// Default number of simulation draws for the null distribution.
inline constexpr std::size_t DEFAULT_NSIM = 2500;

struct McResult {
  // z-score of the observed statistic against the simulated null;
  // nullopt when the null has zero variance.
  std::optional<double> effectSize;
  // one-sided p-value, (#draws >= observed + 1) / (nsim + 1).
  double pValue = 1.0;
};

// Compare `observed` against a null built by repeatedly drawing a
// background sample (`draw`) and reducing it to a scalar (`stat`).
// Single pass, Welford accumulation — no per-draw storage. The p-value is
// one-sided and Laplace-smoothed; the effect size is a z-score against the
// null's population standard deviation.
template <class DrawFn, class StatFn>
McResult monte_carlo_pvalue (
    double observed, DrawFn draw, StatFn stat,
    std::size_t nsim = DEFAULT_NSIM
)
{
  assert (nsim > 0);
  std::size_t countGe = 0;
  double meanAcc = 0.0;
  double m2 = 0.0;
  for (std::size_t k = 1; k <= nsim; ++k) {
    const double s = static_cast<double> (stat (draw()));
    if (s >= observed) {
      ++countGe;
    }
    // Welford online mean/variance
    const double delta = s - meanAcc;
    meanAcc += delta / static_cast<double> (k);
    m2 += delta * (s - meanAcc);
  }

  McResult res;
  res.pValue = static_cast<double> (countGe + 1) /
               static_cast<double> (nsim + 1);
  const double popSd =
      std::sqrt (m2 / static_cast<double> (nsim));
  if (popSd > 0.0) {  // otherwise z-score is undefined
    res.effectSize = (observed - meanAcc) / popSd;
  }
  return res;
}

// Uniform sample of size n drawn without replacement from obs, via a
// partial Fisher-Yates shuffle of an index vector. Precondition: n <= obs
// size.
template <class T>
std::vector<T> subsample_wo_replace (
    const std::vector<T>& obs, std::size_t n, std::mt19937& rng
)
{
  const std::size_t nObs = obs.size();
  assert (n <= nObs);
  std::vector<std::size_t> idx (nObs);
  std::iota (idx.begin(), idx.end(), std::size_t{0});
  std::vector<T> out;
  out.reserve (n);
  for (std::size_t i = 0; i < n; ++i) {
    std::uniform_int_distribution<std::size_t> dist (
        i, nObs - 1
    );
    const std::size_t j = dist (rng);
    std::swap (idx[i], idx[j]);
    out.push_back (obs[idx[i]]);
  }
  return out;
}
