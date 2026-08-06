#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <optional>
#include <random>
#include <utility>
#include <vector>

#include "expos/pileup_features.hpp"

// --- query position spatial clustering (Ripley's K) --- //

// Ripley's K in its general sums, over every point, that point's
// neighbour count within the radius. Here the unordered pair count is returned
// directly, because the constants (intensity normalisation, etc.)
// all cancel when the count is compared against a same-size Monte-Carlo null.

// Count unordered pairs whose 1-D separation is strictly < radius.
std::size_t count_pairs_within_1d (
    std::vector<int32_t> pts, uint64_t radius
);

// --- template overlap --- //

size_t overlap (
    const TemplateEndpoints& a, const TemplateEndpoints& b
);

size_t size (const TemplateEndpoints& t);

// Fraction of the union of two templates that they share. Unlike a
// min(len) denominator this does not saturate when one template is nested
// in a longer one, so the null keeps headroom for real coincidence to show
// against. Reaches 1.0 only for identical intervals.
double jaccard (
    const TemplateEndpoints& a, const TemplateEndpoints& b
);

double pairwise_jaccard_sum (
    const std::vector<TemplateEndpoints>& obs
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

// one-sided pval; the effect size is a z-score against the
// null's population standard deviation.
template <typename DrawFn, typename StatFn>
McResult run_monte_carlo (
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
// partial Fisher-Yates shuffle of an index vector.
// Precondition: n <= obs.size.
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
