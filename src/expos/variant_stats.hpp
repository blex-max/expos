#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <optional>
#include <random>
#include <span>
#include <utility>
#include <vector>

#include "expos/pileup_features.hpp"

// --- query position spatial clustering (Ripley's K) --- //

// Ripley's K in its general form sums, over every point, that point's
// neighbour count within the radius. Here the unordered pair count is returned
// directly, because the constants (intensity normalisation, etc.)
// all cancel when the count is compared against a same-size Monte-Carlo null.

// Count unordered pairs whose 1-D separation is strictly < radius.
// Sorts pts in place, so callers holding data that must survive the call
// pass a copy. Taking a mutable view rather than a value is what lets the
// Monte-Carlo draws reuse one buffer instead of copying per draw.
std::size_t count_pairs_within_1d (
    std::span<int32_t> pts, uint64_t radius
);

// --- template overlap --- //

size_t overlap (
    const TemplateEndpoints& a, const TemplateEndpoints& b
);

size_t size (const TemplateEndpoints& t);

double jaccard (
    const TemplateEndpoints& a, const TemplateEndpoints& b
);

double pairwise_jaccard_sum (
    std::span<const TemplateEndpoints> obs
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

// Reusable buffers for repeated draws from one population.
template <class T>
struct SubsampleScratch {
  std::vector<std::size_t> idxBuf;
  std::vector<T> subsampleBuf;
};

// Per-run mutable state for the Monte-Carlo statistics: the shared RNG stream
// and one scratch per statistic. Built once for the whole record loop, so the
// buffers reach a high-water mark and stop reallocating. A new Monte-Carlo
// statistic adds a scratch here alongside its compute_ function and registry
// row.
struct McState {
  std::mt19937 rng;
  SubsampleScratch<int32_t> qPos;
  SubsampleScratch<TemplateEndpoints> endpoints;
};

// Uniform sample of size n drawn without replacement from obs, via a partial
// Fisher-Yates shuffle of scratch's persistent index permutation. The returned
// span views scratch.sampleBuf and is valid until the next draw from the same
// Precondition: n <= obs.size.
template <class T>
std::span<T> subsample_wo_replace (
    const std::vector<T>& obs, std::size_t n, std::mt19937& rng,
    SubsampleScratch<T>& scratch
)
{
  const std::size_t nObs = obs.size();
  assert (n <= nObs);
  if (scratch.idxBuf.size() != nObs) {
    scratch.idxBuf.resize (nObs);
    std::iota (
        scratch.idxBuf.begin(), scratch.idxBuf.end(),
        std::size_t{0}
    );
  }
  if (scratch.subsampleBuf.size() < n) {
    scratch.subsampleBuf.resize (n);
  }
  for (std::size_t i = 0; i < n; ++i) {
    std::uniform_int_distribution<std::size_t> dist (
        i, nObs - 1
    );
    const std::size_t j = dist (rng);
    std::swap (scratch.idxBuf[i], scratch.idxBuf[j]);
    scratch.subsampleBuf[i] = obs[scratch.idxBuf[i]];
  }
  // Only the first n entries were written this call
  // (any further entries stale)
  return std::span<T>{scratch.subsampleBuf.data(), n};
}
