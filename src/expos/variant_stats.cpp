#include "variant_stats.hpp"

#include <htslib/vcf.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>

#include "expos/pileup_features.hpp"

// --- stats --- //

std::size_t count_pairs_within_1d (
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

size_t overlap (
    const TemplateEndpoints& a, const TemplateEndpoints& b
)
{
  const auto overlapLEdge = std::max (a.first, b.first);
  const auto overlapREdge = std::min (a.second, b.second);
  return (overlapREdge > overlapLEdge)
             ? static_cast<size_t> (overlapREdge - overlapLEdge)
             : 0;
}

size_t size (const TemplateEndpoints& t)
{
  return static_cast<size_t> (t.second - t.first);
}

// Fraction of the union of two templates that they share. Unlike a
// min(len) denominator this does not saturate when one template is nested
// in a longer one, so the null keeps headroom for real coincidence to show
// against. Reaches 1.0 only for identical intervals.
double jaccard (
    const TemplateEndpoints& a, const TemplateEndpoints& b
)
{
  const auto olap = static_cast<double> (overlap (a, b));
  const auto unionSz =
      (static_cast<double> (size (a) + size (b))) - olap;
  return olap / unionSz;
}

double pairwise_jaccard_sum (
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
