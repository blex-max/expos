#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <optional>
#include <string_view>
#include <vector>

// --- reference complexity --- //

// Lempel-Ziv 76 entropy rate (bits per char). Used as a reference
// sequence-complexity measure.
double entropy_lz76 (std::string_view s);

// --- summary statistics --- //

// Arithmetic mean; nullopt for an empty sample.
std::optional<double> mean (const std::vector<double>& v);

// Linear-interpolated percentile at fraction pt (0 < pt < 1); nullopt for
// an empty sample.
static constexpr double PERCENTILE_MEDIAN = 0.5;
template <class T>
  requires std::integral<T> || std::floating_point<T>
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

template <std::ranges::input_range R>
  requires std::integral<std::ranges::range_value_t<R>> ||
           std::floating_point<std::ranges::range_value_t<R>>
std::optional<double> percentile (R&& r, double pt)
{
  using T = std::ranges::range_value_t<R>;
  return percentile (
      std::vector<T> (
          std::ranges::begin (r), std::ranges::end (r)
      ),
      pt
  );
}
