#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <optional>
#include <vector>

namespace summary {

constexpr inline std::optional<double> mean (
    const std::vector<double>& v
)
{
  if (v.empty()) {
    return std::nullopt;
  }

  long double sum = 0.0L;
  for (const auto& x : v) {
    sum += static_cast<long double> (x);
  }

  return static_cast<double> (sum) /
         static_cast<double> (v.size());
}

template <typename T>
  requires std::unsigned_integral<T> || std::floating_point<T>
constexpr inline std::optional<double> percentile (
    std::vector<T> obs, double pt
)
{
  assert (pt > 0 && pt < 1);

  if (obs.empty()) {
    return std::nullopt;
  }

  if (obs.size() == 1) {
    return obs[0];
  }

  std::sort (begin (obs), end (obs));

  double pi = static_cast<double> (obs.size() - 1) *
              pt;     // 0 indexed rank
  auto lower = floor (pi);
  auto frac = pi - lower;
  auto upper_i = static_cast<size_t> (ceil (pi));
  auto lower_i = static_cast<size_t> (lower);

  if (lower_i == upper_i) {
    return obs[lower_i];
  }

    // linear interpolation
  // (NOTE: in principle, precision loss. Unlikely in this domain)
  return static_cast<double> (obs[lower_i]) +
         (frac * (static_cast<double> (obs[upper_i]) -
                  static_cast<double> (obs[lower_i])));
}

}  // namespace summary
