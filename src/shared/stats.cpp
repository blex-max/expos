#include "stats.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>
#include <vector>

double entropy_lz76 (std::string_view s)
{
  const auto nChar = s.size();
  if (nChar == 0) {
    return 0;
  }
  if (nChar == 1) {
    return 1;
  }

  size_t i = 0;           // index into prefix substring
  size_t lzCmplx = 1;     // phrase count
  size_t prefixLen =
      1;  // length of substring which has been assessed and is now memory
  size_t compLen =
      1;  // length of the candidate substring component
  size_t cycleMax = compLen;  // longest match for a cycle
  while (prefixLen + compLen <= nChar) {
    if (s[i + compLen - 1] ==
        s[prefixLen + compLen -
          1]) {  // char at prefix == char at frontier
      ++compLen;
    }
    else  // mismatch
    {
      cycleMax = std::max (compLen, cycleMax);
      ++i;
      if (i == prefixLen) {  // memory completely searched
        ++lzCmplx;
        prefixLen = prefixLen + cycleMax;
        compLen = 1;
        i = 0;
        cycleMax = compLen;
      }
      else {
        compLen = 1;
      }
    }
  }

  // final phrase on loop exit
  if (compLen != 1) {
    ++lzCmplx;
  }

  // length normalise
  // result is bits (entropy) per character
  return (static_cast<double> (lzCmplx) * log2 (nChar)) /
         static_cast<double> (nChar);
}

std::optional<double> mean (const std::vector<double>& v)
{
  if (v.empty()) {
    return std::nullopt;
  }
  long double sum = 0.0L;
  for (const double x : v) {
    sum += static_cast<long double> (x);
  }
  return static_cast<double> (
      sum / static_cast<long double> (v.size())
  );
}

std::size_t count_pairs_within_1d (
    std::vector<int32_t> pts, uint64_t radius
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

double ripley_k (
    std::size_t pairsWithin, std::size_t n, double intensity
)
{
  // each unordered pair is two ordered neighbour relations; mean per
  // point is /n; intensity normalisation is 1/intensity.
  constexpr double k_orderedPerUnordered = 2.0;
  return (1.0 / intensity) * (k_orderedPerUnordered *
                              static_cast<double> (pairsWithin) /
                              static_cast<double> (n));
}
