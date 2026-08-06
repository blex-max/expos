#include "stats.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
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
