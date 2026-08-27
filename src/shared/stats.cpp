#include "stats.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string_view>

double entropy_lz76 (std::string_view s)
{
  const auto nChar = s.size();
  if (nChar == 0) {
    return 0;
  }
  if (nChar == 1) {
    return 1;
  }

  uint64_t i = 0;           // index into prefix substring
  uint64_t lzCmplx = 1;     // phrase count
  uint64_t prefixLen =
      1;  // length of substring which has been assessed and is now memory
  uint64_t compLen =
      1;  // length of the candidate substring component
  uint64_t cycleMax = compLen;  // longest match for a cycle
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
