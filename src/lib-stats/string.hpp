#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string_view>
#include <vector>

namespace string_stats {
// lempel-ziv 76 entropy rate (bits per char)
// used for calculating reference complexity
inline double entropy_lz76 (std::string_view s) {
    const auto nchar = s.size();
    if (nchar == 0) return 0;
    if (nchar == 1) return 1;

    size_t i = 0;            // index into prefix substring
    size_t lz_cmplx = 1;     // phrase count
    size_t prefix_len = 1;   // length of substring which has been assessed and is now memory
    size_t comp_len = 1;     // length of the candidate substring component
    size_t cycle_max = comp_len;  // longest match for a cycle
    while (prefix_len + comp_len <= nchar) {
        if (s[i + comp_len - 1] == s[prefix_len + comp_len - 1])
        {  // char at prefix == char at frontier
            ++comp_len;
        }
        else  // mismatch
        {
            cycle_max = std::max(comp_len, cycle_max);
            ++i;
            if (i == prefix_len)
            {  // memory completely searched
                ++lz_cmplx;
                prefix_len = prefix_len + cycle_max;
                comp_len = 1;
                i = 0;
                cycle_max = comp_len;
            }
            else
            {
                comp_len = 1;
            }
        }
    }

    // final phrase on loop exit
    if (comp_len != 1) {
        ++lz_cmplx;
    }

    // length normalise
    // result is bits (entropy) per character
    return (static_cast<double> (lz_cmplx) * log2 (nchar))
           / static_cast<double> (nchar);
}


inline std::vector<size_t> rle (std::string_view s) {
  std::vector<size_t> out;

  size_t nchar = s.size();
  size_t i = 0;
  size_t run_len = 1;
  for (;(i + run_len) < nchar;) {
    if (s[i + run_len] == s[i]) {
      ++run_len;
    } else {
      out.push_back (run_len);
      i += run_len; // next char
      run_len = 1;
    }
  }

  // handle final
  out.push_back (run_len);

  return out;
}

} // end namespace

