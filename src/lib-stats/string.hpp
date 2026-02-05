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


/*
periodic_rle: greedy segmentation of a string into maximal exact tandem-repeat
segments with bounded period.

While it's hard to coerce this into a direct metric on which we can flag variants
high microhomology content as assessed by this approach is very strongly correlated
with LZ plummeting, so it's a nice verification.

Intention:
- Partition the string left-to-right into contiguous segments.
- At each segment start i, consider exact tandem repeats whose period k is in
  [1, max_k].
- Choose the period that maximizes contiguous repeat length in bases
  (run_period * k), and emit that length as the segment size.
- Length-1 segments are intentional and represent explicit breaks in periodic
  structure at the tested scale.

Assumptions / scope:
- Repeats are exact and phase-anchored at the segment start i.
- Only low-period structure (k ≤ max_k) is considered; longer-period duplication
  is ignored by design.
- Greedy tiling is used; overlapping or alternative coverings are not explored.

Invariants:
- Segment lengths are positive and sum to ≤ string length.
- Candidate comparisons are bounds-safe: full k-length units are always tested.
- “Best” is defined in base pairs, not repeat units, avoiding bias toward small k.
*/
inline std::vector<size_t> periodic_rle(const std::string& s, size_t max_k)
{
    std::vector<size_t> seg_lens;
    const size_t n = static_cast<size_t>(s.size());
    if (n == 0) return seg_lens;

    size_t i = 0;
    while (i < n) {
        size_t best_run_bp = 1;                 // default: 1 bp break
        const size_t k_limit = std::min(max_k, n - i);

        for (size_t k = 1; k <= k_limit; ++k) {
            const std::string unit = s.substr(i, k);

            size_t run_period = 1;
            while (i + (run_period + 1) * k <= n) {
                const std::string cand = s.substr(i + run_period * k, k);
                if (cand == unit) ++run_period;
                else break;
            }

            if (run_period >= 2) {
                const auto run_bp = run_period * k;
                if (run_bp > best_run_bp) best_run_bp = run_bp;
            }
        }

        seg_lens.push_back(best_run_bp);
        i += best_run_bp;
    }

    return seg_lens;
}





} // end namespace

