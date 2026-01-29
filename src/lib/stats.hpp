#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <math.h>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>

template <std::signed_integral T>
constexpr inline uint64_t as_uint (const T &a) {
    if (a < 0) {
        throw std::runtime_error (
            "cannot convert negative value to uint"
        );
    }
    return static_cast<uint64_t> (a);
}


constexpr inline std::optional<double>
mean (const std::vector<double> &v) {
    if (v.empty())
        return std::nullopt;

    long double sum = 0.0L;
    for (const auto &x : v) {
        sum += static_cast<long double> (x);
    }

    return static_cast<double> (sum) / static_cast<double> (v.size());
}


template <typename T>
requires std::unsigned_integral<T> || std::floating_point<T>
constexpr inline std::optional<double> percentile (
    std::vector<T> obs,
    double         pt
) {
    assert (pt > 0 && pt < 1);

    if (obs.empty())
        return std::nullopt;

    if (obs.size() == 1)
        return obs[0];

    std::sort (begin (obs), end (obs));

    double pi = static_cast<double> (obs.size() - 1)
                * pt;     // 0 indexed rank
    auto lower   = floor (pi);
    auto frac    = pi - lower;
    auto upper_i = static_cast<size_t> (ceil (pi));
    auto lower_i = static_cast<size_t> (lower);

    if (lower_i == upper_i) {
        return obs[lower_i];
    }

    // linear interpolation
    // TODO those casts could lose precision
    // need safe double conversion
    return static_cast<double> (obs[lower_i])
           + (frac
              * (static_cast<double> (obs[upper_i])
                 - static_cast<double> (obs[lower_i])));
}


// object described by two
// coordinates on the same axis
struct line_seg {
    uint64_t lmost;
    uint64_t rmost;

    line_seg() = delete;
    line_seg (
        uint64_t a,
        uint64_t b
    ) {
        lmost = a > b ? b : a;
        rmost = a > b ? a : b;
    }

    uint64_t diff () const { return rmost - lmost; }
};


// 2D symmetric square matrix via vector
// rows are contiguous in vector
class PairMatrix {
  private:
    std::vector<uint64_t> mat;
    const size_t          dim_;

    PairMatrix (
        std::vector<uint64_t> v,
        size_t                dim
    )
        : mat (v),
          dim_ (dim) {
        assert (dim > 1);
        assert (v.size() == (dim * dim));
    }

  public:
    PairMatrix() = delete;

    auto dim () const noexcept { return dim_; }
    auto size () const noexcept { return mat.size(); }

    auto get (
        size_t i,
        size_t j
    ) const {
        assert (i < dim_);
        assert (j < dim_);
        assert (!mat.empty());
        assert (mat.size() == (dim() * dim()));
        return mat[(i * dim_) + j];
    }

    const auto &get1D () const noexcept { return mat; }
    auto        copy1D () const noexcept { return mat; }

    // assumes symmetric distance
    // clang-format off
    template <typename U, typename F>
        requires std::invocable<F &, const U &, const U &>
                 && std::same_as<std::invoke_result_t<F &,const U &,const U &>, uint64_t>
    // clang-format on
    static std::optional<PairMatrix> from_sample (
        const std::vector<U> &obs,
        F                   &&sym_pairfn
    ) {
        assert (!obs.empty());
        const auto dim = obs.size();
        if (dim < 2)
            return std::nullopt;
        const auto            nel = dim * dim;
        std::vector<uint64_t> in (nel);     // nel-long vector
        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < (i + 1); ++j) {
                const auto val    = sym_pairfn (obs[i], obs[j]);
                in[(i * dim) + j] = val;
                if (j != i)
                    in[(j * dim) + i] = val;
            }
        }
        return PairMatrix{in, dim};
    }
};


// Two-stage trimmed mean to assess local spatial structure.
// What is the average local density among the densest <second_cutoff_pt> of points?
// Where local density is defined as the mean of distances to the 0.25*n closest
// neigbhours of a given point (excluding self-self).
// The cutoff parameter allows tuning the weighting of the least crowded observations
inline double trimmed_mean_of_lower_tails (const PairMatrix &pwd, double second_cutoff_pt) {
    const auto n = pwd.dim();  // square matrix
    assert(n > 1);
    assert(cutoff_pt > 0.0 && cutoff_pt <= 1.0);

    const auto tail_k = (n + 2) / 4;  // 25% smallest observations of a vector of size n - 1

    std::vector<double> lower_tail_means;
    std::vector<uint64_t> row_dists;
    row_dists.reserve(n - 1);
    uint64_t dist_sum= 0;
    for (size_t row = 0; row < n; ++row) {
        row_dists.clear();
        for (size_t col = 0; col < n; ++col) {
            if (row == col) {
                continue;     // skip self-self
            }
            row_dists.push_back(pwd.get(row, col));
        }
        std::nth_element(
            begin(row_dists),
            begin(row_dists) + (tail_k - 1),
            end(row_dists)
        );

        dist_sum = 0;
        for (size_t i = 0; i < tail_k; ++i) {
            dist_sum += row_dists[i];
        }
        lower_tail_means.push_back(
            static_cast<double>(dist_sum)
            / static_cast<double>(tail_k)
        );
    }

    {
        const auto cutoff_k =
            static_cast<size_t>(
                ceil(
                    static_cast<double>(n * second_cutoff_pt)
                ));
        std::nth_element(
            begin(lower_tail_means),
            begin(lower_tail_means) + (cutoff_k - 1),
            end(lower_tail_means)
        );

        double sum = 0;
        for (size_t i = 0; i < cutoff_k; ++i) {
            sum += lower_tail_means[i];
        }
        return sum / static_cast<double>(cutoff_k);
    }
}


struct stat_eval_s {
    std::optional<double> eff_sz = std::nullopt;
    std::optional<double> pval   = std::nullopt;
    std::string           err    = "";
};
// are event supporting observations
// meaningfully different than total observations
// for a given statistic?
// compare via simulation: if we take random samples
// of the total observations of the same size
// as the true event supporting sample,
// how often do we get a lt/gt value for the
// statistic in question (pvalue), and
// how large is the effect size.
template <
    typename StatT,
    typename StatFn,
    typename DrawFn,
    typename EffFn,
    typename ObsT = typename std::invoke_result_t<DrawFn&>::value_type
>
requires
    std::same_as<
        std::invoke_result_t<DrawFn&>,
        std::vector<ObsT>
    > &&
    std::same_as<
        std::invoke_result_t<StatFn&, const std::vector<ObsT> &>,
        StatT
    > &&
    std::same_as<
        std::invoke_result_t<EffFn&, StatT, const std::vector<StatT> &>,
        double
    >
inline stat_eval_s sim_to_bg (
    StatT    ev_stat,
    DrawFn &&drawfn,
    StatFn &&statfn,
    EffFn  &&efffn,
    size_t   nsim=2500
) {
    stat_eval_s res;
    assert (nsim > 0);
    std::vector<ObsT> sim_obs;

    std::vector<StatT> draw_results;
    draw_results.reserve(nsim);

    size_t sim_count_le = 0;  // count of simulated stats that are <= observed stat
    for (size_t i = 0; i < nsim; ++i) {
        sim_obs.clear();
        sim_obs = drawfn();

        const auto draw_stat = statfn (sim_obs);
        if (draw_stat <= ev_stat) {
            ++sim_count_le;
        }
        draw_results.push_back (draw_stat);
    }

    // report effect size
    // if eff_sz is large then we can
    // get away with a low number of samples
    // if not it's just noise
    res.eff_sz = efffn (ev_stat, draw_results);

    // TODO "power analysis"

    // two sided p val
    res.pval = 2
               * static_cast<double> (std::min(sim_count_le + 1, nsim - sim_count_le + 1))
               / static_cast<double> (nsim + 1);

    // one sided
    // res.pval = static_cast<double> (sim_count_le + 1)
    //            / static_cast<double> (nsim + 1);
    return res;
}


// +1 removes confusing values when 0,
// log makes effect size symmetric around 0
// log2 means -1 = half the size of background
// +1 = double the size of background
inline auto log2_effsz (const double &ev, const std::vector<double> &simv) {
    return log2 ((ev + 1) / (*mean (simv) + 1));
}

// get a random sample without replacement of
// size n from input obs vector.
template <typename T>
inline std::vector<T> subsample_wo_replace (
    const std::vector<T>& obs,
    size_t n,
    std::mt19937 &rng
) {
    size_t nobs = obs.size();
    assert (n < nobs);

    // partial Fisher–Yates shuffle an index vector so that
    // idx[0..n] is a uniform sample without replacement.
    std::vector<size_t> all_idx(nobs);
    std::iota(begin(all_idx), end(all_idx), 0);

    // shuffle
    for (size_t i = 0; i < n; ++i) {
        std::uniform_int_distribution<size_t> dist(i, nobs - 1);
        size_t j = dist(rng);
        std::swap(all_idx[i], all_idx[j]);
    }

    // get obs
    std::vector<T> out;
    for (size_t k = 0; k < n; ++k) {
        out.push_back(obs[all_idx[k]]);
    }
    return out;
}

// FOR REF COMPLEXITY
// lempel-ziv 76 entropy rate (bits per base)
inline double entropy_lz76 (const std::string &s) {
    assert (!s.empty());
    size_t slen = s.size();
    size_t n_phrase = 1;     // number of phrases (complexity), starts at 1 (first char is a phrase)
    size_t frontier_i = 1;     // start index of current phrase in s
    size_t memory_search_i = 0;     // index over already processed length of s
    size_t match_len = 1;
    size_t max_match_len = 1;     // maximum match length for phrase so far
    bool stop = false;

    if (slen < 2)
        stop = true;

    while (!stop) {
        // compare chars
        if (s[frontier_i + match_len - 1]
            != s[memory_search_i + match_len - 1]) {     // mismatch

            if (match_len > max_match_len) {
                max_match_len = match_len;
            }

            ++memory_search_i;     // next memory char

            if (memory_search_i == frontier_i) {     // search complete

                ++n_phrase;
                frontier_i += max_match_len;     // jump to start of next phrase

                if (frontier_i + 1 > slen) {
                    stop = true;
                } else {
                    // reset
                    memory_search_i = 0;
                    match_len       = 1;
                    max_match_len   = 1;
                }
            } else {
                // restart search for matches from next memory position (++memory_search_i)
                match_len = 1;
            }
        } else {             // match
            ++match_len;     // extend current match

            if (frontier_i + match_len > slen) {     // reached end
                ++n_phrase;
                stop = true;
            }
        }
    }

    // *log2(n) == kolmogorov complexity
    // / length normalise
    // result is bits (entropy) per character
    return (static_cast<double> (n_phrase) * log2 (slen))
           / static_cast<double> (slen);
}
