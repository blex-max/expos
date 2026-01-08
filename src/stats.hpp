#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <math.h>
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


inline double medianNN (const PairMatrix &pwd) {
    assert (pwd.dim() > 1);
    std::vector<uint64_t> nndv;
    for (size_t row = 0; row < pwd.dim(); ++row) {
        auto min_nnd = std::numeric_limits<uint64_t>::max();
        for (size_t col = 0; col < pwd.dim(); ++col) {
            if (row == col)
                continue;     // skip self-self
            const auto nnd = pwd.get (row, col);
            if (nnd < min_nnd) {
                min_nnd = nnd;
            }
        }
        nndv.push_back (min_nnd);
    }
    return *percentile (nndv, 0.5);     // median
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
// clang-format off
struct sim_to_bgConfig {
const size_t       nsim              = 2500;
const size_t       event_obs_ext_min = 0;
std::mt19937      rng               = {};
};
template <
    typename ObsT,
    typename StatFn,
    typename EffFn,
    typename CmpFn,
    typename StatT = std::invoke_result_t<StatFn&, const std::vector<ObsT>>
>
requires
    std::invocable<StatFn&, const std::vector<ObsT>> &&
    std::invocable<CmpFn&, StatT, StatT> &&
    std::invocable<EffFn&, StatT, const std::vector<StatT>>
// clang-format on
// TODO SEED
inline stat_eval_s sim_to_bg (
    StatT             ev_stat,
    size_t            n_ev_obs,
    std::vector<ObsT> total_obs,     // intentional copy (or move-in)
    StatFn          &&statfn,
    CmpFn           &&statcmp,
    EffFn           &&efffn,
    sim_to_bgConfig  &conf
) {
    stat_eval_s res;
    assert (nsim > 0);
    if (n_ev_obs < 2 || n_ev_obs < conf.event_obs_ext_min) {
        res.err = "INSUFF_OBS";
        return res;
    }
    if (
        total_obs.size()
        < n_ev_obs
              * 2     // at a bare minimum, we want 2x more total samples than bg
    ) {
        res.err = "INSUFF_BG";
        return res;
    }

    std::vector<StatT> sim_vals;
    size_t             sim_count = 0;
    for (size_t i = 0; i < conf.nsim; ++i) {
        std::shuffle (begin (total_obs), end (total_obs), conf.rng);
        const auto sv = statfn (
            std::vector (
                total_obs.begin(),
                total_obs.begin() + n_ev_obs
            )
        );
        if (statcmp (ev_stat, sv)) {
            ++sim_count;
        }
        sim_vals.push_back (sv);
    }

    // report effect size
    // if eff_sz is large then we can
    // get away with a low number of samples
    // if not it's just noise
    res.eff_sz = efffn (ev_stat, sim_vals);

    // TODO "power analysis"

    res.pval = static_cast<double> (sim_count + 1)
               / static_cast<double> (conf.nsim + 1);
    return res;
}


// FOR REF COMPLEXITY
// length normalised kolgomorov complexity via lempel-ziv 76
inline double nk_lz76 (const std::string &s) {
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
