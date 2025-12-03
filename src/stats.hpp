#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <math.h>
#include <optional>
#include <random>
#include <ranges>
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


//--- GAP BASED STATS ---//
// sparse integer data where the question of interest is:
// is there tight clustering compared to what you would expect
// from a random distribution?
// questions to ask of the data
// what is the shortest span of the range to contain 50%, 90% of the observations
// - easy to calculate and robust to a low number of obs
// compare to simulation
// - (size...? size is a number line in the same way)


// calculate gaps from ascending sorted vector
// Notes:
// - unsigned only because massive opposite signed integral elements will cause overflow
template <typename T>
    requires std::unsigned_integral<T> || std::floating_point<T>
constexpr inline std::vector<T> gaps (std::vector<T> v) {
    std::sort (begin (v), end (v));
    if (v.empty())
        return {};
    std::vector<T> out;
    auto           prev = v[0];
    for (const auto &e :
         std::ranges::subrange (begin (v) + 1, end (v))) {
        assert (e >= prev);     // sorted precondition
        // BUG (and elsewhere) - massive opposite signed elements will cause overflow
        out.emplace_back (e - prev);
        prev = e;
    }
    return out;
}

template <typename T>
    requires std::unsigned_integral<T> || std::floating_point<T>
constexpr inline std::optional<double>
tail_jump (const std::vector<T> &obs) {
    if (obs.size() < 2)
        return std::nullopt;

    T max1 = std::numeric_limits<T>::lowest();
    T max2 = std::numeric_limits<T>::lowest();

    for (const auto &x : obs) {
        if (x > max1) {
            max2 = max1;
            max1 = x;
        } else if (x > max2) {
            max2 = x;
        }
    }

    // could overflow
    // +1 avoids confusing values with obs of 0
    return static_cast<double> (max1 + 1)
           / static_cast<double> (max2 + 1);
}

// Notes:
// - unsigned only because massive opposite signed integral elements will cause overflow
// TODO return boundaries with result
template <typename T>
    requires std::unsigned_integral<T>
constexpr inline std::optional<T> min_span_containing (
    std::vector<T> obs,     // copy for sort, TODO remove
    double         pt
) {
    assert (pt > 0 && pt <= 1);

    if (obs.size() < 2) {
        return std::nullopt;
    }

    // what number of elements constitutes
    // at least pt% of the obs
    auto nel = static_cast<size_t> (
        ceil (static_cast<double> (obs.size()) * pt)
    );
    if (nel <= 1) {
        return std::nullopt;     // single element is not a statistically valid span
    }

    std::sort (begin (obs), end (obs));

    T min_span = std::numeric_limits<T>::max();
    for (size_t i = 0; (i + nel - 1) < obs.size(); ++i) {
        const auto &e0 = obs[i];
        const auto &en = obs[i + nel - 1];
        assert (en >= e0);
        const auto this_span = en - e0;
        if (this_span < min_span)
            min_span = this_span;
    }

    return min_span;
}

// object described by two
// coordinates on the same axis
template <typename T>
    requires std::unsigned_integral<T>
struct line_seg {
    T lmost;
    T rmost;

    line_seg() = delete;
    line_seg (
        T a,
        T b
    ) {
        lmost = a > b ? b : a;
        rmost = a > b ? a : b;
    }

    T diff () const { return rmost - lmost; }
};

// chebyshev distance
template <typename T>
    requires std::unsigned_integral<T>
constexpr inline T ucheb (
    const line_seg<T> &a,
    const line_seg<T> &b
) {
    line_seg upper_pair{a.rmost, b.rmost};
    line_seg lower_pair{a.lmost, b.lmost};
    return std::max<T> (upper_pair.diff(), lower_pair.diff());
}

// 2D symmetric square matrix via vector
// rows are contiguous in vector
template <typename T>
    requires std::integral<T> || std::floating_point<T>
class PairMatrix {
  private:
    std::vector<T> mat;
    const size_t   dim_;

    PairMatrix (
        std::vector<T> v,
        size_t         dim
    )
        : mat (v),
          dim_ (dim) {
        assert (dim > 1);
        assert ((v.size() % dim) == 0);
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
        return mat[(i * dim_) + j];
    }

    const auto &get1D () const noexcept { return mat; }
    auto        copy1D () const noexcept { return mat; }

    // assumes symmetric distance
    // clang-format off
    template <typename U, typename F>
        requires std::invocable<F &, const U &, const U &>
                 && std::same_as<std::invoke_result_t<F &,const U &,const U &>, T>
    // clang-format on
    static PairMatrix<T> from_sample (
        const std::vector<U> &obs,
        F                   &&sym_pairfn
    ) {
        const auto     dim = obs.size();
        const auto     nel = dim * dim;
        std::vector<T> in (nel);     // nel-long vector
        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < (i + 1); ++j) {
                const auto val    = sym_pairfn (obs[i], obs[j]);
                in[(i * dim) + j] = val;
                if (j != i)
                    in[(j * dim) + i] = val;
            }
        }
        return {in, dim};
    }

    // span of smallest cluster of points that
    // contains strictly more than fraction pt of all points
    // assumes symmetric
    std::optional<T> min_span_containing (double pt) const {
        assert (pt > 0 && pt <= 1);
        const auto nobs = dim();

        // what number of elements constitutes
        // at least pt% of the obs
        auto ptel = static_cast<size_t> (
            ceil (static_cast<double> (nobs) * pt)
        );
        if (ptel < 2) {
            return std::nullopt;     // single element is not a statistically valid span
        }

        // find tightest cluster
        struct {
            size_t              query_i;
            std::vector<size_t> member_j;
            T                   rad = std::numeric_limits<T>::max();
        } min_clust;
        std::vector<size_t> rowj (nobs);
        std::iota (begin (rowj), end (rowj), 0);
        for (size_t row = 0; row < nobs; ++row) {
            std::nth_element (
                begin (rowj),
                begin (rowj) + (ptel - 1),     // sort up to i
                end (rowj),
                [&row, this] (const auto &a, const auto &b) {
                    return get (row, a)
                           < get (row, b);     // ascending sort
                }
            );
            const auto pt_dist = get (row, rowj[ptel - 1]);
            if (pt_dist < min_clust.rad) {
                min_clust.query_i  = row;
                min_clust.member_j = std::vector (
                    begin (rowj),
                    begin (rowj) + ptel
                );
                min_clust.rad = pt_dist;
            }
        }

        T span = min_clust.rad;
        // find largest p2p distance
        // within tightest cluster
        for (size_t i = 0; i < ptel; ++i) {
            for (size_t j = 0; j < ptel; ++j) {
                if (i == j)
                    continue;
                const auto ijd = get (
                    min_clust.member_j[i],
                    min_clust.member_j[j]
                );
                if (ijd > span) {
                    span = ijd;
                }
            }
        }
        return span;
    }
};


struct stat_eval_s {
    std::optional<double> eff_sz = std::nullopt;
    std::optional<double> pval   = std::nullopt;
    // int errcode;
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
// TODO templatise better
template <
    typename StatT,
    typename ObsT>
inline stat_eval_s sim_to_bg (
    StatT             ev_stat,
    size_t            n_ev_obs,
    std::vector<ObsT> total_obs,     // intentional copy, TODO remove
    std::function<StatT (const std::vector<ObsT> &)> statfn,
    std::function<bool (StatT)>                      statcmp,
    size_t                                           nsim = 1000,
    size_t event_obs_ext_min                              = 5
) {
    stat_eval_s res;
    if (
        n_ev_obs < 2 || n_ev_obs < event_obs_ext_min
        || total_obs.size()
               < n_ev_obs
                     * 2     // at a bare minimum, we want 2x more total samples than bg
    ) {
        return res;
    }

    // report effect size
    // if eff_sz is large then we can
    // get away with a low number of samples
    // if not it's just noise
    // +1 removes confusing values when 0,
    // log makes effect size symmetric around 0
    // log2 means -1 = half the size of background
    // +1 = doulbe the size of background
    res.eff_sz = log2 (
        static_cast<double> (ev_stat + 1)
        / static_cast<double> (statfn (total_obs) + 1)
    );

    // TODO "power analysis"

    std::random_device rd;
    std::mt19937       g (rd());

    size_t sim_count = 0;
    for (size_t i = 0; i < nsim; ++i) {
        std::shuffle (begin (total_obs), end (total_obs), g);
        const auto sv = statfn (
            std::vector (
                total_obs.begin(),
                total_obs.begin() + n_ev_obs - 1
            )
        );
        if (statcmp (sv)) {
            ++sim_count;
        }
    }

    res.pval = static_cast<double> (sim_count + 1)
               / static_cast<double> (nsim + 1);
    return res;
}


// FOR REF COMPLEXITY
// length normalised kolgomorov complexity via lempel-ziv 76
// TODO test
inline double nk_lz76 (const std::string &s) {
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
