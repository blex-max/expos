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
constexpr inline std::optional<double> tail_jump (
    std::vector<T> obs     // copy for sort
) {
    std::sort (begin (obs), end (obs));
    if (obs.size() < 2 || obs.back() == 0)
        return std::nullopt;

    // can overflow
    // +1 avoids confusing values with obs of 0
    return static_cast<double> (obs.back() + 1)
           / static_cast<double> (obs[obs.size() - 2] + 1);
}

// Notes:
// - unsigned only because massive opposite signed integral elements will cause overflow
// TODO return boundaries with result
template <typename T>
    requires std::unsigned_integral<T> || std::floating_point<T>
constexpr inline std::optional<T> min_span_containing (
    std::vector<T> obs,     // copy for sort
    double         pt
) {
    // TODO consider whether these should be asserts
    // or runtime checks
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
struct endpoints1D {
    T min;
    T max;

    endpoints1D() = delete;
    endpoints1D (
        T a,
        T b
    ) {
        min = a > b ? b : a;
        max = a > b ? a : b;
    }

    T diff () const { return max - min; }
};

// chebyshev distance
template <typename T>
    requires std::unsigned_integral<T>
constexpr inline T ucheb (
    const endpoints1D<T> &a,
    const endpoints1D<T> &b
) {
    endpoints1D upper_pair{a.max, b.max};
    endpoints1D lower_pair{a.min, b.min};
    return std::max<T> (upper_pair.diff(), lower_pair.diff());
}

// for endpoints on a number line
// get the chebyshev distances between connected
// vertices on the mst, discarding the graph itself
// note:
// - implementation suboptimal, but fine for expected n
template <typename T>
    requires std::unsigned_integral<T>
inline std::vector<T> constexpr mst_cheb_dists (
    std::vector<endpoints1D<T>> &obs
) {
    using pwdist = struct {
        endpoints1D<T> p;
        uint64_t       dist;
    };
    using neighbours = std::vector<pwdist>;
    neighbours outside;
    std::transform (
        begin (obs),
        end (obs),
        std::back_inserter (outside),
        [] (const auto &a) {
            return pwdist{a, std::numeric_limits<T>::max()};
        }
    );
    std::vector<endpoints1D<T>> in_tree;
    std::vector<T>              dists;     // n-1 gaps
    while (!outside.empty()) {
        const auto &nn = outside.back();
        if (!in_tree.empty())
            dists.emplace_back (nn.dist);
        in_tree.emplace_back (nn.p);
        const auto &u = in_tree.back();
        outside.pop_back();     // avoid self-self

        // no-ops once empty
        std::transform (
            begin (outside),
            end (outside),
            begin (outside),
            [&u] (const auto &a) {
                const auto d = ucheb (a.p, u);
                return pwdist{
                    a.p,
                    std::min (a.dist, d)
                };     // get chebyshev distance
            }
        );
        sort (
            begin (outside),
            end (outside),
            [] (const auto &a, const auto &b) {
                return a.dist > b.dist;     // smallest to back
            }
        );
    }
    return dists;
}


// smallest chebyshev radius
// such that there exists some centre radius
// contains strictly more than fraction pt of all points
// TODO templatise
// TODO test
inline std::optional<uint64_t> min_cheb_radius_containing (
    const std::vector<endpoints1D<uint64_t>> &obs,
    double                                    pt
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

    uint64_t mind = std::numeric_limits<uint64_t>::max();
    for (size_t i = 0; i < obs.size(); ++i) {
        std::vector<uint64_t> idists;
        const auto           &q = obs[i];
        for (size_t j = 0; j < obs.size(); ++j) {
            if (i == j)
                continue;
            idists.emplace_back (ucheb (obs[j], q));
        }
        std::sort (begin (idists), end (idists));
        if (idists[nel - 2]
            < mind) {     // -2 as we include the query point
            mind = idists[nel - 2];
        }
    }
    return mind;
}


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
template <
    typename StatT,
    typename ObsT>
inline stat_eval_s sim_to_bg (
    StatT             ev_stat,
    size_t            n_ev_obs,
    std::vector<ObsT> total_obs,     // intentional copy
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

    // TODO also compare:
    // i) event obs to uniform model
    // ii) background/total obs to uniform model
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
