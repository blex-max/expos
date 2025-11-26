#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <math.h>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <vector>

template <std::signed_integral T>
constexpr inline uint64_t as_uint (const T &a) {
    if (a < 0) {
        throw std::runtime_error ("cannot convert negative value to uint");
    }
    return static_cast<uint64_t> (a);
}

template <std::signed_integral T>
constexpr inline std::vector<uint64_t>
as_uintv (const std::vector<T> &v) {
    if (v.empty())
        return {};

    return v | std::views::transform ([] (const auto a) {
               if (a < 0) {
                   throw std::runtime_error (
                       "cannot convert negative value to uint"
                   );
               }
               return static_cast<uint64_t> (a);
           })
           | std::ranges::to<std::vector<uint64_t>>();
}

template <std::unsigned_integral T>
constexpr inline std::vector<double>
as_doublev (const std::vector<T> &v) {
    if (v.empty())
        return {};

    constexpr auto dlim = (1ULL << 53);

    return v | std::views::transform ([] (auto a) {
               if (a > dlim) {
                   throw std::runtime_error (
                       "cannot safely convert to double"
                   );
               }
               return static_cast<double> (a);
           })
           | std::ranges::to<std::vector<double>>();
}

// get deviations for each element
// from a value
constexpr inline std::vector<double> devs (
    const std::vector<double> &obs,
    double                     from
) {
    auto vv = obs | std::views::transform ([&from] (auto a) {
                  return fabs (a - from);
              });
    return {begin (vv), end (vv)};
}

// preconditions:
// vector is sorted
// Notes:
// - unsigned only because massive opposite signed integral elements will cause overflow
template <typename T>
    requires std::unsigned_integral<T> || std::floating_point<T>
constexpr inline std::optional<double> percentile_from_sorted (
    const std::vector<T> &obs,
    double                pt
) {
    assert (pt > 0 && pt < 1);
    assert (std::is_sorted (begin (obs), end (obs)));

    if (obs.empty())
        return std::nullopt;

    if (obs.size() == 1)
        return obs[0];

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

// calculate median absolute deviation
// precondition:
// - vector is sorted
template <std::unsigned_integral T>
constexpr inline std::optional<double> mad (const std::vector<T> &obs) {
    assert (std::is_sorted (begin (obs), end (obs)));
    if (obs.empty())
        return std::nullopt;
    auto dv = as_doublev (obs);     // should handle any integral type
    // std::ranges::sort (dv);        // ascending sort
    const auto med  = percentile_from_sorted (dv, 0.5);
    auto       devv = devs (dv, med.value());
    std::ranges::sort (devv);
    return percentile_from_sorted (devv, 0.5);
}


//--- GAP BASED STATS ---//
// sparse integer data where the question of interest is:
// is there tight clustering compared to what you would expect
// from a random distribution?
// questions to ask of the data
// what is the largest number of observations to occur in a small (e.g 6bp) window of the range?
// what is the shortest span of the range to contain 50%, 90% of the observations
// - the above two are really easy to calculate with a sorted vector of observations
// IQR of gap size
// - all very relevant and robust to a low number of obs
// compare all to simulation
// also similar evaluation of qpos as well as template start/end
// - (not size...? because size is not a number line in the same way)


// calculate gaps from ascending sorted vector
// Notes:
// - unsigned only because massive opposite signed integral elements will cause overflow
template <typename T>
    requires std::unsigned_integral<T> || std::floating_point<T>
constexpr inline std::vector<T> gaps (const std::vector<T> &v) {
    assert (std::is_sorted (begin (v), end (v)));
    if (v.empty())
        return {};
    std::vector<T> out;
    auto           prev = v[0];
    for (const auto &e : std::ranges::subrange (begin (v) + 1, end (v))) {
        assert (e >= prev);     // sorted precondition
        // BUG (and elsewhere) - massive opposite signed elements will cause overflow
        out.emplace_back (e - prev);
        prev = e;
    }
    return out;
}

// largest fraction of total observations
// that fall within a particular window width
// preconditions:
// - vector is ascending sorted
// Notes:
// - unsigned only because massive opposite signed integral elements will cause overflow
// TODO return boundaries with result
template <typename T>
    requires std::unsigned_integral<T> || std::floating_point<T>
constexpr inline std::optional<double> max_obs_within (
    const std::vector<T> &obs,
    T                     span
) {
    assert (std::is_sorted (begin (obs), end (obs)));
    if (obs.empty()) {
        return std::nullopt;
    }
    if (obs.size() == 1) {
        return 1;
    }

    size_t max_within  = 1;
    size_t within_this = 1;
    auto   prev        = obs[0];
    auto   spanl       = prev;

    for (const auto &e :
         std::ranges::subrange (begin (obs) + 1, end (obs))) {
        assert (e >= prev);     // sorted precondition
        auto diff = e - spanl;
        if (diff < span) {
            ++within_this;
        } else {
            if (within_this > max_within) {
                max_within = within_this;
            }
            within_this = 1;     // reset;
            spanl       = e;
        }
    }

    return static_cast<double> (max_within)
           / static_cast<double> (obs.size());
}

// preconditions:
// - vector is ascending sorted
// Notes:
// - unsigned only because massive opposite signed integral elements will cause overflow
// TODO return boundaries with result
template <typename T>
    requires std::unsigned_integral<T> || std::floating_point<T>
constexpr inline std::optional<T> min_span_containing (
    const std::vector<T> &obs,
    double                pt
) {
    // TODO consider whether these should be asserts
    // or runtime checks
    assert (pt > 0 && pt <= 1);
    assert (std::is_sorted (begin (obs), end (obs)));

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


// summary stats
// from which we may calculate IQR
// and other measures
template <typename T>
    requires std::integral<T> || std::floating_point<T>
struct summary5_s {
    T      min;
    double q25;
    double q50;
    double q75;
    T      max;
};
template <typename T>
    requires std::unsigned_integral<T>     // uint for percentile func
             || std::floating_point<T>
constexpr inline summary5_s<T> summary_stats (const std::vector<T> &v) {
    assert (std::is_sorted (begin (v), end (v)));
    assert (!v.empty());
    const auto [pvmin, pvmax] = std::minmax_element (begin (v), end (v));
    return {
        *pvmin,
        percentile_from_sorted (v, 0.25).value(),
        percentile_from_sorted (v, 0.5).value(),
        percentile_from_sorted (v, 0.75).value(),
        *pvmax
    };
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
                return pwdist{a.p, std::min (a.dist, d)};     // get chebyshev distance
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

// precondition:
// - vector is ascending sorted
template <typename T>
    requires std::unsigned_integral<T> || std::floating_point<T>
constexpr inline std::optional<double>
tail_jump (const std::vector<T> &obs) {
    assert (std::is_sorted (begin (obs), end (obs)));
    if (obs.size() < 2 || obs.back() == 0)
        return std::nullopt;

    // BUG will fail if divisor >>>
    return static_cast<double> (obs.back())
           / static_cast<double> (obs[obs.size() - 2]);
}

// not useful since won't detect a single cluster (the data would be regular)
// template <std::integral Int>
// inline std::optional<double> cv_gaps (const std::vector<Int> &v) {
//     if (v.empty())
//         return std::nullopt;
//     auto dv = as_double (v);     // should handle any integral type
//     std::ranges::sort (dv);      // ascending sort
//     auto gv    = gaps (dv);
//     auto gmean = mean (gv);
//     auto gsd   = pop_stddev (gv, gmean);
//     return coeff_var (gsd, gmean);
// }
