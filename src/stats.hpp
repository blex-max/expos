#include <algorithm>
#include <cassert>
#include <concepts>
#include <math.h>
#include <ranges>
#include <stdexcept>
#include <vector>

template <std::integral Int>
inline std::vector<double> as_double (const std::vector<Int> &v) {
    if (v.empty())
        return {};

    constexpr auto dlim = (1ULL << 53);

    auto oor = [] (Int a) {
        if constexpr (std::is_signed_v<Int>) {
            return (a > static_cast<Int> (dlim))
                   || (a < -static_cast<Int> (dlim));
        } else {
            return a > dlim;
        }
    };

    if (std::ranges::any_of (v, oor)) {
        throw std::runtime_error ("cannot safely convert to double");
    }

    auto vv = v | std::views::transform ([] (auto a) {
                  return static_cast<double> (a);
              });
    return {begin (vv), end (vv)};
}

// get deviations for each element
// from a value
inline std::vector<double> devs (
    const std::vector<double> &v,
    double                     from
) {
    auto vv = v | std::views::transform ([&from] (auto a) {
                  return fabs (a - from);
              });
    return {begin (vv), end (vv)};
}

// preconditions:
// vector is sorted
// vector is not empty
inline double median_from_sorted (const std::vector<double> &v) {
    assert (!v.empty());
    auto mel = v.size() / 2;
    if (v.size() % 2 != 0) {     // is odd
        return v[mel];
    }
    return (v[mel - 1] + v[mel]) / 2;     // is even
}

// calculate median absolute deviation
template <std::integral Int>
inline std::optional<double> mad (const std::vector<Int> &v) {
    if (v.empty())
        return std::nullopt;
    auto dv = as_double (v);     // should handle any integral type
    std::ranges::sort (dv);      // ascending sort
    const auto med  = median_from_sorted (dv);
    auto       devv = devs (dv, med);
    std::ranges::sort (devv);
    return median_from_sorted (devv);
}
