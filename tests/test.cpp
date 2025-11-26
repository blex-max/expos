#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <set>

#include "stats.hpp"

TEST_CASE ("stats tests") {
    std::vector<uint64_t> empi{};
    std::vector<uint64_t> vi{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    std::vector<double> empd{};
    std::vector<double>
        vd{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};


    SECTION ("percentiles") {
        REQUIRE (percentile_from_sorted (empd, 0.5) == std::nullopt);
        REQUIRE (percentile_from_sorted (vd, 0.5) == 50.0);
        REQUIRE (percentile_from_sorted (vd, 0.25).value() == 25.0);
        REQUIRE (percentile_from_sorted (vd, 0.75) == 75.0);
    }

    SECTION ("MAD") {
        REQUIRE (mad (empi) == std::nullopt);
        REQUIRE (mad (vi) == 30.0);
    }

    SECTION ("gaps") {
        REQUIRE (gaps (empi) == empi);
        REQUIRE (
            gaps (vi)
            == std::vector<uint64_t>{10, 10, 10, 10, 10, 10, 10, 10, 10, 10}
        );
    }

    SECTION ("obs within span") {
        REQUIRE (max_obs_within (empi, 30ULL) == std::nullopt);
        REQUIRE (
            max_obs_within (vi, 30ULL).value()
            == static_cast<double> (3) / static_cast<double> (vi.size())
        );
    }

    SECTION ("min span with at least %obs") {
        REQUIRE (min_span_containing (empi, 0.5) == std::nullopt);
        REQUIRE (min_span_containing (vi, 0.25).value() == 20);
    }
}

TEST_CASE ("chebyshev test") {
    std::vector<endpoints1D<uint64_t>>
        v{{1, 3}, {1, 4}, {2, 6}, {7, 10}, {10, 12}};
    REQUIRE (ucheb (v[0], v[1]) == 1ULL);
    REQUIRE (ucheb (v[1], v[7]) == 4ULL);
    REQUIRE (ucheb (v[0], v[0]) == 0ULL);

    // check mst weights correct
    // where input order must not matter
    // and output order doesn't matter
    std::random_device rd;
    std::mt19937       g (rd());
    std::shuffle (v.begin(), v.end(), g);
    const auto              gaps = mst_cheb_dists (v);
    std::multiset<uint64_t> gapset (begin (gaps), end (gaps));
    REQUIRE (gapset == std::multiset<uint64_t>{3, 5, 2, 1});
}


// TODO adversarial tests
