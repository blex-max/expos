#include <catch2/catch_test_macros.hpp>

#include "stats.hpp"

TEST_CASE ("median even") {
    std::vector<double> v{1.0, 2.0, 3.0, 4.0};
    REQUIRE (median_from_sorted (v) == 2.5);
}

TEST_CASE ("comprehensive MAD test, median odd") {
    std::vector<int64_t> v{1, 2, 3, 4, 5};
    const auto           vd{as_double (v)};
    const auto           med  = median_from_sorted (vd);
    const auto           devv = devs (vd, med);

    REQUIRE (vd == std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0});
    REQUIRE (med == 3.0);
    REQUIRE (devv == std::vector<double>{2.0, 1.0, 0.0, 1.0, 2.0});
    REQUIRE (mad (v) == 1.0);
};
