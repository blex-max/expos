#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <set>

#include "stats.hpp"

TEST_CASE ("stats tests") {
    std::vector<uint64_t> empi{};
    std::vector<uint64_t>
        vi{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    std::vector<double> empd{};
    std::vector<double> vd{
        0.0,
        10.0,
        20.0,
        30.0,
        40.0,
        50.0,
        60.0,
        70.0,
        80.0,
        90.0,
        100.0
    };

    SECTION ("gaps") {
        REQUIRE (gaps (empi) == empi);
        REQUIRE (
            gaps (vi)
            == std::vector<
                uint64_t>{10, 10, 10, 10, 10, 10, 10, 10, 10, 10}
        );
    }

    SECTION ("min span with at least %obs") {
        REQUIRE (min_span_containing (empi, 0.5) == std::nullopt);
        REQUIRE (min_span_containing (vi, 0.25).value() == 20);
    }
}

TEST_CASE ("chebyshev tests") {
    std::vector<line_seg<uint64_t>>
               v{{1, 3}, {1, 4}, {2, 6}, {7, 10}, {10, 12}};
    const auto cheb_self = ucheb (v[0], v[0]);
    const auto cheb_max  = ucheb (v[0], v[4]);

    REQUIRE (cheb_self == 0ULL);
    REQUIRE (cheb_max == 9ULL);

    const auto pwd = PairMatrix<uint64_t>::from_sample (
        v,
        ucheb<uint64_t>
    );
    REQUIRE (pwd.get (0, 4) == 9ULL);
    REQUIRE (pwd.get (1, 3) == pwd.get (3, 1));

    // check mst weights correct
    // where input order must not matter
    // and output order doesn't matter
    // std::random_device rd;
    // std::mt19937       g (rd());
    // std::shuffle (v.begin(), v.end(), g);
    // const auto              gaps = mst_cheb_dists (v);
    // std::multiset<uint64_t> gapset (begin (gaps), end (gaps));
    // REQUIRE (gapset == std::multiset<uint64_t>{3, 5, 2, 1});

    REQUIRE (
        pwd.min_span_containing (0.5).value() == 3     // half
    );
    REQUIRE (
        pwd.min_span_containing (1).value() == 9     // total
    );
}


TEST_CASE ("kolmogorov complexity") {
    std::string tlow ("AAAAAAAAAA");
    std::string tmid ("ACGTACGTAC");     // repeated phrase
    std::string thi ("ACAGTCAGGT");

    const double tsz     = static_cast<double> (tlow.size());
    const auto   exp_cmp = (2 * log2 (tsz)) / (tsz); // TODO verify that this is the correct expectation
    const auto   kclow   = nk_lz76 (tlow);

    REQUIRE (std::fabs (kclow - exp_cmp) < 1e-6);  // within float tolerance

    const auto kcmid = nk_lz76 (tmid);
    const auto kchi  = nk_lz76 (thi);

    // fuzzy test
    REQUIRE (kclow < kcmid);
    REQUIRE (kcmid < kchi);
}

// TODO adversarial tests
