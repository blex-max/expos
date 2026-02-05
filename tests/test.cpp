#include <catch2/catch_test_macros.hpp>
#include <cstdint>

#include "lib-stats/spatial.hpp"
#include "lib-stats/summary.hpp"
#include "lib-stats/string.hpp"

// TODO adversarial tests


TEST_CASE ("percentile") {
    using namespace summary;

    std::vector<uint64_t> empi{};
    std::vector<uint64_t> vi{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    std::vector<double> empd{};
    std::vector<double>
        vd{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};

    SECTION ("percentiles") {
        REQUIRE (percentile(empd, 0.5) == std::nullopt);
        REQUIRE (percentile(vd, 0.5) == 50.0);
        REQUIRE (percentile(vd, 0.25) == 25.0);
        REQUIRE (percentile(vd, 0.75) == 75.0);
    }
}

// TODO test with known binary sequences
TEST_CASE ("entropy") {
    using namespace string_stats;

    std::string t_low ("AAAAAAAAAA");
    std::string t_mid ("ACGTACGTAC");     // repeated phrase
    std::string t_hi ("ACAGTCAGGT");      // more disordered

    const double tsz = static_cast<double> (t_low.size());
    const auto exp_cmp = (2 * log2 (tsz)) / (tsz);
    const auto lz_low = entropy_lz76 (t_low);

    REQUIRE (
        std::fabs (lz_low - exp_cmp) < 1e-6
    );     // within float tolerance

    const auto lz_mid = entropy_lz76 (t_mid);
    const auto lz_hi  = entropy_lz76 (t_hi);

    // approximate test
    REQUIRE (lz_low < lz_mid);
    REQUIRE (lz_mid < lz_hi);
}


TEST_CASE ("periodic rle") {
    using namespace string_stats;

    {
        std::string runs ("CTCTCTAAA");
        std::vector<size_t> expected_res {6, 3};
        REQUIRE (periodic_rle(runs, 6) == expected_res);
    }

    {
        std::string runs ("CTCTCTAAAAAAAATAAATAAATC");
        std::vector<size_t> expected_res {6, 8, 8, 1, 1};
        REQUIRE (periodic_rle(runs, 6) == expected_res);
    }
}


TEST_CASE ("Ripley's K") {
    using namespace spatial;
    
    const double W = 100.0;
    const std::vector<uint64_t> obs{1, 2, 3, 95, 96, 97};
    const auto point_intensity = (static_cast<double> (obs.size()) / W);
    const auto t = 5;  // search radius
    const double exp_k = 2.0 * (1 / point_intensity); // given obs, we expect a mean of 2 points within t

    const auto pwds = PairMatrix::from_sample(
        obs,
        [] (const auto &a, const auto &b)
        { return (a > b) ? (a - b) : (b - a); }
    );
    const auto res = ripley_k(*pwds, t, point_intensity);

    CAPTURE (pwds->get1D());
    CAPTURE (res);
    CAPTURE (exp_k);
    REQUIRE (fabs(res - exp_k) < 1e-6);
}
