#include <catch2/catch_test_macros.hpp>

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

    std::string tlow ("AAAAAAAAAA");
    std::string tmid ("ACGTACGTAC");     // repeated phrase
    std::string thi ("ACAGTCAGGT");      // more disordered

    const double tsz = static_cast<double> (tlow.size());
    const auto exp_cmp = (2 * log2 (tsz)) / (tsz);
    const auto kclow = entropy_lz76 (tlow);

    REQUIRE (
        std::fabs (kclow - exp_cmp) < 1e-6
    );     // within float tolerance

    const auto kcmid = entropy_lz76 (tmid);
    const auto kchi  = entropy_lz76 (thi);

    // approximate test
    REQUIRE (kclow < kcmid);
    REQUIRE (kcmid < kchi);
}


TEST_CASE ("rle") {
    using namespace string_stats;

    std::string runs ("AACCCCCCGTTTT");  // 2, 6, 1, 4
    std::vector<size_t> expected_res {2, 6, 1, 4};

    REQUIRE (rle(runs) == expected_res);
}
