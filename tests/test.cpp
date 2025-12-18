#include <catch2/catch_test_macros.hpp>

#include "stats.hpp"


TEST_CASE ("percentile") {
    std::vector<uint64_t> empi{};
    std::vector<uint64_t> vi{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    std::vector<double> empd{};
    std::vector<double>
        vd{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};

    SECTION ("percentiles") {
        REQUIRE (percentile(empd, 0.5) == std::nullopt);
        REQUIRE (percentile(vd, 0.5) == 50.0);
        REQUIRE (percentile(vd, 0.25).value() == 25.0);
        REQUIRE (percentile(vd, 0.75) == 75.0);
    }
}

TEST_CASE ("kolmogorov complexity") {
    std::string tlow ("AAAAAAAAAA");
    std::string tmid ("ACGTACGTAC");     // repeated phrase
    std::string thi ("ACAGTCAGGT");      // more disordered

    const double tsz = static_cast<double> (tlow.size());
    const auto
        exp_cmp = (2 * log2 (tsz))
                  / (tsz);     // TODO verify that this is the correct expectation
    const auto kclow = nk_lz76 (tlow);

    REQUIRE (
        std::fabs (kclow - exp_cmp) < 1e-6
    );     // within float tolerance

    const auto kcmid = nk_lz76 (tmid);
    const auto kchi  = nk_lz76 (thi);

    // approximate test
    REQUIRE (kclow < kcmid);
    REQUIRE (kcmid < kchi);
}

// TODO adversarial tests
