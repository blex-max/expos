#include <catch2/catch_test_macros.hpp>

#include "stats.hpp"

TEST_CASE ("kolmogorov complexity") {
    std::string tlow ("AAAAAAAAAA");
    std::string tmid ("ACGTACGTAC");     // repeated phrase
    std::string thi ("ACAGTCAGGT");

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

    // fuzzy test
    REQUIRE (kclow < kcmid);
    REQUIRE (kcmid < kchi);
}

// TODO adversarial tests
