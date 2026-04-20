#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cstdint>

#include "lib-stats/monte-carlo.hpp"
#include "lib-stats/spatial.hpp"
#include "lib-stats/summary.hpp"
#include "lib-stats/string.hpp"


TEST_CASE ("percentile") {
    using namespace summary;

    std::vector<uint64_t> empi{};
    std::vector<uint64_t> vi{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    std::vector<double> empd{};
    std::vector<double>
        vd{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};

    SECTION ("percentiles") {
        REQUIRE (percentile(empd, 0.5) == std::nullopt);
        REQUIRE (std::fabs(*percentile(vd, 0.5)  - 50.0) < 1e-9);
        REQUIRE (std::fabs(*percentile(vd, 0.25) - 25.0) < 1e-9);
        REQUIRE (std::fabs(*percentile(vd, 0.75) - 75.0) < 1e-9);
    }

    SECTION ("single element") {
        std::vector<double> single{42.0};
        REQUIRE (percentile(single, 0.5).has_value());
        REQUIRE (std::fabs(*percentile(single, 0.5) - 42.0) < 1e-9);
    }

    SECTION ("duplicates") {
        std::vector<double> dups{5.0, 5.0, 5.0, 5.0, 5.0};
        REQUIRE (std::fabs(*percentile(dups, 0.25) - 5.0) < 1e-9);
        REQUIRE (std::fabs(*percentile(dups, 0.5)  - 5.0) < 1e-9);
        REQUIRE (std::fabs(*percentile(dups, 0.75) - 5.0) < 1e-9);
    }

    SECTION ("monotone") {
        REQUIRE (*percentile(vd, 0.25) < *percentile(vd, 0.5));
        REQUIRE (*percentile(vd, 0.5) < *percentile(vd, 0.75));
    }
}


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

    // approximate ordering test
    REQUIRE (lz_low < lz_mid);
    REQUIRE (lz_mid < lz_hi);
}


TEST_CASE ("entropy ground-truth") {
    using namespace string_stats;
    // LZ76 phrase counts verified by tracing entropy_lz76 source.
    // All n=4: entropy = (phrases * log2(4)) / 4 = phrases * 0.5
    REQUIRE (std::fabs (entropy_lz76 ("AAAA") - 1.0) < 1e-9);   // 2 phrases
    REQUIRE (std::fabs (entropy_lz76 ("ABAB") - 1.5) < 1e-9);   // 3 phrases
    REQUIRE (std::fabs (entropy_lz76 ("ABCD") - 2.0) < 1e-9);   // 4 phrases
}


TEST_CASE ("entropy edge cases") {
    using namespace string_stats;
    REQUIRE (std::fabs (entropy_lz76 ("")) < 1e-9);           // empty string early return
    REQUIRE (std::fabs (entropy_lz76 ("A") - 1.0) < 1e-9);   // single-char early return
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


TEST_CASE ("periodic rle edge cases") {
    using namespace string_stats;
    REQUIRE (periodic_rle("", 6).empty());
    REQUIRE (periodic_rle("A", 6) == std::vector<size_t>{1});
    REQUIRE (periodic_rle("AABB", 100) == std::vector<size_t>{2, 2});
}


TEST_CASE ("Ripley's K") {
    using namespace spatial;

    const double W = 100.0;
    const std::vector<uint64_t> obs{1, 2, 3, 95, 96, 97};
    const auto point_intensity = (static_cast<double> (obs.size()) / W);
    const auto t = 5;  // search radius
    const double exp_k = 2.0 * (1 / point_intensity);

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


TEST_CASE ("subsample_wo_replace") {
    using namespace monte_carlo;

    std::mt19937 rng(42);
    const std::vector<int> obs{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    SECTION ("correct size") {
        auto s = subsample_wo_replace(obs, 5, rng);
        REQUIRE (s.size() == 5);
    }

    SECTION ("no duplicates") {
        auto s = subsample_wo_replace(obs, 5, rng);
        auto sorted = s;
        std::sort(sorted.begin(), sorted.end());
        for (size_t i = 1; i < sorted.size(); ++i)
            REQUIRE (sorted[i] != sorted[i - 1]);
    }

    SECTION ("all elements from input") {
        auto s = subsample_wo_replace(obs, 5, rng);
        for (auto v : s)
            REQUIRE (std::find(obs.begin(), obs.end(), v) != obs.end());
    }

    SECTION ("n equals one") {
        auto s = subsample_wo_replace(obs, 1, rng);
        REQUIRE (s.size() == 1);
        REQUIRE (std::find(obs.begin(), obs.end(), s[0]) != obs.end());
    }

    SECTION ("n equals obs.size() minus one") {
        auto s = subsample_wo_replace(obs, obs.size() - 1, rng);
        REQUIRE (s.size() == obs.size() - 1);
        auto sorted = s;
        std::sort(sorted.begin(), sorted.end());
        for (size_t i = 1; i < sorted.size(); ++i)
            REQUIRE (sorted[i] != sorted[i - 1]);
    }
}


TEST_CASE ("sim_to_bg") {
    using namespace monte_carlo;

    // drawfn always returns the same 5-element vector; statfn = mean (= 3.0)
    auto drawfn = []() { return std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0}; };
    auto statfn = [](const std::vector<double>& v) {
        double s = 0;
        for (auto x : v) s += x;
        return s / static_cast<double>(v.size());
    };
    auto efffn = [](double ev, const std::vector<double>& sv) {
        return log2_effsz(ev, sv);
    };

    SECTION ("result fields populated") {
        auto r = sim_to_bg(3.0, drawfn, statfn, efffn, 10);
        REQUIRE (r.pval.has_value());
        REQUIRE (r.eff_sz.has_value());
    }

    SECTION ("extreme high observed gives minimum p-value") {
        const size_t nsim = 10;
        auto r = sim_to_bg(1000.0, drawfn, statfn, efffn, nsim);
        const double min_pval = 2.0 / static_cast<double>(nsim + 1);
        REQUIRE (std::fabs(*r.pval - min_pval) < 1e-9);
    }

    SECTION ("extreme low observed gives minimum p-value") {
        const size_t nsim = 10;
        auto r = sim_to_bg(0.0, drawfn, statfn, efffn, nsim);
        const double min_pval = 2.0 / static_cast<double>(nsim + 1);
        REQUIRE (std::fabs(*r.pval - min_pval) < 1e-9);
    }

    SECTION ("effect size positive when observed exceeds background") {
        auto r = sim_to_bg(100.0, drawfn, statfn, efffn, 10);
        REQUIRE (*r.eff_sz > 0.0);
    }

    SECTION ("effect size near zero when observed equals background mean") {
        auto r = sim_to_bg(3.0, drawfn, statfn, efffn, 10);
        REQUIRE (std::fabs(*r.eff_sz) < 1e-9);
    }
}
