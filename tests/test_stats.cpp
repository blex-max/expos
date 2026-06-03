#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cstdint>
#include <random>

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


TEST_CASE ("string entropy") {
    using namespace string_stats;

    SECTION ("expected values") {
        // LZ76 phrase counts verified by tracing entropy_lz76 source.
        // All n=4: entropy = (phrases * log2(4)) / 4 = phrases * 0.5
        REQUIRE (std::fabs (entropy_lz76 ("AAAA") - 1.0) < 1e-9);   // 2 phrases
        REQUIRE (std::fabs (entropy_lz76 ("ABAB") - 1.5) < 1e-9);   // 3 phrases
        REQUIRE (std::fabs (entropy_lz76 ("ABCD") - 2.0) < 1e-9);   // 4 phrases
    }

    SECTION ("edge cases") {
        using namespace string_stats;
        REQUIRE (std::fabs (entropy_lz76 ("")) < 1e-9);           // empty string early return
        REQUIRE (std::fabs (entropy_lz76 ("A") - 1.0) < 1e-9);   // single-char early return
    }
}


TEST_CASE ("periodic rle") {
    using namespace string_stats;

    SECTION ("expected values") {
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

    SECTION ("edge cases") {
        REQUIRE (periodic_rle("", 6).empty());
        REQUIRE (periodic_rle("A", 6) == std::vector<size_t>{1});
        REQUIRE (periodic_rle("AABB", 100) == std::vector<size_t>{2, 2});
    }
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

    std::mt19937_64 rng(24601);
    std::uniform_real_distribution<double> dist(1.0, 6.0);

    auto drawfn = [&]() { return std::vector<double> {
        dist(rng),
        dist(rng),
        dist(rng),
        dist(rng),
        dist(rng)
        }; };
    auto statfn = [](const std::vector<double>& v) {
        double s = 0;
        for (auto x : v) s += x;
        return s / static_cast<double>(v.size());
    };
    auto efffn = ses;
    size_t nsim = 100;

    SECTION ("result fields populated") {
        auto r = sim_to_bg(3.0, drawfn, statfn, efffn, nsim);
        REQUIRE (r.pval.has_value());
        REQUIRE (r.eff_sz.has_value());
    }

    SECTION ("effect size positive when observed exceeds background") {
        auto r = sim_to_bg(100.0, drawfn, statfn, efffn, nsim);
        REQUIRE (*r.eff_sz > 0.0);
    }

    SECTION ("effect size negative when background exceeds observed") {
        auto r = sim_to_bg(0.1, drawfn, statfn, efffn, nsim);
        REQUIRE (*r.eff_sz < 0.0);
    }

    SECTION ("extreme high observed gives minimum p-value") {
        // since such a value could not be observed from the
        // simulation population
        auto r = sim_to_bg(1000.0, drawfn, statfn, efffn, nsim);
        // assuming pval calc is extremes + 1 / N + 1
        // NOTE: it would be better to factor this out such
        // that it may be independently tested
        // (ditto effect size)
        const double min_pval = 1.0 / static_cast<double>(nsim + 1);
        REQUIRE (std::fabs(*r.pval - min_pval) < 1e-9);
    }

    SECTION ("extreme low observed gives insignificant p-value") {
        // since one tailed test
        auto r = sim_to_bg(0.1, drawfn, statfn, efffn, nsim);
        REQUIRE (std::fabs(*r.pval - 1.0) < 1e-9);
    }
}
