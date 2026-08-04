// Unit tests for the pure statistical primitives in shared/stats.hpp.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include "shared/stats.hpp"

using Catch::Approx;

static std::size_t brute_pairs_1d (
    const std::vector<int32_t>& v, uint64_t r
)
{
  std::size_t count = 0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    for (std::size_t j = i + 1; j < v.size(); ++j) {
      const uint64_t d =
          v[i] > v[j] ? static_cast<uint64_t> (v[i] - v[j])
                      : static_cast<uint64_t> (v[j] - v[i]);
      if (d < r) {
        ++count;
      }
    }
  }
  return count;
}

TEST_CASE ("mean")
{
  REQUIRE_FALSE (mean ({}).has_value());
  REQUIRE (*mean ({1.0, 2.0, 3.0, 4.0}) == Approx (2.5));
  REQUIRE (*mean ({42.0}) == Approx (42.0));
}

TEST_CASE ("percentile")
{
  REQUIRE_FALSE (
      percentile (std::vector<uint32_t>{}, 0.5).has_value()
  );
  REQUIRE (
      *percentile (std::vector<uint32_t>{7}, 0.5) == Approx (7.0)
  );

  const std::vector<uint64_t> ladder{0,  10, 20, 30, 40, 50,
                                     60, 70, 80, 90, 100};
  REQUIRE (*percentile (ladder, 0.5) == Approx (50.0));
  REQUIRE (*percentile (ladder, 0.25) == Approx (25.0));
  REQUIRE (*percentile (ladder, 0.75) == Approx (75.0));
  REQUIRE (
      *percentile (std::vector<uint32_t>{1, 2, 3, 4}, 0.5) ==
      Approx (2.5)
  );
  // linear interpolation: rank 0.25 between 0 and 10 -> 2.5
  REQUIRE (
      *percentile (std::vector<uint32_t>{0, 10}, 0.25) ==
      Approx (2.5)
  );
  // floating-point sample + duplicates
  REQUIRE (
      *percentile (
          std::vector<double>{5.0, 5.0, 5.0, 5.0}, 0.25
      ) == Approx (5.0)
  );
}

TEST_CASE ("count_pairs_within_1d")
{
  REQUIRE (count_pairs_within_1d ({}, 5) == 0);
  REQUIRE (count_pairs_within_1d ({5}, 10) == 0);
  REQUIRE (count_pairs_within_1d ({1, 2, 3}, 2) == 2);
  REQUIRE (count_pairs_within_1d ({1, 2, 3}, 3) == 3);
  REQUIRE (
      count_pairs_within_1d ({7, 7, 7}, 1) == 3
  );  // dups, dist 0 < 1
  REQUIRE (
      count_pairs_within_1d ({1, 2, 3}, 0) == 0
  );  // radius 0 -> none

  SECTION ("randomised cross-check against a brute-force oracle")
  {
    std::mt19937 rng (123);
    std::uniform_int_distribution<int32_t> valDist (0, 50);
    for (int trial = 0; trial < 300; ++trial) {
      const std::size_t nObs = rng() % 40;
      std::vector<int32_t> vals;
      vals.reserve (nObs);
      for (std::size_t i = 0; i < nObs; ++i) {
        vals.push_back (valDist (rng));
      }
      const uint64_t radius = rng() % 20;
      REQUIRE (
          count_pairs_within_1d (vals, radius) ==
          brute_pairs_1d (vals, radius)
      );
    }
  }
}

TEST_CASE ("count_pairs_within (generic)")
{
  using Point = std::pair<int64_t, int64_t>;
  auto absDiff = [] (int64_t a, int64_t b) -> uint64_t {
    return a > b ? static_cast<uint64_t> (a - b)
                 : static_cast<uint64_t> (b - a);
  };
  auto manhattan =
      [&] (const Point& a, const Point& b) -> uint64_t {
    return absDiff (a.first, b.first) +
           absDiff (a.second, b.second);
  };
  const std::vector<Point> pts{{0, 0}, {1, 1}, {10, 10}};
  REQUIRE (
      count_pairs_within (pts, 3, manhattan) == 1
  );  // (0,0)-(1,1)=2
  REQUIRE (
      count_pairs_within (pts, 2, manhattan) == 0
  );  // 2 not < 2
  REQUIRE (count_pairs_within (pts, 100, manhattan) == 3);
  REQUIRE (
      count_pairs_within (std::vector<Point>{}, 5, manhattan) ==
      0
  );
}

TEST_CASE ("ripley_k")
{
  // Two tight clusters, radius 5: each cluster of 3 contributes 3
  // within-radius pairs (6 total); the clusters are ~92 apart so no
  // cross-cluster pairs. K = (1/intensity) * 2 * pairs / n.
  const std::vector<int32_t> obs{1, 2, 3, 95, 96, 97};
  const double intensity =
      static_cast<double> (obs.size()) / 100.0;
  const std::size_t pairs = count_pairs_within_1d (obs, 5);
  REQUIRE (pairs == 6);
  REQUIRE (
      ripley_k (pairs, obs.size(), intensity) ==
      Approx (2.0 * (1.0 / intensity))
  );
}

TEST_CASE ("monte_carlo_pvalue")
{
  auto identity = [] (double x) { return x; };

  SECTION ("exact against a deterministic null")
  {
    // null {2,4,4,4,5,5,7,9}: mean 5, population SD 2.
    const std::vector<double> seq{2, 4, 4, 4, 5, 5, 7, 9};
    auto run = [&] (double observed) {
      std::size_t idx = 0;
      auto draw = [&]() { return seq[idx++]; };
      return monte_carlo_pvalue (
          observed, draw, identity, seq.size()
      );
    };

    const auto above = run (7.0);
    REQUIRE (above.effectSize.has_value());
    REQUIRE (*above.effectSize == Approx (1.0));  // (7-5)/2
    REQUIRE (
        above.pValue == Approx (3.0 / 9.0)
    );  // #{>=7}=2 -> 3/9

    const auto below = run (3.0);
    REQUIRE (*below.effectSize == Approx (-1.0));  // (3-5)/2

    REQUIRE (
        run (1000.0).pValue == Approx (1.0 / 9.0)
    );  // #{>=}=0
    REQUIRE (run (0.0).pValue == Approx (1.0));  // all >= 0
  }

  SECTION ("zero-variance null yields no effect size")
  {
    auto draw = [] { return 3.0; };
    const auto res =
        monte_carlo_pvalue (3.0, draw, identity, 100);
    REQUIRE_FALSE (res.effectSize.has_value());
    REQUIRE (res.pValue == Approx (1.0));
  }

  SECTION ("deterministic under a fixed seed; p-value in (0,1]")
  {
    std::vector<int32_t> pop (200);
    for (std::size_t i = 0; i < pop.size(); ++i) {
      pop[i] = static_cast<int32_t> (i);
    }
    auto run = [&] (uint32_t seed) {
      std::mt19937 rng (seed);
      auto draw = [&]() {
        return subsample_wo_replace (pop, 20, rng);
      };
      auto stat = [] (const std::vector<int32_t>& s) {
        return static_cast<double> (
            count_pairs_within_1d (s, 5)
        );
      };
      return monte_carlo_pvalue (10.0, draw, stat, 500);
    };
    const auto a = run (99);
    const auto b = run (99);
    REQUIRE (a.pValue == Approx (b.pValue));
    REQUIRE (
        a.effectSize.has_value() == b.effectSize.has_value()
    );
    if (a.effectSize && b.effectSize) {
      REQUIRE (*a.effectSize == Approx (*b.effectSize));
    }
    REQUIRE (a.pValue > 0.0);
    REQUIRE (a.pValue <= 1.0);
  }
}

TEST_CASE ("subsample_wo_replace")
{
  const std::vector<int> obs{10, 20, 30, 40, 50};

  SECTION ("size, membership, no replacement")
  {
    std::mt19937 rng (7);
    const auto s = subsample_wo_replace (obs, 3, rng);
    REQUIRE (s.size() == 3);
    const std::set<int> uniq (s.begin(), s.end());
    REQUIRE (uniq.size() == 3);
    for (const int x : s) {
      REQUIRE (
          std::find (obs.begin(), obs.end(), x) != obs.end()
      );
    }
  }

  SECTION ("deterministic under identical seed")
  {
    std::mt19937 rngA (7);
    std::mt19937 rngB (7);
    REQUIRE (
        subsample_wo_replace (obs, 3, rngA) ==
        subsample_wo_replace (obs, 3, rngB)
    );
  }

  SECTION ("n == nObs yields a full permutation")
  {
    std::mt19937 rng (1);
    const auto full =
        subsample_wo_replace (obs, obs.size(), rng);
    const std::set<int> uniq (full.begin(), full.end());
    REQUIRE (uniq.size() == obs.size());
  }

  SECTION ("n == 0 yields empty")
  {
    std::mt19937 rng (1);
    REQUIRE (subsample_wo_replace (obs, 0, rng).empty());
  }
}

TEST_CASE ("entropy_lz76")
{
  // exact LZ76 values (n=4: entropy = phrases * log2(4)/4 = phrases * 0.5)
  REQUIRE (entropy_lz76 ("AAAA") == Approx (1.0));  // 2 phrases
  REQUIRE (entropy_lz76 ("ABAB") == Approx (1.5));  // 3 phrases
  REQUIRE (entropy_lz76 ("ABCD") == Approx (2.0));  // 4 phrases

  REQUIRE (entropy_lz76 ("") == Approx (0.0));
  REQUIRE (entropy_lz76 ("A") == Approx (1.0));

  // repetitive is less complex per char than diverse
  REQUIRE (
      entropy_lz76 ("AAAAAAAAAAAAAAAA") <
      entropy_lz76 ("ACGTACGATCGGATCA")
  );
}
