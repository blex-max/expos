// Unit tests for the variant statistics layer: the clustering, overlap and
// Monte-Carlo primitives in variant_stats.hpp, then the compute functions in
// compute_info_field.hpp built on them. The latter are exercised through
// the public expos_field_registry(), so they stay internal.

#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <random>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "expos/compute_info_field.hpp"
#include "expos/encode_info_field.hpp"
#include "expos/pileup_features.hpp"
#include "expos/skip.hpp"
#include "expos/variant_stats.hpp"

using Catch::Approx;

// --- primitives (variant_stats.hpp) --- //

namespace {

std::size_t brute_pairs_1d (
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

}  // namespace

namespace {

// count_pairs_within_1d takes a mutable span and sorts through it, so each
// case needs its own storage. Returns the count for a fresh copy of `v`.
std::size_t pairs_of (std::vector<int32_t> v, uint64_t r)
{
  return count_pairs_within_1d (v, r);
}

}  // namespace

TEST_CASE ("count_pairs_within_1d")
{
  REQUIRE (pairs_of ({}, 5) == 0);
  REQUIRE (pairs_of ({5}, 10) == 0);
  REQUIRE (pairs_of ({1, 2, 3}, 2) == 2);
  REQUIRE (pairs_of ({1, 2, 3}, 3) == 3);
  REQUIRE (pairs_of ({7, 7, 7}, 1) == 3);  // dups, dist 0 < 1
  REQUIRE (pairs_of ({1, 2, 3}, 0) == 0);  // radius 0 -> none

  SECTION ("sorts its argument in place")
  {
    // The in-place sort is what lets the Monte-Carlo draws reuse one buffer;
    // callers whose data must survive pass a copy.
    std::vector<int32_t> v{30, 10, 20};
    count_pairs_within_1d (v, 5);
    REQUIRE (v == std::vector<int32_t>{10, 20, 30});
  }

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
      // Brute force first: the call below reorders `vals`.
      const std::size_t expected = brute_pairs_1d (vals, radius);
      REQUIRE (count_pairs_within_1d (vals, radius) == expected);
    }
  }
}

TEST_CASE ("run_monte_carlo")
{
  auto identity = [] (double x) { return x; };

  SECTION ("exact against a deterministic null")
  {
    // null {2,4,4,4,5,5,7,9}: mean 5, population SD 2.
    const std::vector<double> seq{2, 4, 4, 4, 5, 5, 7, 9};
    auto run = [&] (double observed) {
      std::size_t idx = 0;
      auto draw = [&]() { return seq[idx++]; };
      return run_monte_carlo (
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
    const auto res = run_monte_carlo (3.0, draw, identity, 100);
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
      SubsampleScratch<int32_t> scratch;
      auto draw = [&]() {
        return subsample_wo_replace (pop, 20, rng, scratch);
      };
      auto stat = [] (std::span<int32_t> s) {
        return static_cast<double> (
            count_pairs_within_1d (s, 5)
        );
      };
      return run_monte_carlo (10.0, draw, stat, 500);
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
    SubsampleScratch<int> scratch;
    const auto s = subsample_wo_replace (obs, 3, rng, scratch);
    REQUIRE (s.size() == 3);
    const std::set<int> uniq (s.begin(), s.end());
    REQUIRE (uniq.size() == 3);
    for (const int x : s) {
      REQUIRE (
          std::find (obs.begin(), obs.end(), x) != obs.end()
      );
    }
  }

  SECTION (
      "deterministic under identical seed and fresh scratch"
  )
  {
    auto draw = [&] (uint32_t seed) {
      std::mt19937 rng (seed);
      SubsampleScratch<int> scratch;
      const auto s = subsample_wo_replace (obs, 3, rng, scratch);
      return std::vector<int> (s.begin(), s.end());
    };
    REQUIRE (draw (7) == draw (7));
  }

  SECTION ("n == nObs yields a full permutation")
  {
    std::mt19937 rng (1);
    SubsampleScratch<int> scratch;
    const auto full =
        subsample_wo_replace (obs, obs.size(), rng, scratch);
    const std::set<int> uniq (full.begin(), full.end());
    REQUIRE (uniq.size() == obs.size());
  }

  SECTION ("n == 0 yields empty")
  {
    std::mt19937 rng (1);
    SubsampleScratch<int> scratch;
    REQUIRE (
        subsample_wo_replace (obs, 0, rng, scratch).empty()
    );
  }

  SECTION ("one scratch survives a changing population size")
  {
    // The index permutation is only valid for the population it was built
    // for, so a size change must rebuild it rather than resize it.
    std::mt19937 rng (3);
    SubsampleScratch<int32_t> scratch;
    for (const std::size_t nObs : {60U, 17U, 60U, 200U, 5U}) {
      std::vector<int32_t> pop (nObs);
      for (std::size_t i = 0; i < nObs; ++i) {
        pop[i] = static_cast<int32_t> (i);
      }
      const std::size_t n = nObs / 2;
      const auto s = subsample_wo_replace (pop, n, rng, scratch);
      REQUIRE (s.size() == n);
      const std::set<int32_t> uniq (s.begin(), s.end());
      REQUIRE (uniq.size() == n);  // still no replacement
      for (const int32_t v : s) {
        REQUIRE (v >= 0);
        REQUIRE (v < static_cast<int32_t> (nObs));
      }
    }
  }

  SECTION ("repeated draws from one scratch stay uniform")
  {
    // The scratch carries its permutation between draws rather than
    // re-initialising it, which is what makes a draw O(n) instead of
    // O(population). That is only sound if the carried permutation biases
    // nothing: draw k+1 applies a fresh independent partial Fisher-Yates and
    // then maps through a fixed bijection, and a bijection of a uniform
    // n-subset is a uniform n-subset. This is the test for that claim.
    constexpr std::size_t nObs = 40;
    constexpr std::size_t n = 8;
    constexpr std::size_t nDraws = 40000;

    std::vector<int32_t> pop (nObs);
    for (std::size_t i = 0; i < nObs; ++i) {
      pop[i] = static_cast<int32_t> (i);
    }

    std::mt19937 rng (20260806);
    SubsampleScratch<int32_t> scratch;
    std::vector<std::size_t> hits (nObs, 0);
    for (std::size_t d = 0; d < nDraws; ++d) {
      for (const int32_t v :
           subsample_wo_replace (pop, n, rng, scratch)) {
        ++hits[static_cast<std::size_t> (v)];
      }
    }

    const double expected = static_cast<double> (nDraws * n) /
                            static_cast<double> (nObs);
    double chi2 = 0.0;
    for (const std::size_t h : hits) {
      const double d = static_cast<double> (h) - expected;
      chi2 += (d * d) / expected;
    }
    // 39 df; the 0.1% upper tail is 72.05. A biased carry blows well past it.
    REQUIRE (chi2 < 72.05);
  }
}

// --- compute layer (compute_info_field.hpp) --- //

namespace {

const VariantStat& by_id (
    std::span<const VariantStat> stats, std::string_view id
)
{
  for (const auto& s : stats) {
    if (s.field.id == id) {
      return s;
    }
  }
  throw std::runtime_error ("stat not found");
}

constexpr std::string_view REF_PLACEHOLDER{};

}  // namespace

TEST_CASE ("variant_stats registry")
{
  const auto stats = expos_field_registry();
  REQUIRE (
      stats.size() == 4
  );  // QRK, TJAC, MLAS, RCMPLX — the full package
  REQUIRE (by_id (stats, "QRK").field.nValues == 2);
  REQUIRE (by_id (stats, "TJAC").field.nValues == 2);
  REQUIRE (by_id (stats, "MLAS").field.nValues == 2);
  REQUIRE (by_id (stats, "RCMPLX").field.nValues == 1);
}

TEST_CASE ("compute QRK (query-position clustering)")
{
  const auto qrk = by_id (expos_field_registry(), "QRK");
  std::mt19937 rng (1);

  SECTION ("missing with insufficient support")
  {
    PileupFeatures supporting;
    supporting.qPos = {10};
    PileupFeatures all;
    all.qPos = {0, 5, 10, 15, 20, 25};
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = qrk.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (r.error() == StatSkipReason::insufficientSupport);
  }

  SECTION ("missing with insufficient background")
  {
    PileupFeatures supporting;
    supporting.qPos = {10, 11, 12};
    PileupFeatures all;
    all.qPos = {
        10, 11, 12, 40, 90
    };  // 5: under 2*3, and under the floor
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = qrk.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (
        r.error() == StatSkipReason::insufficientBackground
    );
  }

  SECTION ("suppressed when read lengths are heterogeneous")
  {
    PileupFeatures supporting;
    supporting.qPos = {50, 51, 52, 53};
    PileupFeatures all;
    for (int32_t i = 0; i < 200; ++i) {
      all.qPos.push_back (i);
    }
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{
        mc, StatSkipReason::heterogeneousReadLength
    };
    const auto r = qrk.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (
        r.error() == StatSkipReason::heterogeneousReadLength
    );
  }

  // The ratio alone would admit this: 6 background against 2 supporting is
  // comfortably over 2x. But 6 gives C(6,2) = 15 distinct subsamples to
  // build a null from, so the absolute floor has to catch it.
  SECTION (
      "a background over the ratio but under the floor is "
      "refused"
  )
  {
    PileupFeatures supporting;
    supporting.qPos = {50, 51};
    PileupFeatures all;
    all.qPos = {10, 20, 30,
                40, 50, 51};  // 6 >= 2*2, but < MIN_BACKGROUND
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = qrk.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (
        r.error() == StatSkipReason::insufficientBackground
    );
  }

  // The guard fails closed on a primary with too few reads to measure the
  // spread. That is not the same claim as reads of uneven length, and the
  // token has to say which one happened.
  SECTION (
      "unverifiable read lengths are reported as unverified"
  )
  {
    PileupFeatures supporting;
    supporting.qPos = {50, 51, 52, 53};
    PileupFeatures all;
    for (int32_t i = 0; i < 200; ++i) {
      all.qPos.push_back (i);
    }
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{
        mc, StatSkipReason::readLengthUnverified
    };
    const auto r = qrk.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (r.error() == StatSkipReason::readLengthUnverified);
  }

  // A locus too thin to test is reported as thin, whatever the read-length
  // guard concluded: the size check runs first precisely so that a site
  // with no reads never gets a verdict on its read-length distribution.
  SECTION (
      "too little support outranks a read-length suppression"
  )
  {
    const PileupFeatures
        supporting;  // no supporting reads at all
    PileupFeatures all;
    for (int32_t i = 0; i < 200; ++i) {
      all.qPos.push_back (i);
    }
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{
        mc, StatSkipReason::heterogeneousReadLength
    };
    const auto r = qrk.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (r.error() == StatSkipReason::insufficientSupport);
  }

  SECTION ("tight support vs spread background is extreme")
  {
    PileupFeatures supporting;
    supporting.qPos = {50, 51, 52, 53};  // all within radius 5
    PileupFeatures all;
    for (int32_t i = 0; i < 200; ++i) {
      all.qPos.push_back (i);  // spread 0..199
    }
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = qrk.compute (in, ctx);
    REQUIRE (r.has_value());
    REQUIRE (r->size() == 2);  // matches the registry's nValues
    REQUIRE ((*r)[0].value.has_value());  // effect size
    REQUIRE ((*r)[1].value.has_value());  // p-value
    REQUIRE (
        *(*r)[0].value > 0.0
    );  // more clustered than the null
    REQUIRE (*(*r)[1].value < 0.05);  // significant
    REQUIRE (*(*r)[1].value > 0.0);
  }

  SECTION ("whole statistic skipped when the null has no spread")
  {
    // Every background position is identical, so every draw of 4 yields the
    // same C(4,2) = 6 pairs within radius 5 and the null is degenerate. The
    // z-score is then undefined and the p-value is pinned by the degeneracy
    // rather than the data, so neither subfield is reportable.
    PileupFeatures supporting;
    supporting.qPos = {50, 50, 50, 50};
    PileupFeatures all;
    for (int32_t i = 0; i < 100; ++i) {
      all.qPos.push_back (50);
    }
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = qrk.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (r.error() == StatSkipReason::zeroVariance);
  }

  SECTION ("deterministic under a fixed seed")
  {
    PileupFeatures supporting;
    supporting.qPos = {50, 51, 52, 53};
    PileupFeatures all;
    for (int32_t i = 0; i < 200; ++i) {
      all.qPos.push_back (i);
    }
    std::mt19937 rngA (7);
    std::mt19937 rngB (7);
    VariantStatInputs inA{supporting, all, REF_PLACEHOLDER};
    McState mcA{std::move (rngA), {}, {}};
    const StatContext ctxA{mcA, std::nullopt};
    VariantStatInputs inB{supporting, all, REF_PLACEHOLDER};
    McState mcB{std::move (rngB), {}, {}};
    const StatContext ctxB{mcB, std::nullopt};
    const auto a = qrk.compute (inA, ctxA);
    const auto b = qrk.compute (inB, ctxB);
    REQUIRE (a.has_value());
    REQUIRE (b.has_value());
    REQUIRE ((*a)[0].value.has_value());
    REQUIRE (*(*a)[0].value == Approx (*(*b)[0].value));
    REQUIRE (*(*a)[1].value == Approx (*(*b)[1].value));
  }
}

namespace {

// Templates of three very different lengths, all covering ~1000, positions
// jittered so no two coincide. Under a min(len) denominator almost every
// pair here is fully nested and scores 1.0; under Jaccard they do not.
std::vector<TemplateEndpoints> mixed_length_background()
{
  std::vector<TemplateEndpoints> out;
  for (int64_t i = 0; i < 30; ++i) {
    out.push_back ({1000 - 10 - i, 1000 + 10 - i});
    out.push_back ({1000 - 100 - i, 1000 + 100 - i});
    out.push_back ({1000 - 200 - i, 1000 + 200 - i});
  }
  return out;
}

}  // namespace

TEST_CASE ("compute TJAC (graded pairwise template overlap)")
{
  const auto tjac = by_id (expos_field_registry(), "TJAC");

  SECTION ("missing with insufficient support")
  {
    std::mt19937 rng (1);
    PileupFeatures supporting;
    supporting.endpoints = {{100, 300}};
    PileupFeatures all;
    all.endpoints = {
        {100, 300}, {110, 320}, {500, 700}, {900, 1100}
    };
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = tjac.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (r.error() == StatSkipReason::insufficientSupport);
  }

  SECTION ("missing with insufficient background")
  {
    std::mt19937 rng (1);
    PileupFeatures supporting;
    supporting.endpoints = {{100, 300}, {101, 301}, {102, 302}};
    PileupFeatures all;
    all.endpoints = {
        {100, 300},
        {101, 301},
        {102, 302},
        {500, 700},
        {900, 1100}
    };  // 5: under 2*3, and under the floor
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = tjac.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (
        r.error() == StatSkipReason::insufficientBackground
    );
  }

  // TJAC sizes on templates, so the floor is counted in templates here --
  // the same guard applied to a different population.
  SECTION (
      "a background over the ratio but under the floor is "
      "refused"
  )
  {
    std::mt19937 rng (1);
    PileupFeatures supporting;
    supporting.endpoints = {{100, 300}, {101, 301}};
    PileupFeatures all;
    for (int64_t i = 0; i < 6;
         ++i) {  // 6 >= 2*2, but < MIN_BACKGROUND
      all.endpoints.push_back ({i, i + 200});
    }
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = tjac.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (
        r.error() == StatSkipReason::insufficientBackground
    );
  }

  SECTION ("clustered endpoints vs spread background is extreme")
  {
    std::mt19937 rng (2);
    PileupFeatures supporting;
    supporting.endpoints = {
        {100, 300}, {101, 301}, {102, 302}, {100, 300}
    };  // near-identical templates
    PileupFeatures all;
    for (int64_t i = 0; i < 100; ++i) {
      // spread over 0..99, so draws have a spread of pairwise overlap and
      // the null has non-zero variance.
      all.endpoints.push_back ({i, i + 200});
    }
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = tjac.compute (in, ctx);
    REQUIRE (r.has_value());
    REQUIRE ((*r)[0].value.has_value());
    REQUIRE (*(*r)[0].value > 0.0);
    REQUIRE (*(*r)[1].value < 0.05);
  }

  SECTION (
      "a support set typical of the background does not fire"
  )
  {
    std::mt19937 rng (2);
    PileupFeatures all;
    for (int64_t i = 0; i < 100; ++i) {
      all.endpoints.push_back ({i, i + 200});
    }
    PileupFeatures supporting;  // drawn from the same population
    supporting.endpoints = {
        all.endpoints[7], all.endpoints[31], all.endpoints[62],
        all.endpoints[88]
    };
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = tjac.compute (in, ctx);
    REQUIRE (r.has_value());
    REQUIRE ((*r)[0].value.has_value());
    REQUIRE (
        *(*r)[0].value < 0.0
    );  // less overlap than a typical draw
    REQUIRE (*(*r)[1].value > 0.1);
  }

  SECTION ("length disparity alone does not fire")
  {
    // Nested templates of wildly different lengths, against a background of
    // the same length mixture. A min(len) denominator would score nearly
    // every pair here at 1.0; Jaccard must not read that as coincidence.
    std::mt19937 rng (2);
    PileupFeatures supporting;
    supporting.endpoints = {
        {990, 1010}, {900, 1100}, {995, 1015}, {800, 1200}
    };
    PileupFeatures all;
    all.endpoints = mixed_length_background();
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = tjac.compute (in, ctx);
    REQUIRE (r.has_value());
    REQUIRE ((*r)[0].value.has_value());
    REQUIRE (*(*r)[1].value > 0.1);
  }

  SECTION (
      "coincident templates fire against a length-disparate "
      "background"
  )
  {
    // The regression guard for the Jaccard-over-min(len) decision. The
    // background here is mostly nested pairs, so a min(len) denominator
    // leaves the null sitting at ~94% of its ceiling and this genuinely
    // coincident support set is invisible to it (measured z 0.62, p 0.16).
    // Jaccard keeps the null near 39% of ceiling, so the same set stands
    // out sharply.
    std::mt19937 rng (2);
    PileupFeatures supporting;
    supporting.endpoints = {
        {950, 1050}, {950, 1050}, {951, 1051}, {950, 1050}
    };
    PileupFeatures all;
    all.endpoints = mixed_length_background();
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = tjac.compute (in, ctx);
    REQUIRE (r.has_value());
    REQUIRE ((*r)[0].value.has_value());
    REQUIRE (*(*r)[0].value > 2.0);
    REQUIRE (*(*r)[1].value < 0.05);
  }

  SECTION ("whole statistic skipped when the null has no spread")
  {
    // Every background template is identical, so every pair in every draw
    // scores Jaccard 1.0 and every draw sums to the same 6.0. Degenerate null,
    // so neither the z-score nor its p-value is reportable.
    std::mt19937 rng (2);
    PileupFeatures supporting;
    supporting.endpoints = {
        {100, 300}, {100, 300}, {100, 300}, {100, 300}
    };
    PileupFeatures all;
    for (int64_t i = 0; i < 100; ++i) {
      all.endpoints.push_back ({100, 300});
    }
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = tjac.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (r.error() == StatSkipReason::zeroVariance);
  }

  SECTION (
      "not suppressed by heterogeneous read lengths (unlike QRK)"
  )
  {
    std::mt19937 rng (2);
    PileupFeatures supporting;
    supporting.endpoints = {
        {100, 300}, {101, 301}, {102, 302}, {100, 300}
    };
    PileupFeatures all;
    for (int64_t i = 0; i < 100; ++i) {
      all.endpoints.push_back ({i, i + 200});
    }
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{
        mc, StatSkipReason::heterogeneousReadLength
    };
    const auto r = tjac.compute (in, ctx);
    REQUIRE (r.has_value());
    REQUIRE (
        (*r)[0].value.has_value()
    );  // computed despite heterogeneity
  }
}

TEST_CASE ("compute MLAS (median normalised alignment score)")
{
  const auto mlas = by_id (expos_field_registry(), "MLAS");
  std::mt19937 rng (1);

  SECTION ("medians of supporting and all")
  {
    PileupFeatures supporting;
    supporting.normalisedAs = {0.5, 0.7, 0.9};  // median 0.7
    PileupFeatures all;
    all.normalisedAs = {0.2, 0.4, 0.6, 0.8};  // median 0.5
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = mlas.compute (in, ctx);
    REQUIRE (r.has_value());
    REQUIRE (r->size() == 2);
    REQUIRE (*(*r)[0].value == Approx (0.7));
    REQUIRE (*(*r)[1].value == Approx (0.5));
  }

  SECTION ("missing when the supporting group is empty")
  {
    // MLAS's subfields are independent summaries, not one inference, so this
    // is a per-subfield skip and not a whole-statistic one: the background
    // median stays reportable.
    PileupFeatures supporting;  // no reads
    PileupFeatures all;
    all.normalisedAs = {0.3, 0.6};
    VariantStatInputs in{supporting, all, REF_PLACEHOLDER};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = mlas.compute (in, ctx);
    REQUIRE (r.has_value());
    REQUIRE_FALSE ((*r)[0].value.has_value());
    REQUIRE ((*r)[0].reason == StatSkipReason::noSupport);
    REQUIRE ((*r)[1].value.has_value());
  }
}

TEST_CASE ("compute RCMPLX (reference complexity)")
{
  const auto rcmplx = by_id (expos_field_registry(), "RCMPLX");
  std::mt19937 rng (1);
  const PileupFeatures empty;

  SECTION ("missing when the slice is shorter than the window")
  {
    const std::string ref (99, 'A');
    VariantStatInputs in{empty, empty, std::string_view (ref)};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = rcmplx.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (r.error() == StatSkipReason::referenceTooShort);
  }

  SECTION ("missing when the slice contains N")
  {
    std::string ref (150, 'A');
    ref[75] = 'N';
    VariantStatInputs in{empty, empty, std::string_view (ref)};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = rcmplx.compute (in, ctx);
    REQUIRE_FALSE (r.has_value());
    REQUIRE (r.error() == StatSkipReason::referenceHasN);
  }

  SECTION ("homopolymer window has low, known complexity")
  {
    const std::string ref (
        100, 'A'
    );  // exactly one 100-base window
    VariantStatInputs in{empty, empty, std::string_view (ref)};
    McState mc{std::move (rng), {}, {}};
    const StatContext ctx{mc, std::nullopt};
    const auto r = rcmplx.compute (in, ctx);
    REQUIRE (r.has_value());
    REQUIRE ((*r)[0].value.has_value());
    // entropy_lz76 of 100 identical chars = 2*log2(100)/100 ~= 0.1329
    REQUIRE (*(*r)[0].value == Approx (0.13287712));
  }
}

TEST_CASE ("stat_skips deduplicates within a statistic")
{
  const VariantStatField field{"QRK", "", 2};

  SECTION ("both subfields share a reason -> one skip")
  {
    const std::vector<StatValue> bothMissing{
        StatValue{
            std::nullopt, StatSkipReason::insufficientSupport
        },
        StatValue{
            std::nullopt, StatSkipReason::insufficientSupport
        }
    };
    const auto skips = stat_skips (field, bothMissing);
    REQUIRE (skips.size() == 1);
    REQUIRE (
        to_info_format (skips[0]) == "QRK:insufficient_support"
    );
  }

  SECTION ("only the missing subfield contributes a skip")
  {
    // MLAS, since it is the only statistic that can reach a partially-missing
    // result: the joint-inference statistics skip whole or not at all.
    const VariantStatField mlasField{"MLAS", "", 2};
    const std::vector<StatValue> oneMissing{
        StatValue{0.7, std::nullopt},
        StatValue{std::nullopt, StatSkipReason::noBackground}
    };
    const auto skips = stat_skips (mlasField, oneMissing);
    REQUIRE (skips.size() == 1);
    REQUIRE (to_info_format (skips[0]) == "MLAS:no_background");
  }

  SECTION ("no missing subfields -> no skips")
  {
    const std::vector<StatValue> present{
        StatValue{1.0, std::nullopt},
        StatValue{2.0, std::nullopt}
    };
    REQUIRE (stat_skips (field, present).empty());
  }
}

TEST_CASE (
    "stat_all_missing expands a skip to the registry's arity"
)
{
  const auto twoValued = stat_all_missing (
      {"QRK", "", 2}, StatSkipReason::zeroVariance
  );
  REQUIRE (twoValued.size() == 2);
  for (const auto& v : twoValued) {
    REQUIRE_FALSE (v.value.has_value());
    REQUIRE (v.reason == StatSkipReason::zeroVariance);
  }

  const auto oneValued = stat_all_missing (
      {"RCMPLX", "", 1}, StatSkipReason::referenceHasN
  );
  REQUIRE (oneValued.size() == 1);
  REQUIRE (oneValued[0].reason == StatSkipReason::referenceHasN);
}
