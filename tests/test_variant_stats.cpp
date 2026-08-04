// Unit tests for the variant statistics compute layer (variant_stats.hpp).
// Exercised through the public variant_stats() registry, so the compute
// functions stay internal.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "expos/pileup_features.hpp"
#include "expos/variant_stats.hpp"

using Catch::Approx;

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
  const auto stats = variant_stats();
  REQUIRE (
      stats.size() == 4
  );  // QRK, TRK, MLAS, RCMPLX — the full package
  REQUIRE (by_id (stats, "QRK").field.nValues == 2);
  REQUIRE (by_id (stats, "TRK").field.nValues == 2);
  REQUIRE (by_id (stats, "MLAS").field.nValues == 2);
  REQUIRE (by_id (stats, "RCMPLX").field.nValues == 1);
}

TEST_CASE ("compute QRK (query-position clustering)")
{
  const auto qrk = by_id (variant_stats(), "QRK");
  std::mt19937 rng (1);

  SECTION ("missing with insufficient support")
  {
    PileupFeatures supporting;
    supporting.qPos = {10};
    PileupFeatures all;
    all.qPos = {0, 5, 10, 15, 20, 25};
    VariantStatInputs in{
        supporting, all, REF_PLACEHOLDER, rng, true
    };
    const auto r = qrk.compute (in);
    REQUIRE (r.size() == 2);
    REQUIRE_FALSE (r[0].value.has_value());
    REQUIRE (r[0].reason == REASON_INSUFFICIENT_SUPPORT);
    REQUIRE (r[1].reason == REASON_INSUFFICIENT_SUPPORT);
  }

  SECTION ("missing with insufficient background")
  {
    PileupFeatures supporting;
    supporting.qPos = {10, 11, 12};
    PileupFeatures all;
    all.qPos = {10, 11, 12, 40, 90};  // 5 < 2*3
    VariantStatInputs in{
        supporting, all, REF_PLACEHOLDER, rng, true
    };
    const auto r = qrk.compute (in);
    REQUIRE (r[0].reason == REASON_INSUFFICIENT_BACKGROUND);
  }

  SECTION ("suppressed when read lengths are heterogeneous")
  {
    PileupFeatures supporting;
    supporting.qPos = {50, 51, 52, 53};
    PileupFeatures all;
    for (int32_t i = 0; i < 200; ++i) {
      all.qPos.push_back (i);
    }
    VariantStatInputs in{
        supporting, all, REF_PLACEHOLDER, rng, false
    };
    const auto r = qrk.compute (in);
    REQUIRE_FALSE (r[0].value.has_value());
    REQUIRE_FALSE (r[1].value.has_value());
    REQUIRE (r[0].reason == REASON_HETEROGENEOUS_READ_LENGTH);
    REQUIRE (r[1].reason == REASON_HETEROGENEOUS_READ_LENGTH);
  }

  SECTION ("tight support vs spread background is extreme")
  {
    PileupFeatures supporting;
    supporting.qPos = {50, 51, 52, 53};  // all within radius 5
    PileupFeatures all;
    for (int32_t i = 0; i < 200; ++i) {
      all.qPos.push_back (i);  // spread 0..199
    }
    VariantStatInputs in{
        supporting, all, REF_PLACEHOLDER, rng, true
    };
    const auto r = qrk.compute (in);
    REQUIRE (r[0].value.has_value());  // effect size
    REQUIRE (r[1].value.has_value());  // p-value
    REQUIRE (*r[0].value > 0.0);  // more clustered than the null
    REQUIRE (*r[1].value < 0.05);  // significant
    REQUIRE (*r[1].value > 0.0);
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
    VariantStatInputs inA{
        supporting, all, REF_PLACEHOLDER, rngA, true
    };
    VariantStatInputs inB{
        supporting, all, REF_PLACEHOLDER, rngB, true
    };
    const auto a = qrk.compute (inA);
    const auto b = qrk.compute (inB);
    REQUIRE (a[0].value.has_value());
    REQUIRE (*a[0].value == Approx (*b[0].value));
    REQUIRE (*a[1].value == Approx (*b[1].value));
  }
}

TEST_CASE ("compute TRK (template-endpoint clustering)")
{
  const auto trk = by_id (variant_stats(), "TRK");

  SECTION ("missing with insufficient support")
  {
    std::mt19937 rng (1);
    PileupFeatures supporting;
    supporting.endpoints = {{100, 300}};
    PileupFeatures all;
    all.endpoints = {
        {100, 300}, {110, 320}, {500, 700}, {900, 1100}
    };
    VariantStatInputs in{
        supporting, all, REF_PLACEHOLDER, rng, true
    };
    const auto r = trk.compute (in);
    REQUIRE_FALSE (r[0].value.has_value());
    REQUIRE (r[0].reason == REASON_INSUFFICIENT_SUPPORT);
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
      // spread over 0..99 (manhattan 2*|i-j|), so draws have a spread of
      // within-radius counts and the null has non-zero variance.
      all.endpoints.push_back ({i, i + 200});
    }
    VariantStatInputs in{
        supporting, all, REF_PLACEHOLDER, rng, true
    };
    const auto r = trk.compute (in);
    REQUIRE (r[0].value.has_value());
    REQUIRE (*r[0].value > 0.0);
    REQUIRE (*r[1].value < 0.05);
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
    VariantStatInputs in{
        supporting, all, REF_PLACEHOLDER, rng, false
    };
    const auto r = trk.compute (in);
    REQUIRE (
        r[0].value.has_value()
    );  // computed despite heterogeneity
  }
}

TEST_CASE ("compute MLAS (median normalised alignment score)")
{
  const auto mlas = by_id (variant_stats(), "MLAS");
  std::mt19937 rng (1);

  SECTION ("medians of supporting and all")
  {
    PileupFeatures supporting;
    supporting.normalisedAs = {0.5, 0.7, 0.9};  // median 0.7
    PileupFeatures all;
    all.normalisedAs = {0.2, 0.4, 0.6, 0.8};  // median 0.5
    VariantStatInputs in{
        supporting, all, REF_PLACEHOLDER, rng, true
    };
    const auto r = mlas.compute (in);
    REQUIRE (r.size() == 2);
    REQUIRE (*r[0].value == Approx (0.7));
    REQUIRE (*r[1].value == Approx (0.5));
  }

  SECTION ("missing when the supporting group is empty")
  {
    PileupFeatures supporting;  // no reads
    PileupFeatures all;
    all.normalisedAs = {0.3, 0.6};
    VariantStatInputs in{
        supporting, all, REF_PLACEHOLDER, rng, true
    };
    const auto r = mlas.compute (in);
    REQUIRE_FALSE (r[0].value.has_value());
    REQUIRE (r[0].reason == REASON_NO_SUPPORT);
    REQUIRE (r[1].value.has_value());
  }
}

TEST_CASE ("compute RCMPLX (reference complexity)")
{
  const auto rcmplx = by_id (variant_stats(), "RCMPLX");
  std::mt19937 rng (1);
  const PileupFeatures empty;

  SECTION ("missing when the slice is shorter than the window")
  {
    const std::string ref (99, 'A');
    VariantStatInputs in{
        empty, empty, std::string_view (ref), rng, true
    };
    const auto r = rcmplx.compute (in);
    REQUIRE_FALSE (r[0].value.has_value());
    REQUIRE (r[0].reason == REASON_REFERENCE_TOO_SHORT);
  }

  SECTION ("missing when the slice contains N")
  {
    std::string ref (150, 'A');
    ref[75] = 'N';
    VariantStatInputs in{
        empty, empty, std::string_view (ref), rng, true
    };
    const auto r = rcmplx.compute (in);
    REQUIRE_FALSE (r[0].value.has_value());
    REQUIRE (r[0].reason == REASON_REFERENCE_HAS_N);
  }

  SECTION ("homopolymer window has low, known complexity")
  {
    const std::string ref (
        100, 'A'
    );  // exactly one 100-base window
    VariantStatInputs in{
        empty, empty, std::string_view (ref), rng, true
    };
    const auto r = rcmplx.compute (in);
    REQUIRE (r[0].value.has_value());
    // entropy_lz76 of 100 identical chars = 2*log2(100)/100 ~= 0.1329
    REQUIRE (*r[0].value == Approx (0.13287712));
  }
}

TEST_CASE ("stat_skip_tokens deduplicates within a statistic")
{
  const VariantStatField field{"QRK", "", 2};

  SECTION ("both subfields share a reason -> one token")
  {
    const std::vector<StatValue> bothMissing{
        stat_missing (REASON_INSUFFICIENT_SUPPORT),
        stat_missing (REASON_INSUFFICIENT_SUPPORT)
    };
    const auto tokens = stat_skip_tokens (field, bothMissing);
    REQUIRE (tokens.size() == 1);
    REQUIRE (tokens[0] == "QRK:insufficient_support");
  }

  SECTION ("only the missing subfield contributes a token")
  {
    const std::vector<StatValue> oneMissing{
        stat_value (1.0), stat_missing (REASON_ZERO_VARIANCE)
    };
    const auto tokens = stat_skip_tokens (field, oneMissing);
    REQUIRE (tokens.size() == 1);
    REQUIRE (tokens[0] == "QRK:zero_variance");
  }

  SECTION ("no missing subfields -> no tokens")
  {
    const std::vector<StatValue> present{
        stat_value (1.0), stat_value (2.0)
    };
    REQUIRE (stat_skip_tokens (field, present).empty());
  }
}
