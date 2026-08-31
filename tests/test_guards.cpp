// Unit tests for guards.hpp: pre-merge admission decisions for
// additional -b/--background-sample samples, plus the primary sample's own read-length gate.

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <optional>

#include "expos/compute_info_field.hpp"
#include "expos/guards.hpp"
#include "expos/pileup_features.hpp"

namespace {

PileupFeatures make_uniform (
    int32_t readLen, int64_t fragLen, std::size_t n
)
{
  PileupFeatures pf;
  for (std::size_t i = 0; i < n; ++i) {
    pf.readLen.push_back (readLen);
    pf.endpoints.emplace_back (
        static_cast<int64_t> (i),
        static_cast<int64_t> (i) + fragLen
    );
  }
  return pf;
}

}  // namespace

TEST_CASE ("set_qrk_guard (Guard A)")
{
  Mwc192 rng (1);
  McState mc{rng, {}, {}};

  SECTION ("homogeneous reads leave QRK enabled")
  {
    const auto primary = make_uniform (150, 300, 50);
    StatContext ctx{mc, std::nullopt};
    set_qrk_guard (ctx, primary);
    REQUIRE_FALSE (ctx.readLenSuppression.has_value());
  }

  // Each of the sections below seeds ctx with an unrelated reason.

  SECTION ("heterogeneous reads suppress QRK as such")
  {
    PileupFeatures primary;
    for (int32_t i = 0; i < 50; ++i) {
      primary.readLen.push_back (i * 10);  // huge spread
    }
    StatContext ctx{mc, StatSkipReason::zeroVariance};
    set_qrk_guard (ctx, primary);
    REQUIRE (
        ctx.readLenSuppression ==
        StatSkipReason::heterogeneousReadLength
    );
  }

  // Too few reads to judge the spread at all, so fail closed and suppress
  // QRK.
  SECTION ("insufficient data suppresses QRK as unverified")
  {
    PileupFeatures primary;
    primary.readLen = {150, 150};  // fewer than MIN_READS
    StatContext ctx{mc, StatSkipReason::zeroVariance};
    set_qrk_guard (ctx, primary);
    REQUIRE (
        ctx.readLenSuppression ==
        StatSkipReason::readLengthUnverified
    );
  }

  SECTION ("no reads at all suppresses QRK as unverified")
  {
    const PileupFeatures primary;
    StatContext ctx{mc, StatSkipReason::zeroVariance};
    set_qrk_guard (ctx, primary);
    REQUIRE (
        ctx.readLenSuppression ==
        StatSkipReason::readLengthUnverified
    );
  }
}

TEST_CASE ("verify_bg_sample (Guards B + C)")
{
  const auto primary =
      summarise_primary (make_uniform (150, 300, 50));

  SECTION ("a well-matched candidate is admitted")
  {
    const auto candidate = make_uniform (150, 300, 50);
    REQUIRE (verify_bg_sample (primary, candidate).has_value());
  }

  SECTION (
      "candidate with heterogeneous own read lengths is excluded"
  )
  {
    PileupFeatures candidate;
    for (int32_t i = 0; i < 50; ++i) {
      candidate.readLen.push_back (
          i * 10
      );  // huge internal spread
      candidate.endpoints.emplace_back (i, i + 300);
    }
    const auto ret = verify_bg_sample (primary, candidate);
    REQUIRE_FALSE (ret.has_value());
    REQUIRE (
        ret.error() == BgGuardFail::HeterogeneousReadLength
    );
  }

  SECTION (
      "candidate with a mismatched median read length is "
      "excluded"
  )
  {
    const auto candidate =
        make_uniform (100, 300, 50);  // 150 vs 100
    const auto ret = verify_bg_sample (primary, candidate);
    REQUIRE_FALSE (ret.has_value());
    REQUIRE (ret.error() == BgGuardFail::ReadLengthMismatch);
  }

  SECTION (
      "candidate with a mismatched median fragment length is "
      "excluded"
  )
  {
    const auto candidate =
        make_uniform (150, 600, 50);  // 300 vs 600
    const auto ret = verify_bg_sample (primary, candidate);
    REQUIRE_FALSE (ret.has_value());
    REQUIRE (ret.error() == BgGuardFail::FragmentLengthMismatch);
  }

  SECTION ("insufficient reads excludes to be safe")
  {
    PileupFeatures candidate;
    candidate.readLen = {150, 150};  // fewer than MIN_READS
    candidate.endpoints = {{0, 300}, {1, 301}};
    const auto ret = verify_bg_sample (primary, candidate);
    REQUIRE_FALSE (ret.has_value());
    REQUIRE (ret.error() == BgGuardFail::InsufficientReadData);
  }

  SECTION ("insufficient templates excludes to be safe")
  {
    PileupFeatures candidate =
        make_uniform (150, 300, MIN_READS);
    candidate.endpoints.resize (
        MIN_TEMPLATES - 1
    );  // enough reads, too few templates
    const auto ret = verify_bg_sample (primary, candidate);
    REQUIRE_FALSE (ret.has_value());
    REQUIRE (
        ret.error() == BgGuardFail::InsufficientFragmentData
    );
  }

  SECTION ("a thin primary excludes an otherwise-good candidate")
  {
    PileupFeatures thinPrimary;
    thinPrimary.readLen = {150, 150};
    thinPrimary.endpoints = {{0, 300}, {1, 301}};
    const auto candidate = make_uniform (150, 300, 50);
    const auto ret = verify_bg_sample (
        summarise_primary (thinPrimary), candidate
    );
    REQUIRE_FALSE (ret.has_value());
    REQUIRE (ret.error() == BgGuardFail::InsufficientReadData);
  }
}

TEST_CASE ("to_string covers every BgGuardFail")
{
  REQUIRE_FALSE (
      to_string (BgGuardFail::InsufficientReadData).empty()
  );
  REQUIRE_FALSE (
      to_string (BgGuardFail::InsufficientFragmentData).empty()
  );
  REQUIRE_FALSE (
      to_string (BgGuardFail::HeterogeneousReadLength).empty()
  );
  REQUIRE_FALSE (
      to_string (BgGuardFail::ReadLengthMismatch).empty()
  );
  REQUIRE_FALSE (
      to_string (BgGuardFail::FragmentLengthMismatch).empty()
  );
}
