// Unit tests for background_guard.hpp: pre-merge admission decisions for
// additional --bg samples, plus the primary sample's own read-length gate.

#include <catch2/catch_test_macros.hpp>
#include <cstdint>

#include "expos/background_guard.hpp"
#include "expos/pileup_features.hpp"

namespace {

constexpr std::size_t MIN_READS = 10;
constexpr std::size_t MIN_TEMPLATES = 10;
constexpr double REL_IQR_TOL = 0.10;
constexpr double MEDIAN_REL_TOL = 0.10;

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

TEST_CASE ("primary_read_length_homogeneous (Guard A)")
{
  SECTION ("homogeneous reads are admitted")
  {
    const auto primary = make_uniform (150, 300, 50);
    REQUIRE (primary_read_length_homogeneous (
        primary, MIN_READS, REL_IQR_TOL
    ));
  }

  SECTION ("heterogeneous reads fail")
  {
    PileupFeatures primary;
    for (int32_t i = 0; i < 50; ++i) {
      primary.readLen.push_back (i * 10);  // huge spread
    }
    REQUIRE_FALSE (primary_read_length_homogeneous (
        primary, MIN_READS, REL_IQR_TOL
    ));
  }

  SECTION ("insufficient data defaults to homogeneous")
  {
    PileupFeatures primary;
    primary.readLen = {10, 200};  // fewer than MIN_READS
    REQUIRE (primary_read_length_homogeneous (
        primary, MIN_READS, REL_IQR_TOL
    ));
  }
}

TEST_CASE ("evaluate_background (Guards B + C)")
{
  const auto primary = make_uniform (150, 300, 50);

  SECTION ("a well-matched candidate is admitted")
  {
    const auto candidate = make_uniform (150, 300, 50);
    REQUIRE (
        evaluate_background (
            primary, candidate, MIN_READS, REL_IQR_TOL,
            MIN_TEMPLATES, MEDIAN_REL_TOL
        ) == BackgroundGuardReason::Admitted
    );
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
    REQUIRE (
        evaluate_background (
            primary, candidate, MIN_READS, REL_IQR_TOL,
            MIN_TEMPLATES, MEDIAN_REL_TOL
        ) == BackgroundGuardReason::HeterogeneousOwnReadLength
    );
  }

  SECTION (
      "candidate with a mismatched median read length is "
      "excluded"
  )
  {
    const auto candidate =
        make_uniform (100, 300, 50);  // 150 vs 100
    REQUIRE (
        evaluate_background (
            primary, candidate, MIN_READS, REL_IQR_TOL,
            MIN_TEMPLATES, MEDIAN_REL_TOL
        ) == BackgroundGuardReason::ReadLengthMismatch
    );
  }

  SECTION (
      "candidate with a mismatched median fragment length is "
      "excluded"
  )
  {
    const auto candidate =
        make_uniform (150, 600, 50);  // 300 vs 600
    REQUIRE (
        evaluate_background (
            primary, candidate, MIN_READS, REL_IQR_TOL,
            MIN_TEMPLATES, MEDIAN_REL_TOL
        ) == BackgroundGuardReason::FragmentLengthMismatch
    );
  }

  SECTION ("insufficient reads excludes to be safe")
  {
    PileupFeatures candidate;
    candidate.readLen = {150, 150};  // fewer than MIN_READS
    candidate.endpoints = {{0, 300}, {1, 301}};
    REQUIRE (
        evaluate_background (
            primary, candidate, MIN_READS, REL_IQR_TOL,
            MIN_TEMPLATES, MEDIAN_REL_TOL
        ) == BackgroundGuardReason::InsufficientReadData
    );
  }

  SECTION ("insufficient templates excludes to be safe")
  {
    PileupFeatures candidate =
        make_uniform (150, 300, MIN_READS);
    candidate.endpoints.resize (
        2
    );  // enough reads, too few templates
    REQUIRE (
        evaluate_background (
            primary, candidate, MIN_READS, REL_IQR_TOL,
            MIN_TEMPLATES, MEDIAN_REL_TOL
        ) == BackgroundGuardReason::InsufficientFragmentData
    );
  }
}

TEST_CASE ("to_string covers every BackgroundGuardReason")
{
  REQUIRE (
      to_string (BackgroundGuardReason::Admitted) == "admitted"
  );
  REQUIRE_FALSE (
      to_string (BackgroundGuardReason::InsufficientReadData)
          .empty()
  );
  REQUIRE_FALSE (
      to_string (BackgroundGuardReason::InsufficientFragmentData)
          .empty()
  );
  REQUIRE_FALSE (
      to_string (
          BackgroundGuardReason::HeterogeneousOwnReadLength
      )
          .empty()
  );
  REQUIRE_FALSE (
      to_string (BackgroundGuardReason::ReadLengthMismatch)
          .empty()
  );
  REQUIRE_FALSE (
      to_string (BackgroundGuardReason::FragmentLengthMismatch)
          .empty()
  );
}
