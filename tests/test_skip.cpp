// Unit tests for the EXPOS_SKIP vocabulary. These strings are user-facing —
// they are the reason half of every EXPOS_SKIP token and are promised by that
// field's header description — so they are pinned here against silent
// rewording by a refactor.

#include <catch2/catch_test_macros.hpp>
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include "expos/skip.hpp"

namespace {

// Every enumerator, so a new one added without a string fails here as well
// as at the -Wswitch in to_string.
const std::vector<RecordSkipReason> ALL_RECORD_REASONS{
    RecordSkipReason::notBiallelic,
    RecordSkipReason::complex,
};

const std::vector<StatSkipReason> ALL_STAT_REASONS{
    StatSkipReason::insufficientSupport,
    StatSkipReason::insufficientBackground,
    StatSkipReason::heterogeneousReadLength,
    StatSkipReason::readLengthUnverified,
    StatSkipReason::zeroVariance,
    StatSkipReason::noSupport,
    StatSkipReason::noBackground,
    StatSkipReason::referenceTooShort,
    StatSkipReason::referenceHasN,
};

}  // namespace

TEST_CASE (
    "every skip reason renders to a distinct non-empty string"
)
{
  std::set<std::string_view> seen;

  for (const auto reason : ALL_RECORD_REASONS) {
    const auto s = to_string (reason);
    REQUIRE_FALSE (s.empty());
    REQUIRE (seen.insert (s).second);  // no duplicates
  }
  for (const auto reason : ALL_STAT_REASONS) {
    const auto s = to_string (reason);
    REQUIRE_FALSE (s.empty());
    REQUIRE (seen.insert (s).second);
  }

  REQUIRE (
      seen.size() ==
      ALL_RECORD_REASONS.size() + ALL_STAT_REASONS.size()
  );
}

TEST_CASE ("skip reason strings are stable")
{
  // Reword only with a deliberate decision: downstream consumers parse these.
  REQUIRE (
      to_string (RecordSkipReason::notBiallelic) ==
      "not_biallelic"
  );
  REQUIRE (to_string (RecordSkipReason::complex) == "complex");

  REQUIRE (
      to_string (StatSkipReason::insufficientSupport) ==
      "insufficient_support"
  );
  REQUIRE (
      to_string (StatSkipReason::insufficientBackground) ==
      "insufficient_background"
  );
  REQUIRE (
      to_string (StatSkipReason::heterogeneousReadLength) ==
      "heterogeneous_read_length"
  );
  REQUIRE (
      to_string (StatSkipReason::readLengthUnverified) ==
      "read_length_unverified"
  );
  REQUIRE (
      to_string (StatSkipReason::zeroVariance) == "zero_variance"
  );
  REQUIRE (
      to_string (StatSkipReason::noSupport) == "no_support"
  );
  REQUIRE (
      to_string (StatSkipReason::noBackground) == "no_background"
  );
  REQUIRE (
      to_string (StatSkipReason::referenceTooShort) ==
      "reference_too_short"
  );
  REQUIRE (
      to_string (StatSkipReason::referenceHasN) ==
      "reference_has_n"
  );
}

TEST_CASE ("Skip renders as <scope>:<reason>")
{
  SECTION ("record scope")
  {
    const auto skip =
        make_record_skip (RecordSkipReason::complex);
    REQUIRE (skip.scope == SCOPE_RECORD);
    REQUIRE (to_info_format (skip) == "record:complex");
  }

  SECTION ("statistic scope uses the statistic's INFO ID")
  {
    const auto skip =
        make_stat_skip ("QRK", StatSkipReason::zeroVariance);
    REQUIRE (skip.scope == "QRK");
    REQUIRE (to_info_format (skip) == "QRK:zero_variance");
  }

  SECTION ("equality covers both halves")
  {
    const auto a =
        make_stat_skip ("QRK", StatSkipReason::zeroVariance);
    const auto b =
        make_stat_skip ("QRK", StatSkipReason::zeroVariance);
    const auto differentScope =
        make_stat_skip ("TJAC", StatSkipReason::zeroVariance);
    const auto differentReason =
        make_stat_skip ("QRK", StatSkipReason::noSupport);
    REQUIRE (a == b);
    REQUIRE_FALSE (a == differentScope);
    REQUIRE_FALSE (a == differentReason);
  }
}
