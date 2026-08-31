// Unit tests for the pure statistical primitives in shared/stats.hpp.
// Standard-library-only, so they need no expos types. The clustering and
// Monte-Carlo primitives live in tests/test_variant_stats.cpp alongside the
// statistics that consume them.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <vector>

#include "shared/stats.hpp"

using Catch::Approx;

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

TEST_CASE ("lz76")
{
  // exact LZ76 phrase counts
  REQUIRE (lz76 ("AAAA") == 2);
  REQUIRE (lz76 ("ABAB") == 3);
  REQUIRE (lz76 ("ABCD") == 4);

  REQUIRE (lz76 ("") == 0);
  REQUIRE (lz76 ("A") == 1);

  // repetitive is less complex than diverse, at equal length
  REQUIRE (
      lz76 ("AAAAAAAAAAAAAAAA") < lz76 ("ACGTACGATCGGATCA")
  );
}
