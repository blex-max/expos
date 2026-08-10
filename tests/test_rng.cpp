#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <vector>

#include "shared/rng.hpp"

TEST_CASE ("Mwc192 matches the reference stream")
{
  // Mwc192 transcribes prng.di.unimi.it/MWC192.c by hand, because the
  // reference keeps its state in file-scope globals and cannot be
  // instantiated. A wrong shift or a swapped state word would still produce
  // plausible-looking noise, so these are the reference's own first outputs,
  // captured by seeding its globals exactly as the constructor seeds its
  // members. If the transcription ever drifts, this is what says so.
  Mwc192 rng (20260810);

  const std::uint64_t expected[] = {
      0x350222913F7385DBULL, 0x4389165430FCACC9ULL,
      0x13AD9440DF93EEF1ULL, 0xA9167DF10C728ABBULL,
      0xCF5D8EF8D54389E6ULL, 0x115E5F0EF141E7B5ULL,
      0x1B2EC0016D427470ULL, 0x86A81C726F927898ULL,
  };

  for (const std::uint64_t want : expected) {
    REQUIRE (rng() == want);
  }
}

TEST_CASE ("bounded stays inside its range")
{
  Mwc192 rng (7);

  // Loops tally and assert once; a REQUIRE per draw would put a hundred
  // thousand assertions on the counter and say nothing more.
  std::size_t violations = 0;
  for (std::size_t k = 0; k < 1000; ++k) {
    violations += (bounded (rng, 1) != 0) ? 1 : 0;
  }
  REQUIRE (violations == 0);

  for (const std::uint64_t range :
       {2ULL, 3ULL, 17ULL, 79ULL, 1000ULL}) {
    CAPTURE (range);
    violations = 0;
    for (std::size_t k = 0; k < 10000; ++k) {
      violations += (bounded (rng, range) >= range) ? 1 : 0;
    }
    REQUIRE (violations == 0);
  }
}

TEST_CASE ("the half-open range maps onto an inclusive draw")
{
  // subsample_wo_replace wants a value in [i, nObs - 1] inclusive, which is
  // what std::uniform_int_distribution (i, nObs - 1) gave it. Lemire is
  // half-open, so the range is nObs - i. Off by one in either direction is a
  // biased Fisher-Yates that still passes every determinism test, which is
  // why this checks the distribution and not just the bounds.
  constexpr std::size_t nObs = 12;
  Mwc192 rng (3);

  for (std::size_t i = 0; i < nObs; ++i) {
    CAPTURE (i);
    std::vector<std::size_t> seen (nObs, 0);
    std::size_t outOfRange = 0;
    for (std::size_t k = 0; k < 20000; ++k) {
      const std::size_t j = i + bounded (rng, nObs - i);
      if (j < i || j >= nObs) {
        ++outOfRange;
        continue;
      }
      ++seen[j];
    }
    REQUIRE (outOfRange == 0);

    // every value in [i, nObs - 1] appears, and nothing below i does
    const auto split =
        seen.begin() + static_cast<std::ptrdiff_t> (i);
    REQUIRE (
        std::all_of (seen.begin(), split, [] (std::size_t h) {
          return h == 0;
        })
    );
    REQUIRE (std::all_of (split, seen.end(), [] (std::size_t h) {
      return h > 0;
    }));
  }
}
