#pragma once

#include <cassert>
#include <cstdint>
#include <limits>

// Monte-Carlo random source: MWC192 for the bit stream, Lemire for the
// bounded draw. Much faster than std's mersenne twister, which was
// a large cost when profiled.
//
// Sources, both public domain (CC0):
//   MWC192   Vigna, 2023. https://prng.di.unimi.it/MWC192.c
//   splitmix64  Steele, Lea & Flood, OOPSLA 2014.
//            https://prng.di.unimi.it/splitmix64.c
//   Lemire   "Fast Random Integer Generation in an Interval", ACM TOMACS
//            29(1), 2019. arXiv:1805.10941, section 4.
//
// Vigna notes of MWC192 that its modulus "has a particular form, which creates
// some theoretical issues, but at this size a generator of this kind passes
// all known statistical tests".
//
// A Marsaglia multiply-with-carry generator, whatever that means.
// State is three uint64s against mt19937's ~2.5 KB, which may help
// keep the L1 cache hot for expos.
class Mwc192 {
 public:
  using result_type = std::uint64_t;

  static constexpr result_type min() { return 0; }
  static constexpr result_type max()
  {
    return std::numeric_limits<result_type>::max();
  }

  // The reference requires 0 < c < MWC_A2 - 1 and suggests c = 1 with x and y
  // taking a 128-bit seed. splitmix64 expands the seed to fill them.
  explicit Mwc192 (std::uint64_t seed)
  {
    std::uint64_t z = seed;
    const auto splitmix64 = [&z] {
      // Apparently:
      // A Weyl sequence stepped by the 64-bit golden ratio, then a
      // MurmurHash3-style finaliser; the two odd multipliers are the
      // constants that finaliser was tuned to.
      z += 0x9E3779B97F4A7C15ULL;
      std::uint64_t r = z;
      r = (r ^ (r >> 30)) * 0xBF58476D1CE4E5B9ULL;
      r = (r ^ (r >> 27)) * 0x94D049BB133111EBULL;
      return r ^ (r >> 31);
    };
    x = splitmix64();
    y = splitmix64();
    c = 1;
    // Holds trivially for c = 1, and is the invariant to re-check if anyone
    // ever derives c from the seed as well.
    assert (c > 0 && c < MWC_A2 - 1);
  }

  // equivalent to next() in MWC192.c
  result_type operator()()
  {
    __extension__ using u128 =
        unsigned __int128;  // GCC -Wpedantic
    const std::uint64_t result = y;
    const u128 t = static_cast<u128> (MWC_A2) * x + c;
    x = y;
    y = static_cast<std::uint64_t> (t);
    c = static_cast<std::uint64_t> (t >> 64);
    return result;
  }

 private:
  static constexpr std::uint64_t MWC_A2 = 0xFFA04E67B3C95D86ULL;

  std::uint64_t x;
  std::uint64_t y;
  std::uint64_t c;
};

// Lemire's Nearly Divisionless Bounded Draw
//
// Uniform in [0, range). Precondition: range > 0.
inline std::uint64_t bounded (Mwc192& rng, std::uint64_t range)
{
  assert (range > 0);
  __extension__ using u128 =
      unsigned __int128;  // GCC -Wpedantic
  u128 m = static_cast<u128> (rng()) * static_cast<u128> (range);
  auto low = static_cast<std::uint64_t> (m);
  if (low < range) {
    // 2^64 mod range, without a 128-bit intermediate.
    const std::uint64_t rejectBelow =
        (std::numeric_limits<std::uint64_t>::max() - range + 1) %
        range;
    while (low < rejectBelow) {
      m = static_cast<u128> (rng()) * static_cast<u128> (range);
      low = static_cast<std::uint64_t> (m);
    }
  }
  return static_cast<std::uint64_t> (m >> 64);
}
