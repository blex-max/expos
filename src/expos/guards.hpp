#pragma once

#include <cstdint>
#include <expected>
#include <optional>
#include <string_view>
#include <vector>

#include "expos/pileup_features.hpp"
#include "expos/skip.hpp"

bool sufficient_reads (uint64_t nReads);

bool read_lens_within_tol (const std::vector<int32_t>& readLen);

// Whether a Monte-Carlo statistic has enough to work with
// or the reason it does not.
std::optional<StatSkipReason> size_guard (
    uint64_t nObs, uint64_t nBackground
);

// Aggregate stats of primary sample data at one locus,
// for guard calculations
struct PrimaryGuardStats {
  uint64_t nReads = 0;
  uint64_t nTemplates = 0;
  std::optional<double> medianReadLen;  // engaged iff nReads > 0
  std::optional<double>
      medianFragLen;  // engaged iff nTemplates > 0
};
PrimaryGuardStats summarise_primary (
    const PileupFeatures& primaryAll
);


// Test whether a candidate background sample's local read population at
// a given record is compatible enough with the primary sample's to be
// merged into the Monte-Carlo background pool.
enum class BgGuardFail : uint8_t {
  InsufficientReadData,  // can't evaluate read-length guards
  InsufficientFragmentData,  // can't evaluate the fragment-length guard
  HeterogeneousReadLength,  // the candidate's own reads are too spread
  ReadLengthMismatch,  // median read length differs from primary
  FragmentLengthMismatch,  // median fragment length differs from primary
};
std::string_view to_string (BgGuardFail reason);
using VoidOrFail = std::expected<void, BgGuardFail>;
VoidOrFail verify_bg_sample (
    const PrimaryGuardStats& primary,
    const PileupFeatures& candidate
);
