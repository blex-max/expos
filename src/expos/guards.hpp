#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <string_view>
#include <vector>

#include "expos/pileup_features.hpp"
#include "expos/skip.hpp"

bool sufficient_reads (std::size_t nReads);

bool read_lens_within_tol (const std::vector<int32_t>& readLen);

// Whether a Monte-Carlo statistic has enough to work with, or the reason it
// does not. Deliberately takes plain counts rather than a PileupFeatures:
// each statistic sizes on its own feature vector, and those count different
// things -- qPos drops deleted and ref-skipped reads, endpoints counts
// templates and is empty altogether for unpaired data. A caller passing the
// wrong pair would get a verdict on a population its statistic never uses.
std::optional<StatSkipReason> size_guard (
    std::size_t nObs, std::size_t nBackground
);

// Decides whether a candidate background sample's local read population at
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

// Everything the guards need to know about the primary sample at one
// locus. Taken once per record, before any background is merged: the
// primary is the fixed point of reference for every candidate, so the
// medians are computed once here rather than re-derived (each a copy and
// a sort) for every candidate in turn.
struct PrimaryGuardStats {
  std::size_t nReads = 0;
  std::size_t nTemplates = 0;
  std::optional<double> medianReadLen;  // engaged iff nReads > 0
  std::optional<double>
      medianFragLen;  // engaged iff nTemplates > 0
};

PrimaryGuardStats summarise_primary (
    const PileupFeatures& primaryAll
);

// Check whether candidate be merged into the background pool.
// Insufficient data to evaluate any guard excludes the
// candidate.
using VoidOrFail = std::expected<void, BgGuardFail>;
VoidOrFail verify_bg_sample (
    const PrimaryGuardStats& primary,
    const PileupFeatures& candidate
);
