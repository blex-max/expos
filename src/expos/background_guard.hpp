#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

#include "expos/pileup_features.hpp"

// Decides whether a candidate background sample's local read population at
// a given record is compatible enough with the primary sample's to be
// merged into the Monte-Carlo background pool.

enum class BackgroundGuardReason : std::uint8_t {
  Admitted,
  InsufficientReadData,  // can't evaluate read-length guards
  InsufficientFragmentData,  // can't evaluate the fragment-length guard
  HeterogeneousOwnReadLength,  // the candidate's own reads are too spread
  ReadLengthMismatch,  // median read length differs from primary
  FragmentLengthMismatch,  // median fragment length differs from primary
};

std::string_view to_string (BackgroundGuardReason reason);

// Check the primary sample's own local read population internally
// consistent enough for QRK.
bool primary_read_length_homogeneous (
    const PileupFeatures& primaryAll, std::size_t minReads,
    double relIqrTol
);

// Check whether candidate be merged into the background pool.
// Insufficient data to evaluate any guard excludes the
// candidate.
BackgroundGuardReason evaluate_background (
    const PileupFeatures& primaryAll,
    const PileupFeatures& candidate, std::size_t minReads,
    double relIqrTol, std::size_t minTemplates,
    double medianRelTol
);
