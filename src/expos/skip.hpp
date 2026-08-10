#pragma once

#include <cstdint>
#include <string>
#include <string_view>

// --- EXPOS_SKIP vocabulary --- //

// Why expos produced no value for something. Two disjoint sets: a whole
// record can be skipped, or a single statistic can be.
//
// A record that cannot be analysed at all.
enum class RecordSkipReason : uint8_t { notBiallelic, complex };

// A statistic that could not be computed for an otherwise analysable record.
enum class StatSkipReason : uint8_t {
  insufficientSupport,
  insufficientBackground,
  heterogeneousReadLength,
  readLengthUnverified,
  zeroVariance,
  noSupport,
  noBackground,
  referenceTooShort,
  referenceHasN
};

// These strings are user-facing: they are the reason half of every
// EXPOS_SKIP token, and are promised by that field's header description.
std::string_view to_string (RecordSkipReason reason);
std::string_view to_string (StatSkipReason reason);

// The scope used for record-level skips. Statistic-level skips use the
// statistic's INFO ID, which already lives in the registry.
inline constexpr std::string_view SCOPE_RECORD = "record";

// One EXPOS_SKIP entry. Built with make_ factories
struct Skip {
  std::string_view scope;
  std::string_view reason;
};

inline bool operator== (const Skip& a, const Skip& b)
{
  return a.scope == b.scope && a.reason == b.reason;
}

Skip make_record_skip (RecordSkipReason reason);
Skip make_stat_skip (
    std::string_view statId, StatSkipReason reason
);

// format for info field
std::string to_info_format (const Skip& skip);
