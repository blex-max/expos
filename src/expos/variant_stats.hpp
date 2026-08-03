#pragma once

#include <htslib/vcf.h>

#include <optional>
#include <random>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "expos/pileup_features.hpp"
#include "shared/err.hpp"

// A VariantStat is one domain-specific statistic annotated onto a variant
// as a VCF INFO field. The per-variant loop iterates variant_stats() and
// encodes each result. Adding a statistic is a compute_* free function plus
// a registry row in variant_stats.cpp.

struct VariantStatField {
  std::string_view id;  // INFO ID, e.g. "QRK"
  std::string_view headerLine;  // full ##INFO=<...> definition
  int nValues;  // number of Float subfields
};

// One subfield result: either a value, or the reason it is missing (a
// REASON_* vocabulary constant). Exactly one of the two is set. Reasons are
// surfaced to the reader in the EXPOS_SKIP INFO field.
struct StatValue {
  std::optional<double> value;
  std::optional<std::string_view> reason;
};

// Builders keeping the compute_* functions readable.
inline StatValue stat_value (double v)
{
  return {v, std::nullopt};
}
inline StatValue stat_missing (std::string_view reason)
{
  return {std::nullopt, reason};
}
inline StatValue stat_or (
    std::optional<double> v, std::string_view reasonIfMissing
)
{
  return v ? stat_value (*v) : stat_missing (reasonIfMissing);
}

struct VariantStatInputs {
  const PileupFeatures& supporting;  // observed sample reads
  const PileupFeatures&
      all;  // all reads, including additional samples
  std::string_view refSlice;  // reference span
  std::mt19937& rng;
  bool
      readLenHomogeneous;  // false suppresses read-length-sensitive stats
};

// Map inputs to one StatValue per subfield.
using VariantStatFn =
    std::vector<StatValue> (*) (const VariantStatInputs&);

struct VariantStat {
  VariantStatField field;
  VariantStatFn compute;
};

std::span<const VariantStat> variant_stats();

// --- EXPOS_SKIP reason vocabulary --- //
// Record-level (whole record not analysed):
inline constexpr std::string_view REASON_MULTIALLELIC =
    "multiallelic";
inline constexpr std::string_view REASON_COMPLEX = "complex";
// Statistic-level:
inline constexpr std::string_view REASON_INSUFFICIENT_SUPPORT =
    "insufficient_support";
inline constexpr std::string_view
    REASON_INSUFFICIENT_BACKGROUND = "insufficient_background";
inline constexpr std::string_view
    REASON_HETEROGENEOUS_READ_LENGTH =
        "heterogeneous_read_length";
inline constexpr std::string_view REASON_ZERO_VARIANCE =
    "zero_variance";
inline constexpr std::string_view REASON_NO_SUPPORT =
    "no_support";
inline constexpr std::string_view REASON_NO_BACKGROUND =
    "no_background";
inline constexpr std::string_view REASON_REFERENCE_TOO_SHORT =
    "reference_too_short";
inline constexpr std::string_view REASON_REFERENCE_HAS_N =
    "reference_has_n";

// --- VCF output --- //

// Append the ##INFO header line for one statistic.
[[nodiscard]] VoidOrErr register_variant_stat_header (
    bcf_hdr_t* hdr, const VariantStat& stat
);

// Append the ##INFO header line for the EXPOS_SKIP field.
[[nodiscard]] VoidOrErr register_expos_skip_header (
    bcf_hdr_t* hdr
);

// Write a statistic's per-subfield values to `rec` as its Float INFO field.
// Missing subfields become VCF missing values. `values.size()` must equal
// `field.nValues`.
[[nodiscard]] VoidOrErr encode_variant_stat (
    bcf_hdr_t* hdr, bcf1_t* rec, const VariantStatField& field,
    const std::vector<StatValue>& values
);

// Deduplicated "<ID>:<reason>" tokens for the missing subfields of `values`.
std::vector<std::string> stat_skip_tokens (
    const VariantStatField& field,
    const std::vector<StatValue>& values
);

// Write the EXPOS_SKIP INFO field from a set of "<scope>:<reason>" tokens.
// No-op when tokens is empty.
[[nodiscard]] VoidOrErr set_expos_skip (
    bcf_hdr_t* hdr, bcf1_t* rec,
    const std::vector<std::string>& tokens
);
