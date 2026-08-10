#pragma once

#include <expected>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

#include "expos/pileup_features.hpp"
#include "expos/skip.hpp"
#include "expos/variant_stats.hpp"

// --- generic infrastructure for calculating all variant stats --- //

// Per-record features
struct VariantStatInputs {
  const PileupFeatures& supporting;  // observed sample reads
  const PileupFeatures&
      all;  // all reads, including additional samples
  std::string_view refSlice;  // reference span
};

struct StatContext {
  McState& mc;
  // Engaged when read-length-sensitive stats must be suppressed, and
  // carries the reason so the skip token says which of the two it was:
  // the reads disagreed, or there were too few to tell.
  std::optional<StatSkipReason> readLenSuppression;
};

// Check the primary sample's own local read population internally
// consistent enough for QRK.
void set_qrk_guard (
    StatContext& sCtx, const PileupFeatures& primaryAll
);

// A VariantStat is one domain-specific statistic annotated onto a variant
// as a VCF INFO field. The per-variant loop iterates expos_field_registry()
// and encodes each result. Adding a statistic is a compute_* free function
// plus a registry row in compute_info_field.cpp.

struct VariantStatField {
  std::string_view id;  // INFO ID, e.g. "QRK"
  std::string_view headerLine;  // full ##INFO=<...> definition
  int nValues;  // number of Float subfields
};

// One subfield result: either a value, XOR the reason it is missing.
// Reasons stay as enums this far down; they are rendered into EXPOS_SKIP
// tokens at encode time.
struct StatValue {
  std::optional<double> value;
  std::optional<StatSkipReason> reason;
};

// A compute function returns either one StatValue per subfield, or a single
// reason the whole statistic is missing.
using ValuesOrSkip =
    std::expected<std::vector<StatValue>, StatSkipReason>;
using ComputeFn = ValuesOrSkip (*) (
    const VariantStatInputs&, const StatContext&
);

struct VariantStat {
  VariantStatField field;
  ComputeFn compute;
};

// get stat registry
std::span<const VariantStat> expos_field_registry();
