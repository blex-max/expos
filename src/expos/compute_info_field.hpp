#pragma once

#include <cassert>
#include <cstdint>
#include <expected>
#include <numeric>
#include <optional>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include "expos/pileup_features.hpp"
#include "expos/skip.hpp"
#include "shared/rng.hpp"

// --- monte carlo --- //

// one-sided pval; the effect size is a z-score against the
// null's population standard deviation.
template <typename DrawFn, typename StatFn>
double run_monte_carlo (double obsVal, DrawFn draw, StatFn stat)
{
  constexpr uint16_t NSIM_PVAL = 200; // resolution of 0.005
  uint16_t countGe = 0;
  for (uint16_t k = 1; k <= NSIM_PVAL; ++k) {
    const double s = static_cast<double> (stat (draw()));
    if (s >= obsVal) {
      ++countGe;
    }
  }

  const auto pVal = static_cast<double> (countGe + 1) /
                    static_cast<double> (NSIM_PVAL + 1);
  return pVal;
}

// Reusable buffers for repeated draws from one population.
template <class T>
struct SubsampleScratch {
  std::vector<uint64_t> idxBuf;
  std::vector<T> subsampleBuf;
};

// Per-run mutable state for the Monte-Carlo statistics: the shared RNG stream
// and one scratch per statistic. Built once for the whole record loop, so the
// buffers reach a high-water mark and stop reallocating. A new Monte-Carlo
// statistic adds a scratch here alongside its compute_ function and registry
// row.
struct McState {
  Mwc192 rng;
  SubsampleScratch<int32_t> qPos;
  SubsampleScratch<TemplateEndpoints> endpoints;
};

// Uniform sample of size n drawn without replacement from obs, via a partial
// Fisher-Yates shuffle of scratch's persistent index permutation. The returned
// span views scratch.subsampleBuf and is valid until the next draw from the
// same scratch.
// Precondition: n <= obs.size.
template <class T>
std::span<T> subsample_wo_replace (
    const std::vector<T>& obs, uint64_t n, Mwc192& rng,
    SubsampleScratch<T>& scratch
)
{
  const uint64_t nObs = obs.size();
  assert (n <= nObs);
  if (scratch.idxBuf.size() != nObs) {
    scratch.idxBuf.resize (nObs);
    std::iota (
        scratch.idxBuf.begin(), scratch.idxBuf.end(), uint64_t{0}
    );
  }
  if (scratch.subsampleBuf.size() < n) {
    scratch.subsampleBuf.resize (n);
  }
  for (uint64_t i = 0; i < n; ++i) {
    // dist (i, nObs - 1) was inclusive at both ends,
    // so the half-open range is nObs - i.
    const uint64_t j = i + bounded (rng, nObs - i);
    std::swap (scratch.idxBuf[i], scratch.idxBuf[j]);
    scratch.subsampleBuf[i] = obs[scratch.idxBuf[i]];
  }
  // Only the first n entries were written this call
  // (any further entries stale)
  return std::span<T>{scratch.subsampleBuf.data(), n};
}

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
