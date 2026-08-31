#include "guards.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <ranges>
#include <vector>

#include "shared/stats.hpp"

bool sufficient_reads (uint64_t nReads)
{
  return nReads >= MIN_READS;
}

// A statistic needs >= 2 supporting observations, and a background that is
// both at least twice the support (very, very liberal) and no smaller than
// MIN_BACKGROUND in absolute terms.
//
// min_background ensures that there are enough pairs in the
// background set to ensure a sufficiently granular p-val
//
// Below min_obs at 2 there are no pairs at all
std::optional<StatSkipReason> size_guard (
    uint64_t nObs, uint64_t nBackground
)
{
  if (nObs < MIN_OBS) {
    return StatSkipReason::insufficientSupport;
  }
  if (nBackground < MIN_BACKGROUND || nBackground < (2 * nObs)) {
    return StatSkipReason::insufficientBackground;
  }
  return std::nullopt;
}

// Relative IQR of read lengths, (Q3-Q1)/median. Precondition: non-empty.
bool read_lens_within_tol (const std::vector<int32_t>& readLen)
{
  std::vector<int32_t> lens = readLen;
  std::sort (lens.begin(), lens.end());
  const uint64_t n = lens.size();
  const double q1 = lens[n / 4];
  const double median = lens[n / 2];
  const double q3 = lens[3 * n / 4];
  return ((q3 - q1) / median) <=
         READ_LEN_REL_IQR_TOL;  // median >= 1: read length is always positive
}

static std::optional<double> median_read_length (
    const std::vector<int32_t>& readLen
)
{
  return percentile (
      std::views::transform (
          readLen,
          [] (int32_t len) {
            return static_cast<uint32_t> (len);
          }
      ),
      PERCENTILE_MEDIAN
  );
}

static std::optional<double> median_fragment_length (
    const std::vector<TemplateEndpoints>& endpoints
)
{
  return percentile (
      std::views::transform (
          endpoints,
          [] (const auto& span) {
            return static_cast<uint32_t> (
                span.second - span.first
            );
          }
      ),
      PERCENTILE_MEDIAN
  );
}

std::string_view to_string (BgGuardFail reason)
{
  switch (reason) {
    case BgGuardFail::InsufficientReadData:
      return "insufficient reads to check read-length "
             "consistency";
    case BgGuardFail::InsufficientFragmentData:
      return "insufficient templates to check fragment-length "
             "consistency";
    case BgGuardFail::HeterogeneousReadLength:
      return "read lengths are too heterogeneous";
    case BgGuardFail::ReadLengthMismatch:
      return "median read length differs from the primary "
             "sample";
    case BgGuardFail::FragmentLengthMismatch:
      return "median fragment length differs from the "
             "primary "
             "sample";
  }
}

PrimaryGuardStats summarise_primary (
    const PileupFeatures& primaryAll
)
{
  return PrimaryGuardStats{
      primaryAll.readLen.size(), primaryAll.endpoints.size(),
      median_read_length (primaryAll.readLen),
      median_fragment_length (primaryAll.endpoints)
  };
}

VoidOrFail verify_bg_sample (
    const PrimaryGuardStats& primary,
    const PileupFeatures& candidate
)
{
  // failed-closed, i.e. insufficient
  // data to be certain of an outcome
  // = exclusion.
  if (!sufficient_reads (candidate.readLen.size()) ||
      !sufficient_reads (primary.nReads)) {
    return std::unexpected (BgGuardFail::InsufficientReadData);
  }
  if (candidate.endpoints.size() < MIN_TEMPLATES ||
      primary.nTemplates < MIN_TEMPLATES) {
    return std::unexpected (
        BgGuardFail::InsufficientFragmentData
    );
  }

  // candidate read-length spread.
  if (!read_lens_within_tol (candidate.readLen)) {
    return std::unexpected (
        BgGuardFail::HeterogeneousReadLength
    );
  }

  // read length median vs the primary sample. The size checks above put
  // both counts over their minima, so the primary's medians are engaged.
  const auto candidateReadLenMed =
      median_read_length (candidate.readLen);
  const double readLenRelDiff =
      std::abs (*primary.medianReadLen - *candidateReadLenMed) /
      *primary.medianReadLen;
  if (readLenRelDiff > MEDIAN_REL_TOL) {
    return std::unexpected (BgGuardFail::ReadLengthMismatch);
  }

  // fragment length median vs the primary sample.
  const auto candidateFragLenMed =
      median_fragment_length (candidate.endpoints);
  const double fragLenRelDiff =
      std::abs (*primary.medianFragLen - *candidateFragLenMed) /
      *primary.medianFragLen;
  if (fragLenRelDiff > MEDIAN_REL_TOL) {
    return std::unexpected (BgGuardFail::FragmentLengthMismatch);
  }

  return {};
}
