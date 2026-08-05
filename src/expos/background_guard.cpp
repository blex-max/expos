#include "background_guard.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <ranges>
#include <vector>

#include "shared/stats.hpp"

// Relative IQR of read lengths, (Q3-Q1)/median. Precondition: non-empty.
static double read_length_rel_iqr (
    const std::vector<int32_t>& readLen
)
{
  std::vector<int32_t> lens = readLen;
  std::sort (lens.begin(), lens.end());
  const std::size_t n = lens.size();
  const double q1 = lens[n / 4];
  const double median = lens[n / 2];
  const double q3 = lens[3 * n / 4];
  return (q3 - q1) /
         median;  // median >= 1: read length is always positive
}

constexpr double k_median = 0.5;

static std::optional<double> median_read_length (
    const std::vector<int32_t>& readLen
)
{
  return percentile (
      std::views::transform (
          readLen,
          [] (int32_t len) {
            return static_cast<std::size_t> (len);
          }
      ),
      k_median
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
            return static_cast<std::size_t> (
                span.second - span.first
            );
          }
      ),
      k_median
  );
}

std::string_view to_string (BackgroundGuardReason reason)
{
  switch (reason) {
    case BackgroundGuardReason::Admitted:
      return "admitted";
    case BackgroundGuardReason::InsufficientReadData:
      return "insufficient reads to check read-length "
             "consistency";
    case BackgroundGuardReason::InsufficientFragmentData:
      return "insufficient templates to check fragment-length "
             "consistency";
    case BackgroundGuardReason::HeterogeneousOwnReadLength:
      return "its own read lengths are too heterogeneous";
    case BackgroundGuardReason::ReadLengthMismatch:
      return "its median read length differs from the primary "
             "sample's";
    case BackgroundGuardReason::FragmentLengthMismatch:
      return "its median fragment length differs from the "
             "primary "
             "sample's";
  }
  return "unknown";
}

bool primary_read_length_homogeneous (
    const PileupFeatures& primaryAll, std::size_t minReads,
    double relIqrTol
)
{
  if (primaryAll.readLen.size() < minReads) {
    // Insufficient data: assume homogeneous (unchanged from prior
    // behaviour -- there is no background sample to fall back on here).
    return true;
  }
  return read_length_rel_iqr (primaryAll.readLen) <= relIqrTol;
}

BackgroundGuardReason evaluate_background (
    const PileupFeatures& primaryAll,
    const PileupFeatures& candidate, std::size_t minReads,
    double relIqrTol, std::size_t minTemplates,
    double medianRelTol
)
{
  if (candidate.readLen.size() < minReads ||
      primaryAll.readLen.size() < minReads) {
    return BackgroundGuardReason::InsufficientReadData;
  }
  // Guard B: the candidate's own read-length spread.
  if (read_length_rel_iqr (candidate.readLen) > relIqrTol) {
    return BackgroundGuardReason::HeterogeneousOwnReadLength;
  }
  // Guard C, read length: median vs the primary sample's.
  const auto primaryReadLenMed =
      median_read_length (primaryAll.readLen);
  const auto candidateReadLenMed =
      median_read_length (candidate.readLen);
  const double readLenRelDiff =
      std::abs (*primaryReadLenMed - *candidateReadLenMed) /
      *primaryReadLenMed;
  if (readLenRelDiff > medianRelTol) {
    return BackgroundGuardReason::ReadLengthMismatch;
  }

  if (candidate.endpoints.size() < minTemplates ||
      primaryAll.endpoints.size() < minTemplates) {
    return BackgroundGuardReason::InsufficientFragmentData;
  }
  // Guard C, fragment length: median vs the primary sample's.
  const auto primaryFragLenMed =
      median_fragment_length (primaryAll.endpoints);
  const auto candidateFragLenMed =
      median_fragment_length (candidate.endpoints);
  const double fragLenRelDiff =
      std::abs (*primaryFragLenMed - *candidateFragLenMed) /
      *primaryFragLenMed;
  if (fragLenRelDiff > medianRelTol) {
    return BackgroundGuardReason::FragmentLengthMismatch;
  }

  return BackgroundGuardReason::Admitted;
}
