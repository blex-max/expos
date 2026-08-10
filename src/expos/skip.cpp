#include "skip.hpp"

#include <string>
#include <string_view>

std::string_view to_string (RecordSkipReason reason)
{
  switch (reason) {
    case RecordSkipReason::notBiallelic:
      return "not_biallelic";
    case RecordSkipReason::complex:
      return "complex";
  }
}

std::string_view to_string (StatSkipReason reason)
{
  switch (reason) {
    case StatSkipReason::insufficientSupport:
      return "insufficient_support";
    case StatSkipReason::insufficientBackground:
      return "insufficient_background";
    case StatSkipReason::heterogeneousReadLength:
      return "heterogeneous_read_length";
    case StatSkipReason::readLengthUnverified:
      return "read_length_unverified";
    case StatSkipReason::zeroVariance:
      return "zero_variance";
    case StatSkipReason::noSupport:
      return "no_support";
    case StatSkipReason::noBackground:
      return "no_background";
    case StatSkipReason::referenceTooShort:
      return "reference_too_short";
    case StatSkipReason::referenceHasN:
      return "reference_has_n";
  }
}

Skip make_record_skip (RecordSkipReason reason)
{
  return Skip{SCOPE_RECORD, to_string (reason)};
}

Skip make_stat_skip (
    std::string_view statId, StatSkipReason reason
)
{
  return Skip{statId, to_string (reason)};
}

std::string to_info_format (const Skip& skip)
{
  std::string out;
  out.reserve (skip.scope.size() + 1 + skip.reason.size());
  out += skip.scope;
  out += ':';
  out += skip.reason;
  return out;
}
