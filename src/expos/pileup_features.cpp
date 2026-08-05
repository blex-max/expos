#include "pileup_features.hpp"

#include <iterator>

void reset (PileupFeatures& pf) noexcept
{
  pf.qPos.clear();
  pf.endpoints.clear();
  pf.readLen.clear();
  pf.normalisedAs.clear();
}

void merge (PileupFeatures& into, PileupFeatures& from) noexcept
{
  into.qPos.insert (
      into.qPos.end(),
      std::make_move_iterator (from.qPos.begin()),
      std::make_move_iterator (from.qPos.end())
  );
  into.endpoints.insert (
      into.endpoints.end(),
      std::make_move_iterator (from.endpoints.begin()),
      std::make_move_iterator (from.endpoints.end())
  );
  into.readLen.insert (
      into.readLen.end(),
      std::make_move_iterator (from.readLen.begin()),
      std::make_move_iterator (from.readLen.end())
  );
  into.normalisedAs.insert (
      into.normalisedAs.end(),
      std::make_move_iterator (from.normalisedAs.begin()),
      std::make_move_iterator (from.normalisedAs.end())
  );
}
