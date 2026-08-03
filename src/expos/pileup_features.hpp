#pragma once

#include <cstdint>
#include <utility>
#include <vector>

// Per-read / per-template features extracted from a pileup column. Pure
// data with no htslib dependency — the interface between pileup extraction
// (extract_pileup) and the variant statistics (variant_stats).

// Leftmost/rightmost genomic coordinates of a sequenced template.
using TemplateEndpoints = std::pair<int64_t, int64_t>;

struct PileupFeatures {
  std::vector<int32_t> qPos;  // per aligned read (qpos >= 0)
  std::vector<TemplateEndpoints>
      endpoints;  // per template (deduplicated)
  std::vector<int32_t> readLen;  // per read
  std::vector<double> normalisedAs;  // per read: AS tag / l_qseq
};
