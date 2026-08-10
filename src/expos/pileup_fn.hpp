#pragma once

#include <htslib/sam.h>

#include <cstdint>

#include "hts/hts_types.hpp"

extern "C" {
inline int pileup_fn (void* data, bam1_t* b)
{
  constexpr uint16_t EXCLUDE_BITS = BAM_FSECONDARY | BAM_FDUP |
                                    BAM_FSUPPLEMENTARY |
                                    BAM_FQCFAIL;
  const PileupContext* d = (PileupContext*)(data);

  // find the next good read
  int ret = -1;
  while (true) {
    ret = sam_itr_next (d->br_fh, d->it, b);
    if (ret < 0) {
      break;  // EOF/err
    }
    if ((b->core.flag & EXCLUDE_BITS) == 0) {
      break;  // found good read
    }
  }
  return ret;
}
}
