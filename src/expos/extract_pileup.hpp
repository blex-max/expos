#pragma once

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "expos/pileup_features.hpp"
#include "expos/vcf_record.hpp"
#include "hts/hts_types.hpp"
#include "shared/err.hpp"

// extract features from pileup, appending to out.
VoidOrErr extract_features (
    PileupView plp, uint16_t maxFragLen, PileupFeatures& _out
);

// extract features to pileup, apppending to outAll,
// and, when a read supports rec, also appending to
// outSupport
VoidOrErr extract_partition_features (
    PileupView plp, const VcfRec& rec, uint16_t maxFragLen,
    PileupFeatures& _outSupport, PileupFeatures& _outAll
);
