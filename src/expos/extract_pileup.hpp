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
    const PreparedPileup& plp, PileupFeatures& out
);

// extract features to pileup, apppending to outAll,
// and, when a read supports rec, also appending to
// outSupport
VoidOrErr extract_partition_features (
    const PreparedPileup& plp, const VcfRec& rec,
    PileupFeatures& outSupport, PileupFeatures& outAll
);
