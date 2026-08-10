#pragma once

#include <htslib/vcf.h>

#include <span>
#include <string_view>
#include <vector>

#include "expos/compute_info_field.hpp"
#include "expos/skip.hpp"
#include "shared/err.hpp"

// --- rendering computed statistics into VCF INFO fields --- //

// The only layer that talks to htslib's VCF writer. Keeping it separate is
// what lets compute_info_field.hpp stay free of <htslib/vcf.h>, so the
// statistics can be computed and tested without a VCF in play.

// Append the ##INFO header line for one statistic.
VoidOrErr register_variant_stat_header (
    bcf_hdr_t* hdr, const VariantStat& stat
);

// Append the ##INFO header line for the EXPOS_SKIP field.
VoidOrErr register_expos_skip_header (bcf_hdr_t* hdr);

// Expand a whole-statistic skip into one missing StatValue per subfield.
std::vector<StatValue> stat_all_missing (
    const VariantStatField& field, StatSkipReason reason
);

// Write a statistic's per-subfield values to `rec` as its Float INFO field.
// Missing subfields become VCF missing values. `values.size()` must equal
// `field.nValues`.
VoidOrErr encode_variant_stat (
    bcf_hdr_t* hdr, bcf1_t* rec, const VariantStatField& field,
    const std::vector<StatValue>& values
);

// Deduplicated Skips for the missing subfields of `values`, scoped to the
// statistic. The dedupe is what collapses a whole-statistic skip to a single
// entry, since every subfield of a stat_all_missing() result carries the
// same reason.
std::vector<Skip> stat_skips (
    const VariantStatField& field,
    const std::vector<StatValue>& values
);

// Write the EXPOS_SKIP INFO field from a set of skips.
// No-op when skips is empty.
VoidOrErr set_expos_skip (
    bcf_hdr_t* hdr, bcf1_t* rec, std::span<const Skip> skips
);
