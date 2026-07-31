#pragma once

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <string>
#include <vector>

namespace var {

std::vector<bool> bcf_has_filters (
    const bcf_hdr_t* hdr, const bcf1_t* rec,
    const std::vector<std::string>& flt
);

struct VariantRecordFilters {
  std::vector<std::string> must_have;
  std::vector<std::string> must_not_have;
};
bool check_record_filters (
    const bcf1_t* b, const bcf_hdr_t* hdr,
    const VariantRecordFilters& f
);

std::string variant_to_region_str (bcf1_t* v, bcf_hdr_t* h);

}  // namespace var
