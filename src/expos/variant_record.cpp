#include "variant_record.hpp"

#include <algorithm>
#include <format>

namespace var {

std::vector<bool>
bcf_has_filters
(const bcf_hdr_t* hdr, const bcf1_t* rec, const std::vector<std::string>& flt)
{
    std::vector<bool> out;
    for (const auto &f : flt) {
        out.push_back (
            bcf_has_filter (hdr, const_cast<bcf1_t*> (rec), const_cast<char *> (f.c_str()))
            > 0
        );
    }
    return out;
}

bool check_record_filters
(const bcf1_t* b, const bcf_hdr_t* hdr, const VariantRecordFilters& f)
{
  const auto iflt =
    bcf_has_filters (hdr, b, f.must_have);
  if (std::any_of (begin (iflt), end (iflt),
        [] (const auto a) {
          return !a;  // any not present
        }
      )) {
    return false;
  }

  const auto eflt =
    bcf_has_filters (hdr, b, f.must_not_have);
  if (std::any_of (begin (eflt), end (eflt),
        [] (const auto a) {
          return a;  // any present
        }
      )) {
      return false;
  }

  return true;
}

std::string variant_to_region_str
(bcf1_t* v, bcf_hdr_t* h)
{
    const auto rid_name = bcf_hdr_id2name(h, v->rid);
    return std::format ("{}:{}-{}", rid_name, v->pos + 1, v->pos + 1 + v->rlen);
}

}
