#include <vector>
#include <concepts>
#include <cstdint>
#include <format>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <stdexcept>
#include <string>

template <std::signed_integral T>
uint64_t
constexpr inline as_uint (const T &a)
{
    if (a < 0) {
        throw std::runtime_error (
            "cannot convert negative value to uint"
        );
    }
    return static_cast<uint64_t> (a);
}

std::string
inline variant_to_region_str
(bcf1_t* v, bcf_hdr_t* h)
{
    const auto rid_name = bcf_hdr_id2name(h, v->rid);
    return std::format ("{}:{}-{}", rid_name, v->pos + 1, v->pos + 1 + v->rlen);
}

std::vector<bool>
inline bcf_has_filters (
    bcf_hdr_t                *hdr,
    bcf1_t                   *rec,
    std::vector<std::string> &flt
)
{
    std::vector<bool> out;
    for (const auto &f : flt) {
        out.push_back (
            bcf_has_filter (hdr, rec, const_cast<char *> (f.c_str()))
            > 0
        );
    }
    return out;
}
