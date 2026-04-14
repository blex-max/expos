#include <concepts>
#include <cstdint>
#include <format>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <stdexcept>
#include <string>

template <std::signed_integral T>
constexpr inline uint64_t as_uint (const T &a) {
    if (a < 0) {
        throw std::runtime_error (
            "cannot convert negative value to uint"
        );
    }
    return static_cast<uint64_t> (a);
}

inline std::string variant_to_region_str
(bcf1_t* v, bcf_hdr_t* h)
{
    const auto rid_name = bcf_hdr_id2name(h, v->rid);
    return std::format ("{}:{}-{}", rid_name, v->pos + 1, v->pos + 1 + v->rlen);
}
