#pragma once

#include <cstdint>
#include <format>
#include <htslib/faidx.h>
#include <memory>
#include <string>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>


using htsFile_upt = std::unique_ptr<htsFile, decltype (&hts_close)>;

using hts_idx_upt = std::
    unique_ptr<hts_idx_t, decltype (&hts_idx_destroy)>;
using hts_itr_upt = std::
    unique_ptr<hts_itr_t, decltype (&hts_itr_destroy)>;

using bcf1_upt = std::unique_ptr<bcf1_t, decltype (&bcf_destroy)>;
using bam1_upt = std::unique_ptr<bam1_t, decltype (&bam_destroy1)>;

using sam_hdr_upt = std::
    unique_ptr<sam_hdr_t, decltype (&sam_hdr_destroy)>;
using bcf_hdr_upt = std::
    unique_ptr<bcf_hdr_t, decltype (&bcf_hdr_destroy)>;

// may need tweaking, plp is weird
using bam_plp_upt = std::
    unique_ptr<bam_plp_s, decltype (&bam_plp_destroy)>;

using fai_upt = std::unique_ptr<faidx_t, decltype (&fai_destroy)>;


inline std::string fai_autofetch (
    faidx_t    *f,
    const char *rid_name,
    uint64_t    start,
    uint64_t    end
) {
    if (!f) {
        throw std::runtime_error ("fetch failed, nullptr");
    }
    int64_t out_len = 0;
    char   *p       = faidx_fetch_seq64 (
        f,
        rid_name,
        static_cast<int64_t> (start),
        static_cast<int64_t> (end),
        &out_len
    );
    if (!p || out_len < 0)
        throw std::runtime_error (
            std::format (
                "fetch failed for region {}:{}-{}",
                rid_name,
                start,
                end
            )
        );

    std::string s (p, static_cast<size_t> (out_len));
    free (p);
    return s;
}
