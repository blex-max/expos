#pragma once

#include <memory>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>


using htsFile_upt = std::unique_ptr<htsFile, decltype (&hts_close)>;

using hts_idx_upt = std::unique_ptr<hts_idx_t, decltype (&hts_idx_destroy)>;
using hts_itr_upt = std::unique_ptr<hts_itr_t, decltype (&hts_itr_destroy)>;

using bcf1_upt = std::unique_ptr<bcf1_t, decltype (&bcf_destroy)>;
using bam1_upt = std::unique_ptr<bam1_t, decltype (&bam_destroy1)>;

using sam_hdr_upt = std::unique_ptr<sam_hdr_t, decltype (&sam_hdr_destroy)>;
using bcf_hdr_upt = std::unique_ptr<bcf_hdr_t, decltype (&bcf_hdr_destroy)>;

// may need tweaking, plp is weird
using bam_plp_upt = std::unique_ptr<bam_plp_s, decltype (&bam_plp_destroy)>;
