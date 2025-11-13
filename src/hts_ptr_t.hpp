#pragma once

#include <memory>

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>


using htsFile_upt = std::unique_ptr<htsFile, decltype (&hts_close)>;
using hts_idx_upt =
    std::unique_ptr<hts_idx_t, decltype (&hts_idx_destroy)>;
using vrec_upt = std::unique_ptr<bcf1_t, decltype (&bcf_destroy)>;
using bcf_hdr_upt =
    std::unique_ptr<bcf_hdr_t, decltype (&bcf_hdr_destroy)>;
using fai_upt = std::unique_ptr<faidx_t, decltype (&fai_destroy)>;
using aln_rec_upt = std::unique_ptr<bam1_t, decltype (&bam_destroy1)>;
using sam_hdr_upt =
    std::unique_ptr<sam_hdr_t, decltype (&sam_hdr_destroy)>;
// using bplp1_upt = std::unique_ptr<bam_pileup1_t, decltype (&bam_)
