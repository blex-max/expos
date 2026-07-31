#pragma once

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <cstdint>
#include <functional>
#include <vector>

#include "lib-stats/spatial.hpp"
#include "util.hpp"

struct PileupMetrics {
  std::vector<uint64_t> query_position;
  std::vector<double> normalised_as;
  std::vector<spatial::line_seg> template_endpoints;
  std::vector<uint32_t> read_lengths;
};
PileupMetrics merge_pileup_metrics (
    PileupMetrics a, const PileupMetrics& b
);

PileupMetrics pileup_analyse (
    htsFile* aln_fh, hts_idx_t* aln_idx, bcf1_t* var,
    int sam_flag_include, int sam_flag_exclude,
    const bcf_hdr_t* vcf_hdr
);

struct GroupedMetrics {
  PileupMetrics supporting;
  PileupMetrics all;
};
GroupedMetrics pileup_group_and_anaylse (
    htsFile* aln_fh, hts_idx_t* aln_idx, bcf1_t* var,
    int sam_flag_include, int sam_flag_exclude,
    std::function<bool (const bcf1_t*, const bam_pileup1_t*)>
        support_fn,
    const bcf_hdr_t* vcf_hdr
);
