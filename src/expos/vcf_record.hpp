#pragma once

#include <htslib/sam.h>
#include <htslib/vcf.h>

#include <cstdint>
#include <expected>
#include <string>
#include <string_view>

#include "expos/skip.hpp"
#include "hts/hts_types.hpp"
#include "shared/err.hpp"

// --- expos-specific VCF record handling --- //

struct VcfRec;

using SupportFn = bool (*) (const VcfRec&, const bam_pileup1_t&);

bool eval_support_snp (const VcfRec& r, const bam_pileup1_t& p);

bool eval_support_mnp (const VcfRec& r, const bam_pileup1_t& p);

bool eval_support_ins (const VcfRec& r, const bam_pileup1_t& p);

bool eval_support_del (const VcfRec& r, const bam_pileup1_t& p);

struct VcfRec {
  bcf1_t* ptr = nullptr;  // owned
  bcf_hdr_t* br_hdr = nullptr;
  SupportFn eval_read_support = nullptr;

  // unpacked data. Valid only until the next bcf_read into ptr.
  std::string_view contig;

  ~VcfRec()
  {
    if (ptr != nullptr) {
      bcf_destroy (ptr);
    }
  }

  VcfRec() = default;
  VcfRec (const VcfRec&) = delete;
  VcfRec (VcfRec&&) = delete;
  VcfRec& operator= (const VcfRec&) = delete;
  VcfRec& operator= (VcfRec&&) = delete;
};

// precondition: r is unpacked; specifically r.contig.
enum class Rec2LocusErr : int {
  notPresent = -1,
  parseFail = -2
};
std::expected<GenomicLocus, Rec2LocusErr> record_to_aln_locus (
    const VcfRec& r, const AlnFile& aln
);

std::string stringify_rec (const VcfRec& r);

// True iff the record failed one or more FILTERs, i.e. FILTER is neither
// PASS nor the missing value '.'. Used by --skip-filtered to pass already-
// rejected records through untouched; PASS and unfiltered ('.') records are
// not "filtered" and are still analysed. PASS is the guaranteed internal
// FILTER id 0 (VCF spec / htslib). Precondition: r is unpacked to at least
// BCF_UN_FLT (covered by the BCF_UN_ALL unpack in the read loop).
bool record_is_filtered (const VcfRec& r);

// Check for any htslib level error state which may render the
// record unwritable, and resolve r.contig from the header.
// Does not unpack; call after bcf_unpack.
VoidOrErr check_integrity (VcfRec& r);

// Whether the record can be analysed, or the reason it can't.
// On success r.eval_read_support is set to the predicate for the
// record's mutation type; on failure it is cleared.
// Precondition: r is unpacked.
std::expected<void, RecordSkipReason> validate_for_analysis (
    VcfRec& r
);
