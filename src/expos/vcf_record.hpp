#pragma once

#include <htslib/sam.h>
#include <htslib/vcf.h>

#include <cstdint>
#include <string>
#include <string_view>

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
  SupportFn eval_read_support;

  // Views into ptr->d.allele[...]; set once per record by
  // type_record. Valid only until the next bcf_read into ptr.
  std::string_view refAllele;
  std::string_view altAllele;
  // insertion bases beyond the shared anchor base (INS only)
  std::string_view insertedBases;
  int indelLen = 0;  // negative for del, positive for ins

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

std::string stringify_rec (const VcfRec& r);

// Whether a record can be analysed, or must be passed through unannotated.
enum class RecordClass : std::uint8_t {
  Analysable,  // biallelic and one of SNP/MNP/INS/DEL
  Unanalysable,  // complex or untypeable
};

// Determine the mutation type and set the read-support predicate, returning
// whether the record is analysable. Complex/untypeable records are
// Unanalysable rather than an error. Precondition: r is biallelic and
// unpacked to at least BCF_UN_STR.
RecordClass type_record (VcfRec& r);

// True iff the record failed one or more FILTERs, i.e. FILTER is neither
// PASS nor the missing value '.'. Used by --skip-filtered to pass already-
// rejected records through untouched; PASS and unfiltered ('.') records are
// not "filtered" and are still analysed. PASS is the guaranteed internal
// FILTER id 0 (VCF spec / htslib). Precondition: r is unpacked to at least
// BCF_UN_FLT (covered by the BCF_UN_ALL unpack in the read loop).
bool record_is_filtered (const VcfRec& r);

// Check for a htslib level error in the input
// vcf which may (??) render the record
// unwritable.
VoidOrErr check_record_integrity (const VcfRec& r);
