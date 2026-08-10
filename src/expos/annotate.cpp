#include "expos/annotate.hpp"

#include <htslib/sam.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "expos/background_guard.hpp"
#include "expos/compute_field.hpp"
#include "expos/extract_pileup.hpp"
#include "expos/pileup_features.hpp"
#include "expos/variant_stats.hpp"
#include "expos/vcf_record.hpp"
#include "hts/hts_types.hpp"
#include "shared/err.hpp"

static constexpr std::size_t MIN_READS_FOR_LEN_CHECK = 10;
static constexpr std::size_t MIN_TEMPLATES_FOR_LEN_CHECK = 10;
static constexpr double READ_LEN_REL_IQR_TOL = 0.10;
static constexpr double MEDIAN_REL_TOL = 0.10;
// Reference flank either side of the REF allele, in bases, for RCMPLX.
static constexpr int64_t RCMPLX_FLANK = 400;

// User-facing warning
static void warn (const std::string& msg)
{
  std::cerr << "warning: " << msg << '\n';
}

// The reason a record can't be analysed
// (multiallelic or complex), or nullopt.
// Precondition: r is unpacked.
static std::optional<std::string_view> classify_record (
    VcfRec& r, bool quiet
)
{
  if (r.ptr->n_allele != 2) {
    if (!quiet) {
      warn (
          stringify_rec (r) + " is not biallelic (" +
          std::to_string (r.ptr->n_allele) +
          " alleles); passing through unannotated"
      );
    }
    return REASON_MULTIALLELIC;
  }
  if (type_record (r) == RecordClass::Unanalysable) {
    if (!quiet) {
      warn (
          stringify_rec (r) +
          " is complex or untypeable; passing through "
          "unannotated"
      );
    }
    return REASON_COMPLEX;
  }
  return std::nullopt;
}

// Resolve the record's contig to the alignment file's tid.
// sam_hdr_name2tid: -1 contig not present in alignment, -2 parse failure.
static IntOrErr get_aln_contig_for_record_rid (
    const VcfRec& r, const AlnFile& aln
)
{
  const auto* ridName = bcf_hdr_id2name (r.br_hdr, r.ptr->rid);
  const auto alnTid = sam_hdr_name2tid (aln.hdr, ridName);
  if (alnTid == -1) {
    std::string errMsg = "Contig ";
    errMsg += ridName;
    errMsg += " for variant ";
    errMsg += stringify_rec (r);
    errMsg += " not found";  // TODO: in alignment file at ...
    return std::unexpected (make_err (errMsg));
  }
  if (alnTid == -2) {
    std::string errMsg = "Unable to parse contig ";
    errMsg += ridName;
    errMsg += " for variant ";
    errMsg += stringify_rec (r);
    errMsg +=
        " as contig name for SAM file";  // TODO in alignment file at ...
    return std::unexpected (make_err (errMsg));
  }
  return alnTid;
}

using LocusOrErr = std::expected<GenomicLocus, Err>;
static LocusOrErr make_record_pileup_pos (
    const VcfRec& r, const AlnFile& aln
)
{
  auto tidConvRet = get_aln_contig_for_record_rid (r, aln);
  if (!tidConvRet) {
    return std::unexpected (tidConvRet.error());
  }
  return GenomicLocus{*tidConvRet, r.ptr->pos};
}

// Reference bases around the variant, upper-cased: the REF allele span padded
// by RCMPLX_FLANK either side. Deliberately independent of the alignments —
// the span of the supporting templates is defined by the very reads whose
// reliability RCMPLX helps assess, and one mismapped mate stretches it
// arbitrarily (megabases, in practice). Short near a contig edge, in which
// case RCMPLX reports reference_too_short.
static RefSliceOrErr variant_ref_slice (
    const VcfRec& r, const FastaFile& ref
)
{
  const int64_t sliceStart =
      std::max<int64_t> (0, r.ptr->pos - RCMPLX_FLANK);
  const int64_t sliceEnd =
      r.ptr->pos + r.ptr->rlen +
      RCMPLX_FLANK;  // faidx clamps to contig end

  const std::string_view contig =
      bcf_hdr_id2name (r.br_hdr, r.ptr->rid);
  auto sliceRet =
      fetch_region (ref, contig, sliceStart, sliceEnd);
  if (!sliceRet) {
    return std::unexpected (sliceRet.error());
  }
  std::string slice = std::move (*sliceRet);
  std::transform (
      slice.begin(), slice.end(), slice.begin(),
      [] (unsigned char c) {
        return static_cast<char> (std::toupper (c));
      }
  );
  return slice;
}

// Compute and encode the expos statistics onto an analysable record in place.
// nBackgroundExcluded is incremented whenever a background sample is
// excluded from the merge because its read- or fragment-length
// distribution looks inconsistent with the primary sample's.
static VoidOrErr annotate_record (
    const VcfRec& r, const ExposCtx& ctx, McState& mc,
    std::size_t& nBackgroundExcluded
)
{
  const auto posRet = make_record_pileup_pos (r, ctx.aln);
  if (!posRet) {
    return std::unexpected (posRet.error());
  }
  auto prepRet = prepare_pileup (ctx.aln, *posRet);
  if (!prepRet) {
    return std::unexpected (prepRet.error());
  }
  const auto plp = std::move (*prepRet);

  PileupFeatures supporting;
  PileupFeatures all;
  const auto extractRet =
      extract_partition_features (plp, r, supporting, all);
  if (!extractRet) {
    return std::unexpected (extractRet.error());
  }

  // guard QRK
  const bool readLenHomogeneous =
      primary_read_length_homogeneous (
          all, MIN_READS_FOR_LEN_CHECK, READ_LEN_REL_IQR_TOL
      );
  if (!readLenHomogeneous && !ctx.quiet) {
    warn (
        std::format (
            "{}: primary sample's read lengths are "
            "heterogeneous; QRK suppressed",
            stringify_rec (r)
        )
    );
  }

  // evaluate each background source against the primary
  // before merging.
  PileupFeatures ru_bg;
  for (const auto& bg : ctx.backgrounds) {
    reset (ru_bg);
    const auto bgPosRet = make_record_pileup_pos (r, bg);
    if (!bgPosRet) {
      return std::unexpected (bgPosRet.error());
    }
    const auto bgPlpRet = prepare_pileup (bg, *bgPosRet);
    if (!bgPlpRet) {
      return std::unexpected (bgPlpRet.error());
    }
    const auto bgExtractRet =
        extract_features (*bgPlpRet, ru_bg);
    if (!bgExtractRet) {
      return std::unexpected (bgExtractRet.error());
    }
    const auto reason = evaluate_background (
        all, ru_bg, MIN_READS_FOR_LEN_CHECK,
        READ_LEN_REL_IQR_TOL, MIN_TEMPLATES_FOR_LEN_CHECK,
        MEDIAN_REL_TOL
    );
    if (reason != BackgroundGuardReason::Admitted) {
      ++nBackgroundExcluded;
      if (!ctx.quiet) {
        warn (
            std::format (
                "{}: background sample excluded ({})",
                stringify_rec (r), to_string (reason)
            )
        );
      }
      continue;
    }
    merge (all, ru_bg);
  }

  auto sliceRet = variant_ref_slice (r, ctx.ref);
  if (!sliceRet) {
    return std::unexpected (sliceRet.error());
  }
  const std::string refSlice = std::move (*sliceRet);

  const VariantStatInputs inputs{supporting, all, refSlice};
  const StatContext statCtx{mc, readLenHomogeneous};
  std::vector<std::string> skipTokens;
  for (const auto& stat : variant_stats()) {
    auto result = stat.compute (inputs, statCtx);
    const std::vector<StatValue> values =
        result ? std::move (*result)
               : stat_all_missing (stat.field, result.error());
    const auto encRet = encode_variant_stat (
        ctx.vcfOut.hdr, r.ptr, stat.field, values
    );
    if (!encRet) {
      return std::unexpected (encRet.error());
    }
    for (auto& token : stat_skip_tokens (stat.field, values)) {
      skipTokens.push_back (std::move (token));
    }
  }

  // set skip field, if any skips from compute_*
  if (auto skipRet =
          set_expos_skip (ctx.vcfOut.hdr, r.ptr, skipTokens);
      !skipRet) {
    return std::unexpected (skipRet.error());
  }

  return {};
}

VoidOrErr analyse_records (const ExposCtx& ctx)
{
  // One RNG and set of draw buffers.
  // Buffers should reach high-water mark over
  // the first few records and stop reallocating.
  McState mc{std::mt19937 (ctx.seed), {}, {}};

  VcfRec ru_rec;
  ru_rec.ptr = bcf_init();
  ru_rec.br_hdr = ctx.vcfIn.hdr;

  // Tallies for the end-of-run summary (stderr, unless --quiet).
  std::size_t nRecords = 0;
  std::size_t nSkipped = 0;
  std::size_t nBackgroundExcluded = 0;

  while (bcf_read (ctx.vcfIn.fh, ctx.vcfIn.hdr, ru_rec.ptr) ==
         0) {
    if (bcf_unpack (ru_rec.ptr, BCF_UN_ALL) != 0) {
      return std::unexpected (
          make_err ("Failed to unpack a VCF record")
      );
    }

    const auto integrityRet = check_record_integrity (ru_rec);
    if (!integrityRet) {
      return std::unexpected (integrityRet.error());
    }

    ++nRecords;

    // Unanalysed records passthrough. filter-gated ones untouched (FILTER
    // already documents them), multiallelic/complex ones with an EXPOS_SKIP.
    if (ctx.skipFiltered && record_is_filtered (ru_rec)) {
      ++nSkipped;
    }
    else if (const auto skipReason =
                 classify_record (ru_rec, ctx.quiet)) {
      ++nSkipped;
      const auto skipRet = set_expos_skip (
          ctx.vcfOut.hdr, ru_rec.ptr,
          {"record:" + std::string (*skipReason)}
      );
      if (!skipRet) {
        return std::unexpected (skipRet.error());
      }
    }
    else {
      const auto annotateRet =
          annotate_record (ru_rec, ctx, mc, nBackgroundExcluded);
      if (!annotateRet) {
        return std::unexpected (annotateRet.error());
      }
    }

    if (bcf_write (ctx.vcfOut.fh, ctx.vcfOut.hdr, ru_rec.ptr) !=
        0) {
      return std::unexpected (make_err (
          "Failed to write record " + stringify_rec (ru_rec)
      ));
    }
  }

  if (!ctx.quiet) {
    std::cerr << "expos: " << nRecords << " record(s), "
              << nSkipped << " skipped, " << nBackgroundExcluded
              << " background sample(s) excluded\n";
  }

  return {};
}
