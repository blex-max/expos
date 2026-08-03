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

#include "expos/extract_pileup.hpp"
#include "expos/pileup_features.hpp"
#include "expos/variant_stats.hpp"
#include "expos/vcf_record.hpp"
#include "hts/hts_types.hpp"
#include "shared/err.hpp"

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
  const auto* ridName = bcf_hdr_id2name (r.hdr, r.ptr->rid);
  const auto alnTid = sam_hdr_name2tid (aln.o_hdr, ridName);
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

using PileupPosOrErr = std::expected<PileupPosition, Err>;
static PileupPosOrErr make_record_pileup_pos (
    const VcfRec& r, const AlnFile& aln
)
{
  auto tidConvRet = get_aln_contig_for_record_rid (r, aln);
  if (!tidConvRet) {
    return std::unexpected (tidConvRet.error());
  }
  return PileupPosition{*tidConvRet, r.ptr->pos};
}

// Reference bases spanned by the supporting templates, upper-cased. Empty
// when there are no supporting templates (RCMPLX then reports missing).
static RefSliceOrErr supporting_ref_slice (
    const VcfRec& r, const FastaFile& ref,
    const PileupFeatures& supporting
)
{
  if (supporting.endpoints.empty()) {
    return std::string{};
  }
  int64_t spanStart = supporting.endpoints.front().first;
  int64_t spanEnd = supporting.endpoints.front().second;
  for (const auto& [lo, hi] : supporting.endpoints) {
    spanStart = std::min (spanStart, lo);
    spanEnd = std::max (spanEnd, hi);
  }

  const std::string_view contig =
      bcf_hdr_id2name (r.hdr, r.ptr->rid);
  auto sliceRet = fetch_region (ref, contig, spanStart, spanEnd);
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

// Relative IQR of read lengths, (Q3-Q1)/median. Precondition: non-empty.
static double read_length_rel_iqr (
    const std::vector<int32_t>& readLen
)
{
  std::vector<int32_t> lens = readLen;
  std::sort (lens.begin(), lens.end());
  const std::size_t n = lens.size();
  const double q1 = lens[n / 4];
  const double median = lens[n / 2];
  const double q3 = lens[3 * n / 4];
  return (q3 - q1) /
         median;  // median >= 1: read length is always positive
}

// Compute and encode the expos statistics onto an analysable record in place.
static VoidOrErr annotate_record (
    const VcfRec& r, const ExposCtx& ctx, std::mt19937& rng
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

  // Merge any additional background samples into all.
  for (const auto& bg : ctx.backgrounds) {
    const auto bgPosRet = make_record_pileup_pos (r, bg);
    if (!bgPosRet) {
      return std::unexpected (bgPosRet.error());
    }
    const auto bgPlpRet = prepare_pileup (bg, *bgPosRet);
    if (!bgPlpRet) {
      return std::unexpected (bgPlpRet.error());
    }
    const auto bgExtractRet = extract_features (*bgPlpRet, all);
    if (!bgExtractRet) {
      return std::unexpected (bgExtractRet.error());
    }
  }

  // Read-length homogeneity gates QRK
  constexpr std::size_t k_minReadsForHomogeneity = 10;
  constexpr double k_relIqrTol = 0.10;
  bool readLenHomogeneous = true;
  if (all.readLen.size() >= k_minReadsForHomogeneity) {
    const double relIqr = read_length_rel_iqr (all.readLen);
    if (relIqr > k_relIqrTol) {
      readLenHomogeneous = false;
      if (!ctx.quiet) {
        warn (
            std::format (
                "{}: read lengths are heterogeneous (relative "
                "IQR {:.2f} > "
                "{:.2f}); QRK suppressed",
                stringify_rec (r), relIqr, k_relIqrTol
            )
        );
      }
    }
  }

  auto sliceRet = supporting_ref_slice (r, ctx.ref, supporting);
  if (!sliceRet) {
    return std::unexpected (sliceRet.error());
  }
  const std::string refSlice = std::move (*sliceRet);

  const VariantStatInputs inputs{
      supporting, all, refSlice, rng, readLenHomogeneous
  };
  std::vector<std::string> skipTokens;
  for (const auto& stat : variant_stats()) {
    const auto values = stat.compute (inputs);
    const auto encRet = encode_variant_stat (
        ctx.vcfOut.o_hdr, r.ptr, stat.field, values
    );
    if (!encRet) {
      return std::unexpected (encRet.error());
    }
    for (auto& token : stat_skip_tokens (stat.field, values)) {
      skipTokens.push_back (std::move (token));
    }
  }
  if (auto skipRet =
          set_expos_skip (ctx.vcfOut.o_hdr, r.ptr, skipTokens);
      !skipRet) {
    return std::unexpected (skipRet.error());
  }
  return {};
}

VoidOrErr analyse_records (const ExposCtx& ctx)
{
  std::mt19937 rng (ctx.seed);

  VcfRec ru_rec;
  ru_rec.ptr = bcf_init();
  ru_rec.hdr = ctx.vcfIn.o_hdr;

  // Tallies for the end-of-run summary (stderr, unless --quiet).
  std::size_t nRecords = 0;
  std::size_t nSkipped = 0;

  while (
      bcf_read (ctx.vcfIn.o_fh, ctx.vcfIn.o_hdr, ru_rec.ptr) == 0
  ) {
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
          ctx.vcfOut.o_hdr, ru_rec.ptr,
          {"record:" + std::string (*skipReason)}
      );
      if (!skipRet) {
        return std::unexpected (skipRet.error());
      }
    }
    else {
      const auto annotateRet =
          annotate_record (ru_rec, ctx, rng);
      if (!annotateRet) {
        return std::unexpected (annotateRet.error());
      }
    }

    if (bcf_write (
            ctx.vcfOut.o_fh, ctx.vcfOut.o_hdr, ru_rec.ptr
        ) != 0) {
      return std::unexpected (make_err (
          "Failed to write record " + stringify_rec (ru_rec)
      ));
    }
  }

  if (!ctx.quiet) {
    std::cerr << "expos: " << nRecords << " record(s), "
              << nSkipped << " skipped\n";
  }

  return {};
}
