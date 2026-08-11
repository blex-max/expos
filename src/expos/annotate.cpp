#include "expos/annotate.hpp"

#include <fmt/format.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "expos/compute_info_field.hpp"
#include "expos/encode_info_field.hpp"
#include "expos/extract_pileup.hpp"
#include "expos/guards.hpp"
#include "expos/pileup_features.hpp"
#include "expos/vcf_record.hpp"
#include "hts/hts_types.hpp"
#include "shared/err.hpp"
#include "shared/rng.hpp"
#include "shared/warn.hpp"

// Reference flank either side of the REF allele, in bases, for RCMPLX.
static constexpr int64_t RCMPLX_FLANK = 400;

// Reference bases around the variant, upper-cased: the REF allele span padded
// by RCMPLX_FLANK either side. Deliberately independent of the alignments —
// the span of the supporting templates is defined by the very reads whose
// reliability RCMPLX helps assess, and one mismapped mate stretches it
// arbitrarily (megabases, in practice). Short near a contig edge, in which
// case RCMPLX reports reference_too_short.
static RefSliceOrErr get_variant_ref_slice (
    const VcfRec& r, const FastaFile& ref
)
{
  const int64_t sliceStart =
      std::max<int64_t> (0, r.ptr->pos - RCMPLX_FLANK);
  const int64_t sliceEnd =
      r.ptr->pos + r.ptr->rlen +
      RCMPLX_FLANK;  // faidx clamps to contig end

  auto sliceRet =
      fetch_region (ref, r.contig, sliceStart, sliceEnd);
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

// Merge every admitted background sample's reads at this locus into
// _mergeInto, which arrives holding the primary's reads alone.
static VoidOrErr extract_bg_samples (
    ExposCtx& ctx, const VcfRec& r, PileupFeatures& _mergeInto
)
{
  if (ctx.backgrounds.empty()) {
    return {};
  }

  // Snapshot the primary before the first merge below: it is the fixed
  // point of reference every candidate is judged against, so the outcome
  // does not depend on which candidates were admitted before it. Taking it
  // here rather than at the call site keeps that ordering unbreakable, and
  // keeps a run without --bg from paying for it.
  const PrimaryGuardStats primaryStats =
      summarise_primary (_mergeInto);

  PileupFeatures ru_bg;
  for (auto& bgSamp : ctx.backgrounds) {
    reset (ru_bg);
    const auto bgLocRet = record_to_aln_locus (r, bgSamp.handle);
    if (!bgLocRet) {
      // Both outcomes are fatal for now. Under consideration: hoisting the
      // rid -> tid resolution onto the AlnBundle, at which point notPresent
      // becomes "drop this background for this contig, warn once" rather
      // than killing the run over a background that simply wasn't aligned
      // against this contig. parseFail stays fatal either way.
      switch (bgLocRet.error()) {
        case Rec2LocusErr::notPresent:
          return std::unexpected (make_err (
              fmt::format (
                  "Contig {} (record {}) is not present in the "
                  "background sample at {}.",
                  r.contig, stringify_rec (r), bgSamp.handle.path
              )
          ));
        case Rec2LocusErr::parseFail:
          return std::unexpected (make_err (
              fmt::format (
                  "Could not parse contig {} against the header "
                  "of the background sample at {}.",
                  r.contig, bgSamp.handle.path
              )
          ));
      }
    }
    const auto bgLocus = *bgLocRet;
    const auto bgAdvRet =
        try_advance_pileup (bgSamp.plpIt, bgLocus);
    if (!bgAdvRet) {
      // TODO: stringify error codes
      return std::unexpected (make_err ("failed to pileup!"));
    }
    if (*bgAdvRet == PileupAdvResCode::success) {
      const auto bgExtractRet = extract_features (
          *read (bgSamp.plpIt, bgLocus), ru_bg
      );
      if (!bgExtractRet) {
        return std::unexpected (bgExtractRet.error());
      }
    }
    const auto verifyRet =
        verify_bg_sample (primaryStats, ru_bg);
    if (!verifyRet) {
      if (!ctx.quiet) {
        expos_warn (
            fmt::format (
                "{}: background sample excluded ({})",
                stringify_rec (r), to_string (verifyRet.error())
            )
        );
      }
      continue;
    }
    merge (_mergeInto, ru_bg);
  }
  return {};
}

// Compute and encode the expos statistics onto an analysable record in place.
// nBackgroundExcluded is incremented whenever a background sample is
// excluded from the merge because its read- or fragment-length
// distribution looks inconsistent with the primary sample.
static VoidOrErr annotate_record (
    const VcfRec& r, ExposCtx& ctx, McState& mc
)
{
  // filled by feature extraction for all
  // relevant reads.
  PileupFeatures
      supporting;  // all supporting reads from primary sample
  PileupFeatures all;  // all reads covering the variant locus,
  // including any merged in from bg samples.

  const auto posRet = record_to_aln_locus (r, ctx.aln.handle);
  if (!posRet) {
    // Under consideration alongside the background case above: whether a
    // contig the primary alignment has never heard of should stay fatal or
    // become a record-level skip. Fatal for now.
    switch (posRet.error()) {
      case Rec2LocusErr::notPresent:
        return std::unexpected (make_err (
            fmt::format (
                "Contig {} (record {}) is not present in the "
                "alignment at {}.",
                r.contig, stringify_rec (r), ctx.aln.handle.path
            )
        ));
      case Rec2LocusErr::parseFail:
        return std::unexpected (make_err (
            fmt::format (
                "Could not parse contig {} against the header "
                "of "
                "the alignment at {}.",
                r.contig, ctx.aln.handle.path
            )
        ));
    }
  }
  const auto locus = *posRet;
  const auto advRet = try_advance_pileup (ctx.aln.plpIt, locus);
  if (!advRet) {
    // TODO: stringify error codes
    return std::unexpected (make_err ("failed to pileup!"));
  }
  if (*advRet == PileupAdvResCode::success) {
    const auto extractRet = extract_partition_features (
        *read (ctx.aln.plpIt, locus), r, supporting, all
    );
    if (!extractRet) {
      return std::unexpected (extractRet.error());
    }
  }

  StatContext statCtx{mc, std::nullopt};
  // guard QRK
  // must be done before any background samples merged
  set_qrk_guard (statCtx, all);
  if (statCtx.readLenSuppression && !ctx.quiet) {
    expos_warn (
        fmt::format (
            "{}: primary sample fails the QRK read-length "
            "homogeneity guard "
            "({})",
            stringify_rec (r),
            to_string (*statCtx.readLenSuppression)
        )
    );
  }

  auto bgRet = extract_bg_samples (ctx, r, all);
  if (!bgRet) {
    return std::unexpected (bgRet.error());
  }

  auto sliceRet = get_variant_ref_slice (r, ctx.ref);
  if (!sliceRet) {
    return std::unexpected (sliceRet.error());
  }
  const std::string refSlice = std::move (*sliceRet);

  const VariantStatInputs inputs{supporting, all, refSlice};

  std::vector<Skip> recordSkips;
  for (const auto& stat : expos_field_registry()) {
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
    for (const auto& skip : stat_skips (stat.field, values)) {
      recordSkips.push_back (skip);
    }
  }

  // set skip field, if any skips from compute_*
  if (auto skipRet =
          set_expos_skip (ctx.vcfOut.hdr, r.ptr, recordSkips);
      !skipRet) {
    return std::unexpected (skipRet.error());
  }

  // ready for write
  return {};
}

VoidOrErr analyse_records (ExposCtx& ctx)
{
  // One RNG and set of draw buffers.
  // Buffers should reach high-water mark over
  // the first few records and stop reallocating.
  McState mc{Mwc192 (ctx.seed), {}, {}};

  VcfRec ru_rec;
  ru_rec.ptr = bcf_init();
  ru_rec.br_hdr = ctx.vcfIn.hdr;

  while (bcf_read (ctx.vcfIn.fh, ctx.vcfIn.hdr, ru_rec.ptr) ==
         0) {
    if (bcf_unpack (ru_rec.ptr, BCF_UN_ALL) != 0) {
      return std::unexpected (
          make_err ("Failed to unpack a VCF record")
      );
    }

    const auto integrityRet = check_integrity (ru_rec);
    if (!integrityRet) {
      return std::unexpected (integrityRet.error());
    }

    // Unanalysed records passthrough. filter-gated ones untouched (FILTER
    // already documents them), not-biallelic/complex ones with an
    // EXPOS_SKIP.
    if (ctx.skipFiltered && record_is_filtered (ru_rec)) {
    }
    else {
      const auto valRet = validate_for_analysis (ru_rec);
      if (!valRet) {
        if (!ctx.quiet) {
          expos_warn (
              fmt::format (
                  "Record {} cannot be analysed ({}), passing "
                  "through unannotated",
                  stringify_rec (ru_rec),
                  to_string (valRet.error())
              )
          );
        }
        const Skip skip = make_record_skip (valRet.error());
        const auto skipRet = set_expos_skip (
            ctx.vcfOut.hdr, ru_rec.ptr, {&skip, 1}
        );
        if (!skipRet) {
          return std::unexpected (skipRet.error());
        }
      }
      else {
        const auto annotateRet =
            annotate_record (ru_rec, ctx, mc);
        if (!annotateRet) {
          return std::unexpected (annotateRet.error());
        }
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
    std::cerr << "complete" << std::endl;
  }

  return {};
}
