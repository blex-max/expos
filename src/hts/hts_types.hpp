#pragma once

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include <cstdint>
#include <expected>
#include <optional>
#include <span>
#include <string>

#include "shared/err.hpp"

struct GenomicSpan {
  hts_pos_t start;
  hts_pos_t end;
};
struct GenomicLocus {
  int32_t tid;
  hts_pos_t pos;
};

inline bool operator== (
    const GenomicLocus& a, const GenomicLocus& b
)
{
  return a.tid == b.tid && a.pos == b.pos;
}

struct AlnFile {
  htsFile* fh = NULL;
  sam_hdr_t* hdr = NULL;
  hts_idx_t* idx = NULL;
  std::string path;

  AlnFile() = default;

  // no copy
  AlnFile (const AlnFile&) = delete;
  AlnFile& operator= (const AlnFile&) = delete;

  // move ok
  AlnFile (AlnFile&& o) noexcept;
  AlnFile& operator= (AlnFile&& o) noexcept;

  ~AlnFile() noexcept;
};

using AlnOrErr = std::expected<AlnFile, Err>;
AlnOrErr load_aln (const std::string& path);

struct FastaFile {
  faidx_t* fai = nullptr;
  std::string path;

  operator faidx_t*() const noexcept { return fai; }

  FastaFile() = default;

  // no copy
  FastaFile (const FastaFile& o) = delete;
  FastaFile& operator= (const FastaFile& o) = delete;

  // move ok
  FastaFile (FastaFile&& o) noexcept;
  FastaFile& operator= (FastaFile&& o) noexcept;

  ~FastaFile();
};

using FastaOrErr = std::expected<FastaFile, Err>;
FastaOrErr load_fasta (const std::string& path);

using RefSliceOrErr = std::expected<std::string, Err>;
RefSliceOrErr fetch_region (
    const FastaFile& ff, const std::string_view contigName,
    hts_pos_t regStart, hts_pos_t regEnd
);

struct VcfFile {
  htsFile* fh = NULL;
  bcf_hdr_t* hdr = NULL;
  std::string path;

  VcfFile() = default;

  // no copy
  VcfFile (const VcfFile&) = delete;
  VcfFile& operator= (const VcfFile&) = delete;

  // move ok
  VcfFile (VcfFile&& o) noexcept;
  VcfFile& operator= (VcfFile&& o) noexcept;

  ~VcfFile() noexcept;
};

using VcfOrErr = std::expected<VcfFile, Err>;
VcfOrErr load_vcf (const std::string& path);

// pileup iterator type:
// stateful forward-per-contig iterator
// for retrieving pileups in a performant
// manner - only rebuilds iterator at large
// gaps between requested positions.

// Forward-only iteration: In context of expos,
// a request behind the previous one on the same
// contig means the caller's records are not
// coordinate-sorted. I.e. this mandates
// coordinate-sorted VCF (which the VCF spec
// mandates anyway).

// Context object recieved by bam_plp_auto_f (the
// htslib pileup function type) as a void*.
struct PileupContext {
  htsFile* br_fh;  // br_ = borrowed, must not close
  hts_idx_t* br_fhIdx;
  hts_itr_t* it = nullptr;  // owned, must destroy

  explicit PileupContext (const AlnFile& aln);

  PileupContext (const PileupContext&) = delete;
  PileupContext& operator= (const PileupContext&) = delete;

  // virtual: ForwardPileupIterator destroys the context
  // through this type, so a derived context must get its
  // own destructor run.
  virtual ~PileupContext();
};

struct ForwardPileupIterator {
  // The cost of a seek is ~fixed,
  // the cost of streaming through reads
  // scales linearly with sample depth.
  // ~4000-8000 was found to be optimal on a 100x
  // depth BAM. Since streaming should be ~2x
  // cheaper on a 50x bam (closer to intended
  // samples for expos) set at 10000. May be
  // suboptimal but it should be close enough.
  static constexpr hts_pos_t SEEK_GAP = 10000;

  struct {
    // View into plpIt's buffer; invalidated by the next
    // bam_plp64_auto or bam_plp_reset on plpIt.
    const bam_pileup1_t* plpArr = nullptr;
    uint32_t nReads = 0;
    GenomicLocus locus{
        -1, -1
    };  // current position. Start invalid.
  } curs;
  GenomicLocus lastRequest{-1, -1};

  PileupContext* ctx =
      nullptr;  // runs own destructor. ptr as requires stable address.
  bam_plp_s* plpIt =
      nullptr;  // owned, must destroy on destruction

  ~ForwardPileupIterator();

  // no copy
  ForwardPileupIterator (const ForwardPileupIterator&) = delete;
  ForwardPileupIterator& operator= (
      const ForwardPileupIterator&
  ) = delete;

  // move ok
  ForwardPileupIterator (ForwardPileupIterator&&) noexcept;
  ForwardPileupIterator& operator= (
      ForwardPileupIterator&&
  ) noexcept;

  ForwardPileupIterator() = default;
};

// Takes ownership of `snk_ctx`, which therefore should
// probably be heap allocated. snk_ctx may be any derived
// type of PileupContext. Using a derived type allows
// additional context fields. If using a derived type,
// `fn` must cast in two steps, first to PileupContext*,
// and then to the derived type ptr.
ForwardPileupIterator init_pileup_iterator (
    PileupContext* snk_ctx, bam_plp_auto_f fn
);

enum class PileupAdvResCode : uint8_t {
  success,
  exhausted,
  noCoverage
};
enum class PileupAdvErrCode : uint8_t {
  invalidLocus,
  locusBehindItr,
  samItrFailed,
  pileupEngineFailed
};
// Move cursor to `to`.
// Separate to accessing the result as the
// outcomes are quite rich.
using AdvResOrErr =
    std::expected<PileupAdvResCode, PileupAdvErrCode>;
[[nodiscard]] AdvResOrErr try_advance_pileup (
    ForwardPileupIterator& pc, const GenomicLocus& to
);

// Non-owning view of the pileup the cursor is currently
// parked on.
//
// unexpected(nullopt) == the cursor is not parked on
// expectedLocus, which is the case after any
// try_advance_pileup that did not return `success`: a
// noCoverage park sits on the next covered position, so its
// span is typically non-empty but belongs to a different
// locus; an exhausted park sits on the end sentinel. Only a
// `success` guarantees a view of the locus asked for.
//
// The span aliases plpIt's internal buffer, and so is
// invalidated by the next try_advance_pileup on the same
// iterator.
using PileupView = std::span<const bam_pileup1_t>;
using PileupViewOrNone =
    std::expected<PileupView, std::nullopt_t>;
PileupViewOrNone read (
    const ForwardPileupIterator& pc,
    const GenomicLocus& expectedLocus
);
