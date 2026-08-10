#pragma once

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include <expected>
#include <functional>
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

// resolve tid to name
using Tid2StrFn = std::function<const char*(int)>;

struct AlnFile {
  htsFile* fh = NULL;
  sam_hdr_t* hdr = NULL;
  hts_idx_t* idx = NULL;

  AlnFile() = default;
  AlnFile (const AlnFile&) = delete;
  AlnFile& operator= (const AlnFile&) = delete;

  ~AlnFile() noexcept
  {
    if (idx != nullptr) {
      hts_idx_destroy (idx);
    }
    if (hdr != nullptr) {
      sam_hdr_destroy (hdr);
    }
    if (fh != nullptr) {
      hts_close (fh);
    }
  }

  AlnFile (AlnFile&& o) noexcept
      : fh{o.fh}, hdr{o.hdr}, idx{o.idx}
  {
    o.fh = NULL;
    o.hdr = NULL;
    o.idx = NULL;
  }

  AlnFile& operator= (AlnFile&& o) noexcept
  {
    if (this != &o) {
      if (idx != nullptr) {
        hts_idx_destroy (idx);
      }
      if (hdr != nullptr) {
        sam_hdr_destroy (hdr);
      }
      if (fh != nullptr) {
        hts_close (fh);
      }
      fh = o.fh;
      hdr = o.hdr;
      idx = o.idx;
      o.fh = NULL;
      o.hdr = NULL;
      o.idx = NULL;
    }
    return *this;
  }
};

using AlnOrErr = std::expected<AlnFile, Err>;
AlnOrErr load_aln (const char* fn);

struct FastaFile {
  faidx_t* fai = nullptr;

  operator faidx_t*() const noexcept { return fai; }

  FastaFile() = default;
  FastaFile (const FastaFile& o) = delete;
  FastaFile& operator= (const FastaFile& o) = delete;

  FastaFile (FastaFile&& o) noexcept : fai{o.fai}
  {
    o.fai = nullptr;
  }

  FastaFile& operator= (FastaFile&& o) noexcept
  {
    if (this != &o) {
      if (fai != nullptr) {
        fai_destroy (fai);
      }
      fai = o.fai;
      o.fai = nullptr;
    }
    return *this;
  }

  ~FastaFile()
  {
    if (fai != nullptr) {
      fai_destroy (fai);
    }
  }
};

using FastaOrErr = std::expected<FastaFile, Err>;
FastaOrErr load_fasta (const char* fn);

using RefSliceOrErr = std::expected<std::string, Err>;
RefSliceOrErr fetch_region (
    const FastaFile& ff, const std::string_view contigName,
    hts_pos_t regStart, hts_pos_t regEnd
);

struct VcfFile {
  htsFile* fh = NULL;
  bcf_hdr_t* hdr = NULL;

  VcfFile() = default;
  VcfFile (const VcfFile&) = delete;
  VcfFile& operator= (const VcfFile&) = delete;

  ~VcfFile() noexcept
  {
    if (hdr != nullptr) {
      bcf_hdr_destroy (hdr);
    }
    if (fh != nullptr) {
      hts_close (fh);
    }
  }

  VcfFile (VcfFile&& o) noexcept : fh{o.fh}, hdr{o.hdr}
  {
    o.fh = NULL;
    o.hdr = NULL;
  }

  VcfFile& operator= (VcfFile&& o) noexcept
  {
    if (this != &o) {
      if (hdr != nullptr) {
        bcf_hdr_destroy (hdr);
      }
      if (fh != nullptr) {
        hts_close (fh);
      }
      fh = o.fh;
      hdr = o.hdr;
      o.fh = NULL;
      o.hdr = NULL;
    }
    return *this;
  }
};

using VcfOrErr = std::expected<VcfFile, Err>;
VcfOrErr load_vcf (const char* fn);

struct PreparedPileup {
  bam_plp_t plpBacking = nullptr;   // owned
  // View into plpBacking's buffer; invalidated by the next
  // bam_plp64_auto or bam_plp_reset on plpBacking.
  const bam_pileup1_t* plpArr = nullptr;
  size_t nPlp = 0;

  ~PreparedPileup()
  {
    if (plpBacking != nullptr) {
      bam_plp_destroy (plpBacking);
    }
    plpArr = nullptr;
  }
  PreparedPileup() = default;
  PreparedPileup (PreparedPileup&) = delete;
  PreparedPileup& operator= (PreparedPileup&) = delete;
  PreparedPileup (PreparedPileup&& o) noexcept
      : plpBacking (o.plpBacking),
        plpArr (o.plpArr),
        nPlp (o.nPlp)
  {
    o.plpBacking = nullptr;
    o.plpArr = nullptr;
    o.nPlp = 0;
  };
  PreparedPileup& operator= (PreparedPileup&&) = delete;
};

using PileupOrErr = std::expected<PreparedPileup, Err>;
PileupOrErr prepare_pileup (
    const AlnFile& aln, const GenomicLocus& pos
);

// default pileup function: yields primary, non-dup,
// non-qcfail, non-supplementary reads. Consumes a
// PileupContext, so usable with any derived context which
// does not need extra filtering.
extern "C" int pileup_func (void* data, bam1_t* b);

// SUPERCEDING PreparedPileup

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

  explicit PileupContext (const AlnFile& aln)
      : br_fh{aln.fh}, br_fhIdx{aln.idx}
  {
  }

  PileupContext (const PileupContext&) = delete;
  PileupContext& operator= (const PileupContext&) = delete;

  // virtual: ForwardPileupIterator destroys the context
  // through this type, so a derived context must get its
  // own destructor run.
  virtual ~PileupContext()
  {
    if (it != nullptr) {
      hts_itr_destroy (it);
    }
  }
};

struct ForwardPileupIterator {
  static constexpr hts_pos_t SEEK_GAP = 10000;

  struct {
    // View into plpIt's buffer; invalidated by the next
    // bam_plp64_auto or bam_plp_reset on plpIt.
    const bam_pileup1_t* plpArr = nullptr;
    size_t nReads = 0;
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

  // constructor
  friend ForwardPileupIterator init_pileup_iterator (
      PileupContext* snk_ctx, bam_plp_auto_f fn
  );

 private:
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

enum class PileupAdvErr : uint8_t {
  invalidLocus,
  locusBehindItr,
  samItrFailed,
  pileupEngineFailed
};
enum class PileupAdvResult : uint8_t {
  success,
  exhausted,
  noCoverage
};
// Move cursor to `to`.
// Separate to accessing the result as the
// outcomes are quite rich.
using AdvResultOrErr =
    std::expected<PileupAdvResult, PileupAdvErr>;
[[nodiscard]] AdvResultOrErr try_advance_pileup (
    ForwardPileupIterator& pc, const GenomicLocus& to
);

// get non-owning view into pileup at current position.
// `at` enforces that you are retrieving the position
// you expected.
// empty span == no coverage
// nullopt == cursor not at `at`
using PileupView = std::span<const bam_pileup1_t>;
std::optional<PileupView> read (
    const ForwardPileupIterator& pc, const GenomicLocus& at
);
