#include "hts_types.hpp"

#include <fmt/format.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <plog/Log.h>

#include <optional>
#include <utility>

#include "shared/err.hpp"

AlnFile::AlnFile (AlnFile&& o) noexcept
    : fh{o.fh}, hdr{o.hdr}, idx{o.idx}, path{std::move (o.path)}
{
  o.fh = NULL;
  o.hdr = NULL;
  o.idx = NULL;
}

AlnFile& AlnFile::operator= (AlnFile&& o) noexcept
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
    path = std::move (o.path);
    o.fh = NULL;
    o.hdr = NULL;
    o.idx = NULL;
  }
  return *this;
}

AlnFile::~AlnFile() noexcept
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

AlnOrErr load_aln (const std::string& path)
{
  AlnFile aln;
  aln.fh = hts_open (path.c_str(), "r");
  if (aln.fh == nullptr) {
    return std::unexpected (make_err (
        fmt::format ("Could not open alignment file at {}", path)
    ));
  }
  aln.hdr = sam_hdr_read (aln.fh);
  if (aln.hdr == nullptr) {
    return std::unexpected (make_err (
        "Could not read header "
        "from alignment file"
    ));
  }
  aln.idx = sam_index_load (aln.fh, path.c_str());
  if (aln.idx == nullptr) {
    return std::unexpected (make_err (
        fmt::format ("Could not load index for {}", path)
    ));
  }

  aln.path = path;

  return aln;
}

VcfFile::VcfFile (VcfFile&& o) noexcept
    : fh{o.fh}, hdr{o.hdr}, path{std::move (o.path)}
{
  o.fh = NULL;
  o.hdr = NULL;
}

VcfFile& VcfFile::operator= (VcfFile&& o) noexcept
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
    path = std::move (o.path);
    o.fh = NULL;
    o.hdr = NULL;
  }
  return *this;
}

VcfFile::~VcfFile() noexcept
{
  if (hdr != nullptr) {
    bcf_hdr_destroy (hdr);
  }
  if (fh != nullptr) {
    hts_close (fh);
  }
}

VcfOrErr load_vcf (const std::string& path)
{
  VcfFile vcf;
  vcf.fh = hts_open (path.c_str(), "r");
  if (vcf.fh == nullptr) {
    return std::unexpected (make_err (
        fmt::format ("Could not open VCF file at {}", path)
    ));
  }
  vcf.hdr = bcf_hdr_read (vcf.fh);
  if (vcf.hdr == nullptr) {
    return std::unexpected (
        make_err ("Could not read header from VCF file")
    );
  }

  vcf.path = path;

  return vcf;
}

FastaFile::FastaFile (FastaFile&& o) noexcept
    : fai{o.fai}, path{std::move (o.path)}
{
  o.fai = nullptr;
}

FastaFile& FastaFile::operator= (FastaFile&& o) noexcept
{
  if (this != &o) {
    if (fai != nullptr) {
      fai_destroy (fai);
    }
    fai = o.fai;
    path = std::move (o.path);
    o.fai = nullptr;
  }
  return *this;
}

FastaFile::~FastaFile()
{
  if (fai != nullptr) {
    fai_destroy (fai);
  }
}

FastaOrErr load_fasta (const std::string& path)
{
  FastaFile ff;

  ff.fai = fai_load3_format (
      path.c_str(), NULL, NULL, 0, fai_format_options::FAI_FASTA
  );

  if (ff.fai == nullptr) {
    return std::unexpected (make_err (
        fmt::format ("Could not open fasta file at {}", path)
    ));
  }

  ff.path = path;

  return ff;
}

RefSliceOrErr fetch_region (
    const FastaFile& ff, const std::string_view contigName,
    hts_pos_t regStart, hts_pos_t regEnd
)
{
  hts_pos_t rc;
  auto* seq = faidx_fetch_seq64 (
      ff, contigName.data(), regStart, regEnd - 1, &rc
  );
  if (seq == NULL) {
    std::string msg{"Could not retrieve region from fasta; "};
    if (rc == -2) {
      msg += fmt::format ("contig {} not found", contigName);
    }
    else {
      msg += "unspecified htslib error";
    }
    return std::unexpected (make_err (msg));
  }

  std::string out{seq, static_cast<size_t> (rc)};

  free (seq);

  return out;
}

PileupContext::PileupContext (const AlnFile& aln)
    : br_fh{aln.fh}, br_fhIdx{aln.idx}
{
}

PileupContext::~PileupContext()
{
  if (it != nullptr) {
    hts_itr_destroy (it);
  }
}

ForwardPileupIterator::~ForwardPileupIterator()
{
  if (plpIt != nullptr) {
    bam_plp_destroy (plpIt);
  }
  delete ctx;
}

ForwardPileupIterator::ForwardPileupIterator (
    ForwardPileupIterator&& o
) noexcept
    : curs (o.curs),
      lastRequest (o.lastRequest),
      ctx (o.ctx),
      plpIt (o.plpIt)
{
  o.ctx = nullptr;
  o.plpIt = nullptr;
  o.curs.plpArr = nullptr;
  o.curs.nReads = 0;
  o.curs.locus = {-1, -1};
  o.lastRequest = {-1, -1};
}

ForwardPileupIterator& ForwardPileupIterator::operator= (
    ForwardPileupIterator&& o
) noexcept
{
  if (this != &o) {
    if (plpIt != nullptr) {
      bam_plp_destroy (plpIt);
    }
    delete ctx;
    curs = o.curs;
    lastRequest = o.lastRequest;
    ctx = o.ctx;
    plpIt = o.plpIt;
    o.ctx = nullptr;
    o.plpIt = nullptr;
    o.curs.plpArr = nullptr;
    o.curs.nReads = 0;
    o.curs.locus = {-1, -1};
    o.lastRequest = {-1, -1};
  }
  return *this;
}

ForwardPileupIterator init_pileup_iterator (
    PileupContext* snk_ctx, bam_plp_auto_f fn
)
{
  ForwardPileupIterator out;

  out.ctx = snk_ctx;
  out.plpIt = bam_plp_init (fn, out.ctx);  // can't fail

  return out;
}

AdvResOrErr try_advance_pileup (
    ForwardPileupIterator& pc, const GenomicLocus& to
)
{
  auto& ctx = pc.ctx;
  auto& curs = pc.curs;
  auto& curLoc = curs.locus;

  if (to.pos < 0 || to.tid < 0 || to.pos == HTS_POS_MAX) {
    return std::unexpected (PileupAdvErrCode::invalidLocus);
  }
  if (pc.lastRequest.tid == to.tid &&
      pc.lastRequest.pos > to.pos) {
    return std::unexpected (PileupAdvErrCode::locusBehindItr);
  }
  const bool unInit = (curLoc.tid < 0);  // first advance request
  pc.lastRequest = to;
  const bool sameTid = (curLoc.tid == to.tid);

  if (!unInit && sameTid) {
    if (curLoc.pos == to.pos) {
      // Can't be exhausted, since exhausted
      // sentinel is not a valid input for `to`.
      // Can't be noCoverage, since cursor
      // never parks on a 0 coverage location
      // (htslib never returns non-null when _n_plp = 0).
      return PileupAdvResCode::success;
    }
    if (curLoc.pos > to.pos) {
      // Assuming correct use (forward iteration only)
      // then the position MUST have had no coverage
      // to get to this state.
      return PileupAdvResCode::noCoverage;
    }
  }

  const bool needSeek =
      unInit || !sameTid ||
      (to.pos - curLoc.pos) > ForwardPileupIterator::SEEK_GAP;

  if (needSeek) {
    auto* alnIter = sam_itr_queryi (
        ctx->br_fhIdx, to.tid, to.pos, HTS_POS_MAX
    );
    if (alnIter == NULL) {
      return std::unexpected (PileupAdvErrCode::samItrFailed);
    }
    if (ctx->it != nullptr) {
      hts_itr_destroy (ctx->it);
    }
    ctx->it = alnIter;
    bam_plp_reset (pc.plpIt);
  }

  int64_t plpPos = -1;
  int plpTid = -1;
  int nPlp = -1;
  const bam_pileup1_t* plpArr;
  while ((plpArr = bam_plp64_auto (
              pc.plpIt, &plpTid, &plpPos, &nPlp
          )) != 0) {
    if (plpPos == to.pos) {
      curs.plpArr = plpArr;
      curs.nReads = static_cast<size_t> (nPlp);
      curLoc = to;
      return PileupAdvResCode::success;
    }
    if (plpPos > to.pos) {
      curs.plpArr = plpArr;
      curs.nReads = static_cast<size_t> (nPlp);
      // parked at
      curLoc = {plpTid, plpPos};
      return PileupAdvResCode::noCoverage;
    }
  }
  // nPlp alone carries error info
  if (nPlp < 0) {
    return std::unexpected (
        PileupAdvErrCode::pileupEngineFailed
    );
  }
  // null return + _n_plp = 0: EOF,
  // iterator exhausted before `to` found.
  curs.plpArr = nullptr;
  curs.nReads = 0;
  // parked at (EOF sentinel)
  curLoc = {to.tid, HTS_POS_MAX};
  return PileupAdvResCode::exhausted;
}

PileupViewOrNone read (
    const ForwardPileupIterator& pc,
    const GenomicLocus& expectedLocus
)
{
  if (pc.curs.locus != expectedLocus) {
    return std::unexpected (std::nullopt);
  }
  return std::span (pc.curs.plpArr, pc.curs.nReads);
}
