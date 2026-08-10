#include "hts_types.hpp"

#include <fmt/format.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <plog/Log.h>

#include <optional>

#include "shared/err.hpp"

AlnOrErr load_aln (const char* fn)
{
  AlnFile aln;
  aln.fh = hts_open (fn, "r");
  if (aln.fh == nullptr) {
    return std::unexpected (make_err (
        fmt::format ("Could not open alignment file at {}", fn)
    ));
  }
  aln.hdr = sam_hdr_read (aln.fh);
  if (aln.hdr == nullptr) {
    return std::unexpected (make_err (
        "Could not read header "
        "from alignment file"
    ));
  }
  aln.idx = sam_index_load (aln.fh, fn);
  if (aln.idx == nullptr) {
    return std::unexpected (make_err (
        fmt::format ("Could not load index for {}", fn)
    ));
  }

  return aln;
}

VcfOrErr load_vcf (const char* fn)
{
  VcfFile vcf;
  vcf.fh = hts_open (fn, "r");
  if (vcf.fh == nullptr) {
    return std::unexpected (make_err (
        fmt::format ("Could not open VCF file at {}", fn)
    ));
  }
  vcf.hdr = bcf_hdr_read (vcf.fh);
  if (vcf.hdr == nullptr) {
    return std::unexpected (
        make_err ("Could not read header from VCF file")
    );
  }

  return vcf;
}

FastaOrErr load_fasta (const char* fn)
{
  FastaFile ff;

  ff.fai = fai_load3_format (
      fn, NULL, NULL, 0, fai_format_options::FAI_FASTA
  );

  if (ff.fai == nullptr) {
    return std::unexpected (make_err (
        fmt::format ("Could not open fasta file at {}", fn)
    ));
  }

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

extern "C" {
int pileup_func (void* data, bam1_t* b)
{
  constexpr uint16_t EXCLUDE_BITS = 3840;
  const PileupContext* d = (PileupContext*)(data);
  int ret = -1;

  // find the next good read
  while (true) {
    ret = sam_itr_next (d->br_fh, d->it, b);
    if (ret < 0) {
      break;     // EOF/err
    }
    if ((b->core.flag & EXCLUDE_BITS) == 0) {
      break;  // primary, non-dup, non-qcfail, non-supplementary read
    }
  }
  return ret;
}
}

PileupOrErr prepare_pileup (
    const AlnFile& aln, const GenomicLocus& pos
)
{
  PreparedPileup out;

  auto* alnIter =
      sam_itr_queryi (aln.idx, pos.tid, pos.pos, pos.pos + 1);
  if (alnIter == NULL) {
    return std::unexpected{make_err (
        "sam_itr_queryi: failed "
        "to create iterator"
    )};
  }

  PileupContext pc{aln};
  pc.it = alnIter;
  auto* plp = bam_plp_init (pileup_func, &pc);
  if (plp == NULL) {
    return std::unexpected{make_err (
        "bam_plp_init: failed to initialise "
        "pileup engine"
    )};
  }
  out.plpBacking = plp;

  int64_t plpPos = -1;
  int plpTid = -1;
  int nPlp = -1;
  const bam_pileup1_t* plpArr;
  PLOGD << "Iterating pileup";
  while ((plpArr = bam_plp64_auto (
              out.plpBacking, &plpTid, &plpPos, &nPlp
          )) != 0) {
    if (nPlp < 0 || plpTid < 0 || plpPos < 0) {
      return std::unexpected{
          make_err ("bam_plp64_auto: pileup failed")
      };
    }
    if (plpPos < pos.pos) {
      continue;  // doesn't cover variant
    }
    PLOGD << "Position found";
    out.plpArr = plpArr;
    out.nPlp = static_cast<size_t> (nPlp);
    return out;
  }
  PLOGD << "Position not covered by alignment file";
  return {};
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

AdvResultOrErr try_advance_pileup (
    ForwardPileupIterator& pc, const GenomicLocus& to
)
{
  auto& ctx = pc.ctx;
  auto& curs = pc.curs;
  auto& curLoc = curs.locus;

  if (to.pos < 0 || to.tid < 0) {
    return std::unexpected (PileupAdvErr::invalidLocus);
  }
  if (pc.lastRequest.tid == to.tid &&
      pc.lastRequest.pos > to.pos) {
    return std::unexpected (PileupAdvErr::locusBehindItr);
  }
  bool unInit =
      (curLoc.tid < 0 &&
       pc.lastRequest.tid < 0);  // first advance request
  pc.lastRequest = to;
  bool sameTid = (curLoc.tid == to.tid);

  if (!unInit && sameTid) {
    if (curLoc.pos == to.pos) {
      return PileupAdvResult::success;
    }
    if (curLoc.pos > to.pos) {
      // Assuming correct use (forward iteration only)
      // then the position MUST have had no coverage
      // to get to this state.
      return PileupAdvResult::noCoverage;
    }
  }

  bool needSeek =
      unInit || !sameTid ||
      (to.pos - curLoc.pos) > ForwardPileupIterator::SEEK_GAP;

  if (needSeek) {
    auto* alnIter = sam_itr_queryi (
        ctx->br_fhIdx, to.tid, to.pos, HTS_POS_MAX
    );
    if (alnIter == NULL) {
      return std::unexpected (PileupAdvErr::samItrFailed);
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
    if (nPlp < 0 || plpTid < 0 || plpPos < 0) {
      return std::unexpected (PileupAdvErr::pileupEngineFailed);
    }
    if (plpPos == to.pos) {
      curs.plpArr = plpArr;
      curs.nReads = static_cast<size_t> (nPlp);
      curLoc = to;
      return PileupAdvResult::success;
    }
    if (plpPos > to.pos) {
      curs.plpArr = plpArr;
      curs.nReads = static_cast<size_t> (nPlp);
      // parked at
      curLoc = {plpTid, plpPos};
      return PileupAdvResult::noCoverage;
    }
  }
  // iterator exhausted before `to` found.
  curs.plpArr = nullptr;
  curs.nReads = 0;
  // parked at
  curLoc = {plpTid, plpPos};
  return PileupAdvResult::exhausted;
}

std::optional<PileupView> read (
    const ForwardPileupIterator& pc, const GenomicLocus& at
)
{
  if (pc.curs.locus != at) {
    return std::nullopt;
  }
  return std::span (pc.curs.plpArr, pc.curs.nReads);
}
