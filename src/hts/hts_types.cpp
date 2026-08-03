#include "hts_types.hpp"

#include <fmt/format.h>
#include <htslib/faidx.h>
#include <plog/Log.h>

AlnOrErr load_aln (const char* fn)
{
  AlnFile aln;
  aln.o_fh = hts_open (fn, "r");
  if (aln.o_fh == nullptr) {
    return std::unexpected (
       make_err (fmt::format ("Could not open alignment file at {}", fn))
    );
  }
  aln.o_hdr = sam_hdr_read (aln.o_fh);
  if (aln.o_hdr == nullptr) {
    return std::unexpected (make_err (
        "Could not read header "
        "from alignment file"
    ));
  }
  aln.o_idx = sam_index_load (aln.o_fh, fn);
  if (aln.o_idx == nullptr) {
    return std::unexpected (make_err (
        fmt::format ("Could not load index for {}", fn)
    ));
  }

  return aln;
}

VcfOrErr load_vcf (const char* fn)
{
  VcfFile vcf;
  vcf.o_fh = hts_open (fn, "r");
  if (vcf.o_fh == nullptr) {
    return std::unexpected (
        make_err (fmt::format ("Could not open VCF file at {}", fn))
    );
  }
  vcf.o_hdr = bcf_hdr_read (vcf.o_fh);
  if (vcf.o_hdr == nullptr) {
    return std::unexpected (
        make_err ("Could not read header from VCF file")
    );
  }

  return vcf;
}

FastaOrErr load_fasta (const char* fn)
{
  FastaFile ff;

  ff.o_fai = fai_load3_format (
      fn, NULL, NULL, 0, fai_format_options::FAI_FASTA
  );

  if (ff.o_fai == nullptr) {
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
  auto* p_fetch = faidx_fetch_seq64 (
      ff, contigName.data(), regStart, regEnd - 1, &rc
  );
  if (p_fetch == NULL) {
    std::string msg{"Could not retrieve region from fasta; "};
    if (rc == -2) {
      msg += fmt::format ("contig {} not found", contigName);
    }
    else {
      msg += "unspecified htslib error";
    }
    return std::unexpected (
        make_err (msg)
    );
  }

  std::string out{p_fetch, static_cast<size_t> (rc)};

  free (p_fetch);

  return out;
}


struct PileupCapture {
  htsFile* fh = nullptr;
  hts_itr_t* it = nullptr;

  ~PileupCapture () {
    if (it != nullptr) {
      hts_itr_destroy(it);
    }
  }
};
extern "C" {
int pileup_func (void* data, bam1_t* b)
{
  constexpr uint16_t k_excludeBits = 3840;
  const PileupCapture* d = (PileupCapture*)(data);
  int ret = -1;

  // find the next good read
  while (true) {
    ret = sam_itr_next (d->fh, d->it, b);
    if (ret < 0) {
      break;     // EOF/err
    }
    if ((b->core.flag & k_excludeBits) == 0) {
      break;  // primary, non-dup, non-qcfail, non-supplementary read
    }
  }
  return ret;
}
}

PileupOrErr prepare_pileup (
    const AlnFile& aln, const PileupPosition& pos
)
{
  PLOGD << "Begin prepare_pileup";
  PreparedPileup out;

  PLOGD << "Initalising sam_itr_queryi";
  auto* alnIter =
      sam_itr_queryi (aln.o_idx, pos.tid, pos.pos, pos.pos + 1);
  if (alnIter == NULL) {
    return std::unexpected{make_err (
        "sam_itr_queryi: failed "
        "to create iterator"
    )};
  }

  PLOGD << "Initialising bam_plp_t";
  PileupCapture pc{aln.o_fh, alnIter};
  auto* plp = bam_plp_init (pileup_func, &pc);
  if (plp == NULL) {
    return std::unexpected{make_err (
        "bam_plp_init: failed to initialise "
        "pileup engine"
    )};
  }
  out.o_plp = plp;

  int64_t plpPos = -1;
  int plpTid = -1;
  int nPlp = -1;
  const bam_pileup1_t* plpArr;
  PLOGD << "Iterating pileup";
  while ((plpArr = bam_plp64_auto (
              out.o_plp, &plpTid, &plpPos, &nPlp
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
