#include "vcf_record.hpp"

#include <fmt/format.h>
#include <htslib/vcf.h>

#include <string>
#include <string_view>

static bool check_aln (const bam_pileup1_t& p)
{
  // at refskip/is_del bases
  // qpos will be -1.
  // These cannot support a variant
  // since all variants are anchored by a real base
  return (!p.is_del) && (!p.is_refskip) && p.qpos >= 0;
}

// The eval_support_* functions read REF/ALT straight off the record. They
// are only ever reached through VcfRec::eval_read_support, which
// set_supporting_read_evaluator only sets on biallelic, typeable records —
// so d.allele[1] is always present here.
static std::string_view ref_allele (const VcfRec& r)
{
  return r.ptr->d.allele[0];
}

static std::string_view alt_allele (const VcfRec& r)
{
  return r.ptr->d.allele[1];
}

// ALT length minus REF length: negative for a deletion, positive for an
// insertion, matching the sign convention of bam_pileup1_t::indel.
static int indel_len (const VcfRec& r)
{
  return static_cast<int> (alt_allele (r).length()) -
         static_cast<int> (ref_allele (r).length());
}

bool eval_support_snp (const VcfRec& r, const bam_pileup1_t& p)
{
  if (!check_aln (p)) {
    return false;
  }

  const auto qpos = static_cast<size_t> (p.qpos);
  const auto qbase =
      seq_nt16_str[bam_seqi (bam_get_seq (p.b), qpos)];
  const auto abase = alt_allele (r)[0];
  return qbase == abase;
}

bool eval_support_mnp (const VcfRec& r, const bam_pileup1_t& p)
{
  if (!check_aln (p)) {
    return false;
  }

  const auto qpos = static_cast<size_t> (p.qpos);
  const auto alt = alt_allele (r);

  std::string read_bases;
  for (size_t i = 0; i < alt.length(); ++i) {
    read_bases.push_back (
        seq_nt16_str[bam_seqi (bam_get_seq (p.b), qpos + i)]
    );
  }

  return (alt == read_bases);
}

bool eval_support_del (const VcfRec& r, const bam_pileup1_t& p)
{
  return check_aln (p) && (indel_len (r) == p.indel);
}

bool eval_support_ins (const VcfRec& r, const bam_pileup1_t& p)
{
  if (!check_aln (p)) {
    return false;
  }
  kstring_t ins;
  ks_initialize (&ins);
  auto bam_ins_len = bam_plp_insertion (&p, &ins, NULL);
  if (bam_ins_len < 0) {
    throw std::runtime_error (
        fmt::format (
            "failed to retrieve insertion bases for read {}",
            bam_get_qname (p.b)
        )
    );
  }
  if (bam_ins_len != indel_len (r)) {
    ks_free (&ins);
    return false;
  }

  // on the basis that insertions seem to go
  // anchor-insertion-rest when ambiguous. And that this gets the most success
  // discussed with htslib team, but really a question for bcftools team
  // and VCF format spec group.
  // bases beyond the anchor base shared with REF
  const std::string_view ins_bases{ks_c_str (&ins), ins.l};
  auto ins_match = ins_bases == alt_allele (r).substr (1);

  ks_free (&ins);

  return (ins_match);
}

std::string stringify_rec (const VcfRec& r)
{
  return fmt::format (
      "{}:{} (id: {})", r.contig, r.ptr->pos + 1, r.ptr->d.id
  );
}

std::expected<GenomicLocus, Rec2LocusErr> record_to_aln_locus (
    const VcfRec& r, const AlnFile& aln
)
{
  // sam_hdr_name2tid: -1 contig not present in alignment, -2 parse failure.
  const auto alnTid =
      sam_hdr_name2tid (aln.hdr, r.contig.data());
  if (alnTid < 0) {
    return std::unexpected (Rec2LocusErr{alnTid});
  }
  return GenomicLocus{alnTid, r.ptr->pos};
}

// True if the record is one of the four supported mutation types, in which
// case r.eval_read_support is set to the matching predicate. Otherwise the
// record is complex or untypeable and the predicate is cleared.
static bool set_supporting_read_evaluator (VcfRec& r)
{
  auto mtype = bcf_has_variant_type (
      r.ptr, 1, VCF_DEL | VCF_INS | VCF_SNP | VCF_MNP
  );
  switch (mtype) {
    // one and only one of
    case (VCF_SNP):
      r.eval_read_support = eval_support_snp;
      return true;
    case (VCF_MNP):
      r.eval_read_support = eval_support_mnp;
      return true;
    case (VCF_DEL):
      r.eval_read_support = eval_support_del;
      return true;
    case (VCF_INS):
      r.eval_read_support = eval_support_ins;
      return true;
    default:
      // unsupported/complex mutation
      r.eval_read_support = nullptr;
      return false;
  }
}

bool record_is_filtered (const VcfRec& r)
{
  const auto nFlt = r.ptr->d.n_flt;
  if (nFlt == 0) {
    return false;  // '.' — no filtering applied
  }
  if (nFlt == 1 && r.ptr->d.flt[0] == 0) {
    return false;  // PASS (always FILTER id 0)
  }
  return true;  // one or more failing filters (PASS cannot co-occur)
}

VoidOrErr check_integrity (VcfRec& r)
{
  const auto* br_rec = r.ptr;

  // Resolve the contig before anything that reports the record.
  // stringify_rec reads r.contig, and one VcfRec is reused for every
  // record in the run, so leaving it until later would make the messages
  // below name the *previous* record's contig.
  const auto* contig = bcf_hdr_id2name (r.br_hdr, r.ptr->rid);
  r.contig = (contig != nullptr) ? std::string_view{contig}
                                 : std::string_view{};

  if (br_rec->errcode != 0) {
    std::string errMsgOut{"Error while reading VCF record "};
    errMsgOut += stringify_rec (r);
    constexpr int ERR_BUF_SIZE = 300;
    char errBuf[ERR_BUF_SIZE];
    const auto* bcfErrMsg =
        bcf_strerror (br_rec->errcode, errBuf, ERR_BUF_SIZE);
    if (bcfErrMsg == nullptr) {
      errMsgOut += ", failed to recover error message";
    }
    else {
      errMsgOut += ", reporting: ";
      errMsgOut += bcfErrMsg;
    }
    return std::unexpected (make_err (errMsgOut));
  }

  if (contig == nullptr) {
    return std::unexpected (make_err (
        fmt::format (
            "Error while reading VCF record {}, could not find "
            "contig in header. VCF likely corrupt.",
            stringify_rec (r)
        )
    ));
  }

  return {};
}

std::expected<void, RecordSkipReason> validate_for_analysis (
    VcfRec& r
)
{
  if (r.ptr->n_allele != 2) {
    return std::unexpected (RecordSkipReason::notBiallelic);
  }
  if (!set_supporting_read_evaluator (r)) {
    return std::unexpected (RecordSkipReason::complex);
  }
  return {};
}
