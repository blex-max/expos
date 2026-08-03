#include "vcf_record.hpp"

#include <fmt/format.h>

#include <format>
#include <string>

static bool check_aln (const bam_pileup1_t& p)
{
  // at refskip/is_del bases
  // qpos will be -1.
  // These cannot support a variant
  // since all variants are anchored by a real base
  return (!p.is_del) && (!p.is_refskip) && p.qpos >= 0;
}

bool eval_support_snp (const VcfRec& r, const bam_pileup1_t& p)
{
  if (!check_aln (p)) {
    return false;
  }

  const auto qpos = static_cast<size_t> (p.qpos);
  const auto qbase =
      seq_nt16_str[bam_seqi (bam_get_seq (p.b), qpos)];
  const auto abase = r.altAllele[0];
  return qbase == abase;
}

bool eval_support_mnp (const VcfRec& r, const bam_pileup1_t& p)
{
  if (!check_aln (p)) {
    return false;
  }

  const auto qpos = static_cast<size_t> (p.qpos);

  std::string read_bases;
  for (size_t i = 0; i < r.altAllele.length(); ++i) {
    read_bases.push_back (
        seq_nt16_str[bam_seqi (bam_get_seq (p.b), qpos + i)]
    );
  }

  return (r.altAllele == read_bases);
}

bool eval_support_del (const VcfRec& r, const bam_pileup1_t& p)
{
  return check_aln (p) && (r.indelLen == p.indel);
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
        std::format (
            "failed to retrieve insertion bases for read {}",
            bam_get_qname (p.b)
        )
    );
  }
  if (bam_ins_len != r.indelLen) {
    ks_free (&ins);
    return false;
  }

  // on the basis that insertions seem to go
  // anchor-insertion-rest when ambiguous. And that this gets the most success
  // discussed with htslib team, but really a question for bcftools team
  // and VCF format spec group.
  const std::string_view ins_bases{ks_c_str (&ins), ins.l};
  auto ins_match = ins_bases == r.insertedBases;

  ks_free (&ins);

  return (ins_match);
}

std::string stringify_rec (const VcfRec& r)
{
  return fmt::format (
      "{}:{} (id: {})", bcf_hdr_id2name (r.hdr, r.ptr->rid),
      r.ptr->pos + 1, r.ptr->d.id
  );
}

RecordClass type_record (VcfRec& r)
{
  r.refAllele = r.ptr->d.allele[0];
  r.altAllele = r.ptr->d.allele[1];
  r.indelLen = static_cast<int> (r.altAllele.length()) -
               static_cast<int> (r.refAllele.length());
  // valid for INS records only; harmless otherwise since unused
  r.insertedBases =
      r.altAllele.substr (1);  // all bases after anchor

  auto mtype = bcf_has_variant_type (
      r.ptr, 1, VCF_DEL | VCF_INS | VCF_SNP | VCF_MNP
  );
  switch (mtype) {
    // one and only one of
    case (VCF_SNP):
      r.eval_read_support = eval_support_snp;
      return RecordClass::Analysable;
    case (VCF_MNP):
      r.eval_read_support = eval_support_mnp;
      return RecordClass::Analysable;
    case (VCF_DEL):
      r.eval_read_support = eval_support_del;
      return RecordClass::Analysable;
    case (VCF_INS):
      r.eval_read_support = eval_support_ins;
      return RecordClass::Analysable;
    default:
      // unsupported/complex mutation, or could not be typed
      return RecordClass::Unanalysable;
  }
};

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

VoidOrErr check_record_integrity (const VcfRec& r)
{
  const auto* uo_rec = r.ptr;

  if (uo_rec->errcode != 0) {
    std::string errMsgOut{"Error while reading VCF record "};
    errMsgOut += stringify_rec (r);
    constexpr int k_errBufSize = 300;
    char errBuf[k_errBufSize];
    const auto* bcfErrMsg =
        bcf_strerror (uo_rec->errcode, errBuf, k_errBufSize);
    if (bcfErrMsg == nullptr) {
      errMsgOut += ", failed to recover error message";
    }
    else {
      errMsgOut += ", reporting: ";
      errMsgOut += bcfErrMsg;
    }
    return std::unexpected (make_err (errMsgOut));
  }
  return {};
}
