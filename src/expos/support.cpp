#include "support.hpp"

#include <format>
#include <string>

bool
check_aln
(const bam_pileup1_t* p)
{
  // at refskip/is_del bases
  // qpos will be -1.
  // These cannot support a variant
  // since all variants are anchored by a real base
  return p->qpos >= 0;
}

int
get_var_indel_len
(const bcf1_t* m) {
  std::string ref = m->d.allele[0];
  std::string alt = m->d.allele[1];

  // negative for del, pos for ins
  return
    static_cast<int> (alt.length())
    - static_cast<int> (ref.length());

}


bool
eval_support_snp
(const bcf1_t* m, const bam_pileup1_t* p) {
  if (!check_aln (p)) {
    return false;
  }
  
  const auto alt = m->d.allele[1];
  const auto qpos = static_cast<size_t> (p->qpos);
  const auto qbase = seq_nt16_str[bam_seqi (bam_get_seq (p->b), qpos)];
  const auto abase = alt[0];
  if (qbase != abase) {
      return false;
  }
  return true;
}

bool
eval_support_mnp
(const bcf1_t* m, const bam_pileup1_t* p)
{
  if (!check_aln (p)) {
    return false;
  }

  std::string alt = m->d.allele[1];
  const auto qpos = static_cast<size_t> (p->qpos);

  std::string read_bases;
  for (size_t i = 0; i < alt.length(); ++i) {
    read_bases.push_back (
        seq_nt16_str[bam_seqi (bam_get_seq (p->b), qpos + i)]
    );
  }

  return (alt == read_bases);
}


bool
eval_support_del
(const bcf1_t* m, const bam_pileup1_t* p)
{
  return check_aln (p) && (get_var_indel_len (m) == p->indel);
}

bool
eval_support_ins
(const bcf1_t* m, const bam_pileup1_t* p)
{
  if (!check_aln (p)) {
    return false;
  }
  kstring_t ins;
  ks_initialize (&ins);
  auto bam_ins_len = bam_plp_insertion (p, &ins, NULL);
  if (bam_ins_len < 0) {
    throw std::runtime_error (
      std::format ("failed to retrieve insertion bases for read {}", bam_get_qname (p->b))
    );
  }
  if (bam_ins_len != get_var_indel_len (m)) {
      ks_free (&ins);
      return false;
  }

  // on the basis that insertions seem to go
  // anchor-insertion-rest when ambiguous. And that this gets the most success
  // discussed with htslib team, but really a question for bcftools team
  // and VCF format spec group.
  std::string alt = m->d.allele[1];
  const auto var_ins = alt.substr (1, static_cast<size_t> (bam_ins_len));
  auto ins_match =
    strcmp (ks_c_str (&ins), var_ins.c_str()) == 0;  // check bases match

  ks_free (&ins);

  return (ins_match);
}

